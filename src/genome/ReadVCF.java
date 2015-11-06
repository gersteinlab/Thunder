package genome;

import java.sql.PreparedStatement;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import annotation.Chromosome;
import database.DBConnect;
import database.DBConnect_SQLite;


public class ReadVCF {

	public static void createTable(Statement st, String tableName) throws SQLException{

		st.execute("create table "+tableName+"(" +
				"chromosome varchar(30) not null, " +
				"position int not null, " +
				"variantID varchar(30), " +
				"reference varchar(10) not null, " +
				"alternative varchar(10) not null, " +
				"isPhased boolean FALSE, " +
				"genotype varchar(5), " +
				"altAlleleCount int null,"+
				"altAlleleFrequency double null,"+
				"wgs_ReadDepth double null, " +
				"wgs_Quality double null, " +
				"chip_haplotypeRatio double null, " +
				//				"chip_logRRatio double null)");
				"chip_logRRatio double null, " +
				"primary key (chromosome, position))");
	}




	public static void ReadVCFToDB(String vcfFile, String dbPath, String sampleID, boolean forceRefresh) throws Exception{

		String tableName = "genomicVariants";
		DBConnect db = null;
		Statement st = null;
		PreparedStatement ps = null;

		try{
			db = new DBConnect_SQLite(dbPath);
			st = db.createStatement();
			
			if(forceRefresh)
				st.execute("drop table if exists "+tableName);
			
			if(!db.containsTable(tableName)){
				createTable(st, tableName);
				
				ps = db.getPreparedStatement("INSERT INTO "+tableName+" VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
				db.setAutoCommit(false);
				
				VCFparser parser = new VCFparser(vcfFile,sampleID);
				int n_ev = 0,var_nucs = 0;
				while (parser.hasMoreInput()) {
					Variant var = parser.parseLine();
					if (var == null) continue;
					if (var.getMaternal() == 0 && var.getPaternal() == 0) continue;

					/*
				System.out.println("chrom: "+Chromosome.getChromID(var.getChromosome()));
				System.out.println("index: "+var.getPosition());
				System.out.println("ID:    "+var.getID());
				System.out.println("Ref:   "+var.getReference());
				System.out.println("Alt:   "+var.getAlternative());
				System.out.println("Phased:"+var.isPhased());
				System.out.println("GT:    "+var.getMaternal()+"/"+var.getPaternal());
				System.out.println("AF:    "+((var.getMaternal()+var.getPaternal())/2.0));
				System.out.println();*/

					ps.setString(1, Chromosome.getChromID(var.getChromosome()));
					ps.setInt(2, var.getPosition());
					ps.setString(3, var.getID());
					ps.setString(4, var.getReference());
					ps.setString(5, var.getAlternative());
					ps.setBoolean(6, var.isPhased());
					ps.setString(7, var.getMaternal()+"/"+var.getPaternal());
					ps.setInt(8, 2);
					ps.setDouble(9, ((var.getMaternal()+var.getPaternal())/2.0));
					ps.setDouble(10, -1.0);
					ps.setDouble(11, -1.0);

					String tmp;
					if((tmp=var.getAttribute("HR")) != null)
						ps.setDouble(12, Double.valueOf(tmp).doubleValue());
					else
						ps.setDouble(12, 0.0);
					if((tmp=var.getAttribute("LR")) != null)
						ps.setDouble(13, Double.valueOf(tmp).doubleValue());
					else
						ps.setDouble(13, 0.0);

					ps.addBatch();

					n_ev++;
					var_nucs += var.variantBases();
				}

				ps.executeBatch();
				db.commit();
				db.setAutoCommit(true);
				ps.close();
				
				System.out.println(vcfFile + ": " + n_ev + " variants, " +
						var_nucs + " variant bases");
			}
			else{
				/* table already exists, do nothing */
				ResultSet rs = st.executeQuery("SELECT COUNT(*) FROM "+tableName);
				System.out.println("Database already contains information for "+rs.getInt(1)+" variants. To overwrite this data use the '-"+Thunder.OPT_CHOICE_DB_FORCE_REFRESH+"' option.");
			}

		}catch(Exception e){
			e.printStackTrace();
		}finally{
			st.close();
			db.closeConnection();
		}
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("vcf").hasArg().withDescription("path to the file containing the GENCODE annotation data in GTF format").create(Thunder.OPT_PATH_INPUT));
		options.addOption(OptionBuilder.withArgName("databasePath").hasArg().withDescription("create or add to the database at this path\n(e.g. /path/to/db)").create(Thunder.OPT_PATH_DB_SAMPLE));
		//		options.addOption(OptionBuilder.withArgName("tableName").hasArg().withDescription("name of the table to which to add the annotation data").create(Thunder.OPT_DB_TABLE_NAME));
		options.addOption(new Option(Thunder.OPT_CHOICE_DB_FORCE_REFRESH, "reads and adds the annotation to the database, regardless of whether the data already exists"));
		return options;
	}


	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {

		//		ReadVCF.ReadVCFToDB("/Users/robk/Desktop/TEST.VCF", "/Users/robk/Desktop/TEST_SAMPLE_DB", "SAMPLE");
		//		ReadVCF.ReadVCFToDB("/Users/robk/Desktop/HSB_123_SNPchip_ALLSNPS.vcf", "/Users/robk/Desktop/TEST_SAMPLE_DB", "SAMPLE", true);

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		boolean forceNewDB = false;
		if(cmdArgs.hasOption(Thunder.OPT_CHOICE_DB_FORCE_REFRESH)){
			forceNewDB = true;
		}

		if(cmdArgs.hasOption(Thunder.OPT_PATH_INPUT) && cmdArgs.hasOption(Thunder.OPT_PATH_DB_SAMPLE)){
			ReadVCF.ReadVCFToDB(cmdArgs.getOptionValue(Thunder.OPT_PATH_INPUT), cmdArgs.getOptionValue(Thunder.OPT_PATH_DB_SAMPLE), "SAMPLE", forceNewDB); 
		}else{
			HelpFormatter formatter = new HelpFormatter();
			//formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ReadVCF", getCmdLineOptions());
			formatter.printHelp("ReadVCF", getCmdLineOptions());
			System.out.println();
		}


	}

}
