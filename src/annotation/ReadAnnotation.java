package annotation;

import java.io.File;
import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.ArrayList;
import java.util.Iterator;

import main.Thunder;
import objects.GTF;
import objects.GenomicCoordinate;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import database.DBConnect;
import database.DBConnect_SQLite;

public class ReadAnnotation {


	public static void addToDB(GTF gtfContents, String dbPath, String tableName, boolean forceNewDB) throws Exception{
		DBConnect db = new DBConnect_SQLite(dbPath);
		//		DBConnect db = new DBConnect_H2(dbPath);

		//		System.out.println("db.containsTable("+tableName+") = "+db.containsTable(tableName));
		//		System.out.println("gtfContents.getType() = "+gtfContents.getType());


		// If the annotation does not exist in this DB, read the GTF and create it
		if(db.containsTable(tableName)==false || forceNewDB){

			// Create the table to hold all of the GTF info:
			Statement st = db.createStatement();
			st.execute("DROP TABLE IF EXISTS "+tableName);
			st.execute("DROP TABLE IF EXISTS "+tableName+"_geneLevel");
			st.execute("DROP TABLE IF EXISTS "+tableName+"_rawGTF");
			st.close();

			System.out.print("Adding all GTF entries to the database...");
			long time = System.currentTimeMillis();
			addAllGTFData(db, tableName, gtfContents);
			System.out.println("Done ("+((System.currentTimeMillis()-time)/1000.0)+"s)");

			gtfContents = GTF.collapseEntriesByTranscriptID(gtfContents);

			if(gtfContents.getType() == GTF.GTF_TYPE_GENCODE){
				System.out.print("Populating the database...");
				time = System.currentTimeMillis();
				System.out.print("Sorting and collapsing GTF entries by geneID...");
				GTF gtfGenes = GTF.collapseEntriesByGeneID(gtfContents);
				populateDB_GENCODE(db, tableName, gtfContents, gtfGenes);
				System.out.println("Done ("+((System.currentTimeMillis()-time)/1000.0)+"s)");
			}else if(gtfContents.getType() == GTF.GTF_TYPE_CUFFMERGE){
				System.out.print("Sorting and collapsing GTF entries by geneID...");
				GTF gtfGenes = GTF.collapseEntriesByGeneID(gtfContents);
				System.out.println("Done.");

				System.out.print("Populating the database...");
				time = System.currentTimeMillis();
				populateDB_CUFFMERGE(db, tableName, gtfContents, gtfGenes);
				System.out.println("Done ("+((System.currentTimeMillis()-time)/1000.0)+"s)");
			}

		}
		db.closeConnection();
	}



	/**
	 * Add GTF entries from a CUFFLINKS annotation file to the DB
	 * 
	 * @param ps
	 * @param it_2
	 */
	private static void addGTFEntries_Cufflinks(PreparedStatement ps, Iterator<GenomicCoordinate> it_2){
		GenomicCoordinate tmpGC;

		while(it_2.hasNext()){
			try{
				tmpGC = it_2.next();

				ps.setString(1, tmpGC.getChrom());
				ps.setString(2, tmpGC.getSource());
				ps.setString(3, tmpGC.getFeatureType());
				ps.setInt(4, tmpGC.getStart());
				ps.setInt(5, tmpGC.getStop());
				//			ps.setDouble(6, null);
				ps.setString(7, tmpGC.getStrand());
				//			ps.setInt(8, null);
				ps.setString(9, tmpGC.getAttribute("gene_id"));
				ps.setString(10, tmpGC.getAttribute("transcript_id"));				
				ps.setString(11, tmpGC.getAttribute("gene_name"));
				ps.setString(12, tmpGC.getAttribute("oId"));

				ps.addBatch();
				ps.executeBatch();
			}catch(SQLException e){
				e.printStackTrace();
			}
		}
	}

	/**
	 * Add GTF entries from a GENCODE annotation file to the DB
	 * 
	 * @param ps
	 * @param it_2
	 */
	private static void addGTFEntries_Gencode(PreparedStatement ps, Iterator<GenomicCoordinate> it_2){
		GenomicCoordinate tmpGC;

		while(it_2.hasNext()){
			try{
				tmpGC = it_2.next();

				ps.setString(1, tmpGC.getChrom());
				ps.setString(2, tmpGC.getSource());
				ps.setString(3, tmpGC.getFeatureType());
				ps.setInt(4, tmpGC.getStart());
				ps.setInt(5, tmpGC.getStop());
				//			ps.setDouble(6, null);
				ps.setString(7, tmpGC.getStrand());
				//			ps.setInt(8, null);
				ps.setString(9, tmpGC.getAttribute("gene_id"));
				ps.setString(10, tmpGC.getAttribute("transcript_id"));				
				ps.setString(11, tmpGC.getAttribute("gene_name"));
				//				ps.setString(12, tmpGC.getAttribute("oId"));

				ps.addBatch();
				ps.executeBatch();
			}catch(SQLException e){
				e.printStackTrace();
			}
		}
	}

	/**
	 * Create a DB table to hold the 'raw' GTF data (individual exons etc.) and populate using the GTF object
	 * @param db
	 * @param tableName
	 * @param gtfData
	 * @throws SQLException
	 */
	private static void addAllGTFData(DBConnect db, String tableName, GTF gtfData) throws SQLException{
		Statement st = db.createStatement();

		Iterator<String> it;
		String tmpChromosome;

		// Create DB table
		st.execute("create table "+tableName+"_rawGTF(" +
				"chromosome varchar(30) not null, " +
				"source varchar(30) not null, " +
				"feature varchar(30) not null, " +
				"start int not null, " +
				"stop int not null, " +
				"score double null, " +
				"strand char(1) not null, " +
				"frame int null, " +
				"geneID varchar(30) not null, " +
				"transcriptID varchar(30) not null, " +
				"geneName varchar(30), " +
				"transcriptID_ensembl varchar(30) null)");//, " +
		//				"primary key (chromosome,feature,start,stop))");

		PreparedStatement ps = db.getPreparedStatement("INSERT INTO "+tableName+"_rawGTF VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
		db.setAutoCommit(false);

		
		it = gtfData.getCoordinates().keySet().iterator();

		// For each chromosome:
		while(it.hasNext()){
			tmpChromosome = it.next();
			
			// Saves time to only check the GTF source once per chomosome:
			if(gtfData.getCoordinates().get(tmpChromosome).get(0).getSource().equalsIgnoreCase("cufflinks")){
				addGTFEntries_Cufflinks(ps, gtfData.getCoordinates().get(tmpChromosome).iterator());
			}else{
				addGTFEntries_Gencode(ps, gtfData.getCoordinates().get(tmpChromosome).iterator());
			}
		}
		db.commit();
		db.setAutoCommit(true);

		st.execute("DROP INDEX IF EXISTS index_"+tableName+"_rawGTF_chrom");
		st.execute("DROP INDEX IF EXISTS index_"+tableName+"_rawGTF_feature");
		st.execute("DROP INDEX IF EXISTS index_"+tableName+"_rawGTF_start");
		st.execute("DROP INDEX IF EXISTS index_"+tableName+"_rawGTF_stop");
		st.execute("CREATE INDEX index_"+tableName+"_rawGTF_chrom ON "+tableName+"_rawGTF(chromosome)");
		st.execute("CREATE INDEX index_"+tableName+"_rawGTF_feature ON "+tableName+"_rawGTF(feature)");
		st.execute("CREATE INDEX index_"+tableName+"_rawGTF_start ON "+tableName+"_rawGTF(start)");
		st.execute("CREATE INDEX index_"+tableName+"_rawGTF_stop ON "+tableName+"_rawGTF(stop)");

		st.close();
		ps.close();
	}





	/**
	 * 
	 * @param db
	 * @param tableName
	 * @param gtfTranscripts
	 * @throws SQLException
	 */
	private static void populateDB_GENCODE(DBConnect db, String tableName, GTF gtfTranscripts, GTF gtfGenes) throws SQLException{

		Statement st = db.createStatement();
		PreparedStatement ps;


		//
		// Add gene-level data
		//
		addGeneLevelData(db, tableName, gtfGenes, "transcript_id");


		//
		// Add transcript-level data
		//
		// chr8	HAVANA	transcript	101960842	101965111	.	-	.	gene_id "ENSG00000164924.12"; transcript_id "ENST00000437293.1"; gene_type "protein_coding"; gene_status "KNOWN"; gene_name "YWHAZ"; transcript_type "protein_coding"; transcript_status "KNOWN"; transcript_name "YWHAZ-006"; level 2; tag "alternative_5_UTR"; tag "mRNA_end_NF"; tag "cds_end_NF"; havana_gene "OTTHUMG00000134291.4"; havana_transcript "OTTHUMT00000259022.1";
		st.execute("create table "+tableName+"(" +
				"transcriptID varchar(30) not null, " +
				"chromosome varchar(30) not null, " +
				"start int not null, " +
				"stop int not null, " +
				"strand char(1) not null, " +
				"geneID varchar(30) not null, " +
				"geneName varchar(30), " +
				"geneType varchar(30), " +
				"geneStatus varchar(50), " +
				"transcriptName varchar(30), " +
				"transcriptType varchar(50), " +
				"transcriptStatus varchar(50), " +
				"transcriptID_HAVANA varchar(30), " +
				"geneID_HAVANA varchar(30), " +
				"primary key (transcriptID))");


		db.setAutoCommit(false);
		ps = db.getPreparedStatement("INSERT INTO "+tableName+" VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");

		Iterator<String> it = gtfTranscripts.getCoordinates().keySet().iterator();
		String tmpChromosome;
		GenomicCoordinate tmpGC;
		while(it.hasNext()){
			tmpChromosome = it.next();
			Iterator<GenomicCoordinate> it_2 = gtfTranscripts.getCoordinates().get(tmpChromosome).iterator();
			while(it_2.hasNext()){
				tmpGC = it_2.next();

				ps.setString(1, tmpGC.getCoordinateID());
				ps.setString(2, tmpGC.getChrom());
				ps.setInt(3, tmpGC.getStart());
				ps.setInt(4, tmpGC.getStop());
				ps.setString(5, tmpGC.getStrand());
				ps.setString(6, tmpGC.getAttribute("gene_id"));
				ps.setString(7, tmpGC.getAttribute("gene_name"));
				ps.setString(8, tmpGC.getAttribute("gene_type"));
				ps.setString(9, tmpGC.getAttribute("gene_status"));
				ps.setString(10, tmpGC.getAttribute("transcript_name"));
				ps.setString(11, tmpGC.getAttribute("transcript_type"));
				ps.setString(12, tmpGC.getAttribute("transcript_status"));
				ps.setString(13, tmpGC.getAttribute("havana_transcript"));
				ps.setString(14, tmpGC.getAttribute("havana_gene"));
				ps.addBatch();
			}
		}
		ps.executeBatch();
		db.commit();
		db.setAutoCommit(true);

		st.execute("DROP INDEX IF EXISTS index_"+tableName+"_geneID");
		st.execute("CREATE INDEX index_"+tableName+"_geneID ON "+tableName+"(geneID)");

		st.close();
		ps.close();
	}






	/**
	 * Method to create the gene-level summary of each annotation
	 * 
	 * @param db
	 * @param tableName
	 * @param gtfGenes
	 * @param variableName_transcriptID
	 * @throws SQLException
	 */
	private static void addGeneLevelData(DBConnect db, String tableName, GTF gtfGenes, String variableName_transcriptID) throws SQLException{
		Statement st = db.createStatement();
		PreparedStatement ps;

		Iterator<String> it;
		String tmpChromosome;
		GenomicCoordinate tmpGC;

		//
		// Add gene-level data
		//
		st.execute("create table "+tableName+"_geneLevel(" +
				"geneID varchar(30) not null, " +
				"chromosome varchar(30) not null, " +
				"start int not null, " +
				"stop int not null, " +
				"strand char(1) not null, " +
				"geneName varchar(30), " +
				"transcriptID_ensembl varchar(30), " +
				"primary key (geneID))");

		db.setAutoCommit(false);
		ps = db.getPreparedStatement("INSERT INTO "+tableName+"_geneLevel VALUES (?, ?, ?, ?, ?, ?, ?)");
		it = gtfGenes.getCoordinates().keySet().iterator();
		int errCount = 0;
		while(it.hasNext()){
			tmpChromosome = it.next();
			Iterator<GenomicCoordinate> it_2 = gtfGenes.getCoordinates().get(tmpChromosome).iterator();
			while(it_2.hasNext()){
				tmpGC = it_2.next();

				ps.setString(1, tmpGC.getAttribute("gene_id"));
				ps.setString(2, tmpGC.getChrom());
				ps.setInt(3, tmpGC.getStart());
				ps.setInt(4, tmpGC.getStop());
				ps.setString(5, tmpGC.getStrand());
				ps.setString(6, tmpGC.getAttribute("gene_name"));
				ps.setString(7, tmpGC.getAttribute(variableName_transcriptID));				
				ps.addBatch();
				try{
					ps.executeBatch();
				}catch(SQLException e){
					errCount++;
					System.err.println("Warning: Duplicate GeneID ("+errCount+"): "+tmpGC.getAttribute("gene_id"));
					//e.printStackTrace();
				}
			}
		}
		db.commit();
		db.setAutoCommit(true);

		st.execute("DROP INDEX IF EXISTS index_"+tableName+"_geneLevel_transcripID_ensembl");
		st.execute("CREATE INDEX index_"+tableName+"_geneLevel_transcripID_ensembl ON "+tableName+"_geneLevel(transcriptID_ensembl)");

		st.close();
		ps.close();

	}


	/**
	 * 
	 * @param db
	 * @param tableName
	 * @param gtfTranscripts
	 * @param gtfGenes
	 * @throws SQLException
	 */
	private static void populateDB_CUFFMERGE(DBConnect db, String tableName, GTF gtfTranscripts, GTF gtfGenes) throws SQLException{
		Statement st = db.createStatement();
		PreparedStatement ps;

		Iterator<String> it;
		String tmpChromosome;
		GenomicCoordinate tmpGC;

		//
		// Add gene-level data
		//
		addGeneLevelData(db, tableName, gtfGenes, "oId");


		//
		// Add transcript-level data
		//
		// chr8	Cufflinks	exon	101928753	101932980	.	-	.	gene_id "XLOC_044997"; transcript_id "TCONS_00286660"; exon_number "1"; gene_name "YWHAZ"; oId "ENST00000395957.2"; nearest_ref "ENST00000395957.2"; class_code "="; tss_id "TSS117693"; p_id "P74517";
		st.execute("create table "+tableName+"(" +
				"transcriptID varchar(30) not null, " +
				"chromosome varchar(30) not null, " +
				"start int not null, " +
				"stop int not null, " +
				"strand char(1) not null, " +
				"geneID varchar(30) not null, " +
				"geneName varchar(30), " +
				"transcriptID_ensembl varchar(30), " +
				"class_code varchar(5), " +
				"tss_id varchar(30), " +
				"p_id varchar(30), " +
				"primary key (transcriptID))");

		db.setAutoCommit(false);
		// chr8	Cufflinks	exon	101928753	101932980	.	-	.	gene_id "XLOC_044997"; transcript_id "TCONS_00286660"; exon_number "1"; gene_name "YWHAZ"; oId "ENST00000395957.2"; nearest_ref "ENST00000395957.2"; class_code "="; tss_id "TSS117693"; p_id "P74517";
		ps = db.getPreparedStatement("INSERT INTO "+tableName+" VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
		it = gtfTranscripts.getCoordinates().keySet().iterator();
		while(it.hasNext()){
			tmpChromosome = it.next();
			Iterator<GenomicCoordinate> it_2 = gtfTranscripts.getCoordinates().get(tmpChromosome).iterator();
			while(it_2.hasNext()){
				tmpGC = it_2.next();

				ps.setString(1, tmpGC.getAttribute("transcript_id"));
				ps.setString(2, tmpGC.getChrom());
				ps.setInt(3, tmpGC.getStart());
				ps.setInt(4, tmpGC.getStop());
				ps.setString(5, tmpGC.getStrand());
				ps.setString(6, tmpGC.getAttribute("gene_id"));
				ps.setString(7, tmpGC.getAttribute("gene_name"));
				ps.setString(8, tmpGC.getAttribute("oId"));
				ps.setString(9, tmpGC.getAttribute("class_code"));
				ps.setString(10, tmpGC.getAttribute("tss_id"));
				ps.setString(11, tmpGC.getAttribute("p_id"));
				ps.addBatch();
				ps.executeBatch();
			}
		}
		db.commit();
		db.setAutoCommit(true);

		st.execute("DROP INDEX IF EXISTS index_"+tableName+"_geneID");
		st.execute("DROP INDEX IF EXISTS index_"+tableName+"_transcripID_ensembl");
		st.execute("CREATE INDEX index_"+tableName+"_geneID ON "+tableName+"(geneID)");
		st.execute("CREATE INDEX index_"+tableName+"_transcripID_ensembl ON "+tableName+"(transcriptID_ensembl)");

		st.close();
		ps.close();
	}





	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("annotationPath").hasArg().withDescription("path to the file containing the GENCODE annotation data in GTF format").create(Thunder.OPT_PATH_ANNOTATION));
		options.addOption(OptionBuilder.withArgName("databasePath").hasArg().withDescription("create or add to the database at this path\n(e.g. /path/to/db)").create(Thunder.OPT_PATH_DB_ANNOTATION));
		options.addOption(OptionBuilder.withArgName("tableName").hasArg().withDescription("name of the table to which to add the annotation data").create(Thunder.OPT_DB_TABLE_NAME));
		options.addOption(new Option(Thunder.OPT_CHOICE_DB_FORCE_REFRESH, "reads and adds the annotation to the database, regardless of whether the data already exists"));
		return options;
	}


	public static void main(String[] args) throws Exception {


		//		args = new String[]{"ReadAnnotation","-F","-A","/Users/robk/Work/YALE_OFFLINE/ANNOTATIONS/gencode14.sqlite","-a","/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v14.annotation_noSelenocysteine.gtf"};
		//args = new String[]{"ReadAnnotation","-F","-A","/Users/robk/Work/YALE_OFFLINE/ANNOTATIONS/gencode18.sqlite","-a","/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v18.annotation.gtf"};
//		args = new String[]{"ReadAnnotation","-F","-A","/Users/robk/Work/YALE_OFFLINE/ANNOTATIONS/gencode14.sqlite","-a","/Users/robk/Dropbox/Work/YALE/My_Proteomics/HEKcells/mRNA/RNA-seq/CuffdiffResults/merged.gtf"};
		

		
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		boolean forceNewDB = false;
		if(cmdArgs.hasOption(Thunder.OPT_CHOICE_DB_FORCE_REFRESH)){
			forceNewDB = true;
		}

		if(cmdArgs.hasOption(Thunder.OPT_PATH_ANNOTATION) && cmdArgs.hasOption(Thunder.OPT_PATH_DB_ANNOTATION)){
			ArrayList<String> features = new ArrayList<String>();
//			features.add("exon");

			GTF gtfContents = ReadGTF.readGTF(cmdArgs.getOptionValue(Thunder.OPT_PATH_ANNOTATION), true, features, false, false, null);

			File tmp = new File(cmdArgs.getOptionValue(Thunder.OPT_PATH_ANNOTATION));
			String tableName = tmp.getName().replaceAll("\\.", "_");

			//			System.out.println("tableName = "+tableName);

			if(cmdArgs.hasOption(Thunder.OPT_DB_TABLE_NAME))
				tableName = cmdArgs.getOptionValue(Thunder.OPT_DB_TABLE_NAME);

			if(gtfContents.getType() == 0){
				System.out.println("Error: Unable to determine the source/type of this annotation file.  Be sure the header lines are intact"); 
			}else{
				addToDB(gtfContents, cmdArgs.getOptionValue(Thunder.OPT_PATH_DB_ANNOTATION), tableName, forceNewDB);
			}

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ReadAnnotation", getCmdLineOptions());
			System.out.println();
		}
	}



}
