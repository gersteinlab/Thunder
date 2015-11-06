package genome;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import database.DBConnect;
import database.DBConnect_SQLite;

public class CreatePersonalGenome {

	private BufferedReader fastaReader;
	private BufferedWriter fastaWriter;
	private String nextChromosomeID = null;
	private Statement st;

	
	/**
	 * 
	 * @return
	 * @throws IOException
	 */
	private StringBuilder readChromosome() throws IOException{
		String line;
		StringBuilder sb = new StringBuilder();
		while((line=fastaReader.readLine()) != null){
			if(line.startsWith(">")){
				break;
			}else{
				sb.append(line.trim());
			}

		}

		if(line != null)
			nextChromosomeID = line.substring(1);
		else
			nextChromosomeID = null;

		return sb;
	}


	/**
	 * 
	 * @param nt1
	 * @param nt2
	 * @return
	 */
	public static String getAmbiguityCode(String nt1, String nt2){
		String code = "";
		String n1, n2;

		if(nt1.length() == 1 && nt2.length() == 1){
			int comp = nt1.compareToIgnoreCase(nt2);
			if(comp < 0){
				n1 = nt1;
				n2 = nt2;
			}else if(comp > 0){
				n1 = nt2;
				n2 = nt1;
			}else{
				// same base input twice, do nothing
				return null;
			}

			if(n1.equalsIgnoreCase("A")){
				if(n2.equalsIgnoreCase("C"))
					code = "M";
				else if(n2.equalsIgnoreCase("G"))
					code = "R";
				else if(n2.equalsIgnoreCase("T"))
					code = "W";
			}else if(n1.equalsIgnoreCase("C")){
				if(n2.equalsIgnoreCase("G"))
					code = "S";
				else if(n2.equalsIgnoreCase("T"))
					code = "Y";
			}else if(n1.equalsIgnoreCase("G")){
				if(n2.equalsIgnoreCase("T"))
					code = "K";
			}
			return code;

		}else{
			return null;
		}
	}

	private int countHomSNPs = 0;
	private int countHetSNPs = 0;
	
	/**
	 * 
	 * @param refChromosome
	 * @param chromosomeID
	 * @return
	 * @throws SQLException
	 */
	private StringBuilder modifyChromosome(StringBuilder refChromosome, String chromosomeID) throws SQLException{

		ResultSet rs = st.executeQuery("SELECT * FROM genomicVariants WHERE chromosome = '"+chromosomeID+"' ORDER BY position DESC");
		System.out.print("...");

		String genotype, refBase, altBase;
		int pos;
		
		while(rs.next()){
			pos = rs.getInt(2);
			refBase = refChromosome.charAt(pos-1)+"";
			altBase = rs.getString(5).toUpperCase();
			genotype = rs.getString(7);

//			System.out.print(rs.getString(1)+"\t"+pos+"\t"+rs.getString(3)+"\t"+rs.getString(4)+"\t"+rs.getString(5)+"\t"+genotype+"\t"+refChromosome.substring(pos-3, pos+2));
			
			boolean refIsLowerCase = false;
			if(refBase.equals("a") || refBase.equals("c") || refBase.equals("g") || refBase.equals("t") || refBase.equals("n"))
				refIsLowerCase = true;

			if(genotype.equals("1/1")){ // homo SNP
				if(refIsLowerCase)
					altBase = altBase.toLowerCase();
//				refChromosome.replace(pos-1, pos, altBase); // SLOW!!!
				refChromosome.setCharAt(pos-1, altBase.toCharArray()[0]);
				countHomSNPs++;
			}else if(genotype.equals("1/0") || genotype.equals("0/1")){ // het SNP
				/*altBase = getAmbiguityCode(refBase.toUpperCase(), altBase.toUpperCase());
				if(refIsLowerCase)
					altBase = altBase.toLowerCase();
//				refChromosome.replace(pos-1, pos, altBase); // SLOW!!!
				refChromosome.setCharAt(pos-1, altBase.toCharArray()[0]);*/
				countHetSNPs++;
			}
			
//			System.out.println("\t"+refBase+"\t"+altBase+"\t"+refChromosome.substring(pos-3, pos+2));
			
//			count ++;
//			if(count % 1000 == 0)
//				System.out.print(".");
		}

		rs.close();

		return refChromosome;
	}
	
	
	
	
	/**
	 * 
	 * @param refGenome
	 * @param refAnnotation
	 * @param dbPath
	 * @param outputGenome
	 * @throws Exception
	 */
	public CreatePersonalGenome(String refGenome, String refAnnotation, String dbPath, String outputGenome) throws Exception{
		DBConnect db = new DBConnect_SQLite(dbPath);
		st = db.createStatement();
		db.setAutoCommit(false);

		fastaReader = new BufferedReader(new FileReader(refGenome));	
		fastaWriter = new BufferedWriter(new FileWriter(outputGenome));

		// Read first line (should be ID for the 1st fasta entry)
		nextChromosomeID = fastaReader.readLine().substring(1);
		String thisChromosomeID;
		StringBuilder thisChromosome;

		while(nextChromosomeID != null){
			thisChromosomeID = nextChromosomeID; 
			System.out.print("Chromosome '"+thisChromosomeID+"': Reading...");
			thisChromosome = readChromosome();
			System.out.print("modifying");
			thisChromosome = modifyChromosome(thisChromosome, thisChromosomeID);
			System.out.print("writing...");
			writeChromosome(thisChromosome, thisChromosomeID);
			System.out.println("Done.");
		}

		System.out.println();
		System.out.println("Total homozygous SNPs:   "+countHomSNPs);
		System.out.println("Total heterozygous SNPs: "+countHetSNPs);
		System.out.println();
		
		
		// tidy up
		db.setAutoCommit(true);
		st.close();
		db.closeConnection();
		fastaReader.close();
		fastaWriter.flush();
		fastaWriter.close();
	}

	
	
	/**
	 * 
	 * @param thisChromosome
	 * @param chromosomeID
	 * @throws IOException
	 */
	public void writeChromosome(StringBuilder thisChromosome, String chromosomeID) throws IOException{
		fastaWriter.write(">"+chromosomeID+"\n");

		int linelength = 50;
		int nLines = (int) Math.round(Math.floor(thisChromosome.length() / linelength));
		for(int i=0;i<nLines;i++){
			fastaWriter.write(thisChromosome.substring(i*linelength, (i+1)*(linelength))+"\n");	
		}
		fastaWriter.write(thisChromosome.substring(nLines*linelength, thisChromosome.length())+"\n");
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("genome").hasArg().withDescription("path to the genome reference sequence in fasta format").create(Thunder.OPT_PATH_INPUT));
		options.addOption(OptionBuilder.withArgName("annotation").hasArg().withDescription("path to the file containing the GENCODE annotation data in GTF format").create(Thunder.OPT_PATH_ANNOTATION));
		options.addOption(OptionBuilder.withArgName("databasePath").hasArg().withDescription("obtain variants for the sample using the database at this path (e.g. /path/to/db)").create(Thunder.OPT_PATH_DB_SAMPLE));
		options.addOption(OptionBuilder.withArgName("output").hasArg().withDescription("output the personal genome at this path").create(Thunder.OPT_PATH_OUTPUT));
		return options;
	}



	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{

//		args = new String[]{"-i","/Users/robk/WORK/YALE_offline/ANNOTATIONS/hg19.fa","-a","/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v14.annotation_noSelenocysteine.gtf",
//				"-S","/Users/robk/Desktop/TEST_SAMPLE_DB","-o","/Users/robk/Desktop/personalGenome.fa"};


		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption(Thunder.OPT_PATH_INPUT) && cmdArgs.hasOption(Thunder.OPT_PATH_ANNOTATION) && 
				cmdArgs.hasOption(Thunder.OPT_PATH_DB_SAMPLE) && cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT)){
			new CreatePersonalGenome(cmdArgs.getOptionValue(Thunder.OPT_PATH_INPUT), cmdArgs.getOptionValue(Thunder.OPT_PATH_ANNOTATION), cmdArgs.getOptionValue(Thunder.OPT_PATH_DB_SAMPLE), cmdArgs.getOptionValue(Thunder.OPT_PATH_OUTPUT)); 
		}else{
			HelpFormatter formatter = new HelpFormatter();
			//formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" CreatePersonalGenome", getCmdLineOptions());
			formatter.printHelp("CreatePersonalGenome", getCmdLineOptions());
			System.out.println();
		}


	}

}
