package main;

import org.apache.commons.cli.BasicParser;
import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import proteome.ProcessPepDigest;
import proteome.ReadXTandem;
import proteome.ReadmzXML;
import transcriptome.CIGAR_2_PWM;
import transcriptome.ExtractSequencesFromGTF;
import transcriptome.ModifyHeadersInTophatTranscriptome;
import annotation.ReadAnnotation;
import database.ListDBTables;
import fastaTools.Fasta2Fastq;
import fastaTools.FastaHeaderGrep;
import fastaTools.FastaRemoveDuplicates;
import fastqTools.Fastq2Fasta;
import fastqTools.FilterBySequenceLength;
import fastqTools.FilterFastxByHeaderList;
import fastqTools.FindAdapter;
import fastqTools.GetSequenceLengths;
import fastqTools.MatchPairedEndSequences;
import fastqTools.ProcessFastqWithRandomBarcode;
import fastqTools.RemoveHomopolymers;
import footprintAlignments.FootprintFrameCalculator;
import footprintAlignments.ReadCoverage;
import genome.CreatePersonalGenome;
import genome.ReadVCF;

public class Thunder {

	public static final String VERSION = "0.6.1";
	
	public static final String OPT_PATH_DB_ANNOTATION = "A";
	public static final String OPT_PATH_DB_SAMPLE = "S";
	public static final String OPT_DB_TABLE_NAME = "t";
	public static final String OPT_PATH_INPUT = "i";
	public static final String OPT_PATH_OUTPUT = "o";
	public static final String OPT_PATH_SPECTRA_TANDEM = "i";
	public static final String OPT_PATH_SPECTRA_MZXML = "i";
	public static final String OPT_PATH_ANNOTATION = "a";
	public static final String OPT_FORMAT_GENCODE = "g";
	public static final String OPT_FORMAT_CUFFMERGE = "c";
	public static final String OPT_FORMAT_SWISSPROT = "s";
	public static final String OPT_FORMAT_TRINITY = "n";
	public static final String OPT_CHOICE_DB_FORCE_REFRESH = "F";

	/**
	 * Parse the command line arguments
	 * @param args
	 * @return
	 * @throws ParseException
	 */
	public static CommandLine parseArgs(String[] args, Options options) throws ParseException{
		CommandLineParser parser = new BasicParser();
		return parser.parse(options, args);	
	}



	public static void main(String[] args) throws Exception {
		String main = "NULL";
		if(args.length > 0){
			main = args[0].trim().toLowerCase(); 
		}

		//		main = "getsequencelengths";
		//		main = "filtersequencesbylength";
		//		main = "modifytophatheaders";
		//		main = "listdbtables";
		//		main = "extractsequencesfromgtf";
		//		main = "readannotation";
		//		main = "readvcf";
		//		main = "createpersonalgenome";
		//		main = "readmzxml";
		//		main = "readtandem"; 

		if(main.equals("getsequencelengths")){
			GetSequenceLengths.main(args);
		}else if(main.equals("filtersequencesbylength")){
			FilterBySequenceLength.main(args);
		}else if(main.equals("removehomopolymerrepeats")){
			RemoveHomopolymers.main(args);
		}else if(main.equals("modifytophatheaders")){
			ModifyHeadersInTophatTranscriptome.main(args);
		}else if(main.equals("listdbtables")){
			ListDBTables.main(args);
		}else if(main.equals("gtf2fasta")){
			ExtractSequencesFromGTF.main(args);
		}else if(main.equals("readannotation")){
			ReadAnnotation.main(args);
		}else if(main.equals("readvcf")){
			ReadVCF.main(args);
		}else if(main.equals("createpersonalgenome")){
			CreatePersonalGenome.main(args);
		}else if(main.equals("readmzxml")){
			ReadmzXML.main(args);
		}else if(main.equals("parsetandemoutput")){
			ReadXTandem.main(args);
		/*}else if(main.equals("isoformem_footprints")){
			//RibosomeFootprintEM.main(args);
			FootprintEM.main(args);
		}else if(main.equals("isoformem_proteomics")){
			//RibosomeFootprintEM.main(args);
			ProteomicEM.main(args);
		*/}else if(main.equals("cigar_2_pwm")){
			CIGAR_2_PWM.main(args);
		}else if(main.equals("matchpairedendsequences")){
			MatchPairedEndSequences.main(args);
		}else if(main.equals("fastaheadergrep")){
			FastaHeaderGrep.main(args);
		//}else if(main.equals("processendogenousalignments")){
		//	ProcessEndogenousAlignments.main(args);
		}else if(main.equals("filterfastxbyidlist")){
			FilterFastxByHeaderList.main(args);
		}else if(main.equals("processfastqwithrandombarcode")){
			ProcessFastqWithRandomBarcode.main(args);
		}else if(main.equals("findadapter")){
			FindAdapter.main(args);
		}else if(main.equals("footprintframeanaysis")){
			FootprintFrameCalculator.main(args);
		}else if(main.equals("peptidedigest")){
			ProcessPepDigest.main(args);
		}else if(main.equals("readcoverage")){
			ReadCoverage.main(args);
		}else if(main.equals("fastq2fasta")){
			Fastq2Fasta.main(args);
		}else if(main.equals("fasta2fastq")){
			Fasta2Fastq.main(args);
		}else if(main.equals("removefastaduplicates")){
			FastaRemoveDuplicates.main(args);
		}
		
		
		
		
		else{
			System.out.println("Thunder version "+VERSION);
			System.out.println("");
			
			System.out.println("Usage:\t "+Thunder.THUNDER_EXE_COMMAND+" <Command>");
			//System.out.println("Usage:\t java -Xmx2G -jar Thunder.jar <Command>");
			System.out.println("");

			System.out.println("Command: GetSequenceLengths            | Get the distribution of sequence lengths in a FASTA/Q file");
			System.out.println("         FastaHeaderGrep               | Filter fasta sequences based on the sequence ID");
			System.out.println("         RemoveFastaDuplicates         | Filter fasta sequences to remove duplicate sequence headers");
			System.out.println("         FilterFastxByIDList           | Filter fasta/q sequences based on a list of sequence IDs");
			System.out.println("         FilterSequencesByLength       | Filter fasta or fastq sequences based on some maximum sequence length");
			System.out.println("         ProcessFastqWithRandomBarcode | Filter fasta or fastq sequences based on some maximum sequence length");
			System.out.println("         FindAdapter                   | Determine most likely 3' adapter sequence from fastq reads");
			System.out.println("         RemoveHomopolymerRepeats      | Filter fasta or fastq sequences based on sequence composition");
			System.out.println("         MatchPairedEndSequences       | Match paired-end fastq sequences based on readID");
			//System.out.println("         ModifyTophatHeaders        | re-formats fasta headers in the tophat-generated transcriptome for downstream compatibility with Thunder");
			//System.out.println("         ListDBTables               | list the contents of a database");
			System.out.println("         GTF2Fasta                     | Extract GTF coordinates from FASTA sequence(s)");
			System.out.println("         Fasta2Fastq                   | Convert FASTA sequence(s) to FASTQ sequence(s)");
			System.out.println("         Fastq2Fasta                   | Convert FASTQ sequence(s) to FASTA sequence(s)");
			//System.out.println("         ReadAnnotation             | import annotation information in GTF to the database");
			//System.out.println("         ReadVCF                    | import genotype data to the database");
			//System.out.println("         CreatePersonalGenome       | create a personalised reference genome for a sample in the database");
			//System.out.println("         ReadMzXML                  | import MS- and MS/MS-spectra to the database");
			//System.out.println("         ReadTandem_toDB            | import X!Tandem MS/MS-spectra IDs to the database");
			System.out.println("         ParseTandemOutput             | Process X!Tandem MS/MS-spectra alignments and output summary table");
			//System.out.println("         IsoformEM_Footprints          | Infer most likely transcripts from ribosome footprint alignments");
			//System.out.println("         IsoformEM_Proteomics          | Infer most likely isoforms from MS/MS spectra mapping");
			System.out.println("         CIGAR_2_PWM                   | Reads SAM alignments and converts the CIGAR strings to a position-weight matrix");
			System.out.println("         ReadCoverage                  | Reads SAM/BAM alignments to the TRANSCRIPTOME and calculates read coverage consistency");
			//System.out.println("         ProcessEndogenousAlignments   | Process endogenous smallRNA alignments for the exceRpt pipeline");
			System.out.println("         FootprintFrameAnaysis         | Analyse ribosome footprint alignments in terms of fidelity to annotated coding frames");
			System.out.println("         PeptideDigest                 | Enzymatically digest a protein reference");
			System.out.println();
		}


	}

	public static final String THUNDER_EXE_COMMAND = "java -Xmx2G -jar Thunder.jar";

	
	
}
