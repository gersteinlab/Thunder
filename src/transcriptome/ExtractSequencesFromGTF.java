package transcriptome;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;

import main.Thunder;
import objects.GTF;
import objects.GenomicCoordinate;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import annotation.ReadGTF;

public class ExtractSequencesFromGTF {


	private Sextractor sextractor;
	//	private ArrayList<String> transcriptIDsToInclude = new ArrayList<String>();
	private HashSet<String> transcriptIDsToInclude = new HashSet<String>();
	//	private HashMap<Integer, String> gtfDump = new HashMap<Integer, String>();
	private ArrayList<String> featuresToKeep = new ArrayList<String>();
	private ArrayList<String> additionalAttributes = new ArrayList<String>();

	public ExtractSequencesFromGTF(){
		this.sextractor = new Sextractor();
	}

	public void defineFeatureList(String features){
		String[] tmp = features.trim().split(",");
		for(int i=0;i<tmp.length;i++)
			featuresToKeep.add(tmp[i].trim());
	}
	
	public void defineAdditionalAttributes(String attributes){
		String[] tmp = attributes.trim().split(",");
		for(int i=0;i<tmp.length;i++)
			additionalAttributes.add(tmp[i].trim());
	}

	public String trimTranscriptID(String in){
		return in.substring(1, in.length()-2);
	}

	/**
	 * Allows reading of a file containing transcript IDs to use for sequence extraction
	 *  
	 * @param transcriptIDs_toInclude Path of the include list
	 * @throws IOException
	 */
	public void readIncludeList(String transcriptIDs_toInclude) throws IOException{

		System.out.print("Reading transcriptIDs to include...");

		BufferedReader in = new BufferedReader(new FileReader(transcriptIDs_toInclude));
		String line = "";

		while((line=in.readLine()) != null){
			this.transcriptIDsToInclude.add(line.trim());
		}

		in.close();

		System.out.println("Done- read "+this.transcriptIDsToInclude.size()+" IDs.");
	}


	




	/**
	 * Read the GTF formatted file produced by cuffmerge, adds exon coordinate info to the Sextractor
	 * 
	 * @param gtfFile
	 * @throws IOException
	 */
	public void readGTF(String gtfFile, boolean suppressNs) throws Exception{

		GTF thisGTF = ReadGTF.readGTF(gtfFile, true, this.featuresToKeep, false, suppressNs, this.additionalAttributes);

		HashMap<String, ArrayList<GenomicCoordinate>> coords = thisGTF.getCoordinates();
		ArrayList<GenomicCoordinate> thisChrom;
		GenomicCoordinate tmp;
		Iterator<String> it;
		String thisChromName = "";
		//boolean removeEntry = false;

		// If we have a white-list of transcripts to include:
		if(this.transcriptIDsToInclude.size() > 0){
			//ArrayList<GenomicCoordinate> thisChrom;
			//GenomicCoordinate tmp;
			it = coords.keySet().iterator();
			while(it.hasNext()){
				thisChromName = it.next();
				thisChrom = coords.get(thisChromName);
				Iterator<GenomicCoordinate> it2 = thisChrom.iterator();
				while(it2.hasNext()){
					tmp = it2.next();

					if( ! this.transcriptIDsToInclude.contains(tmp.getAttribute("transcript_id"))){
						coords.get(thisChromName).remove(tmp);
						//this.sextractor.addCoord(tmp);
					}
				}
			}
		}

		this.sextractor.setCoords(coords);


	}

	public void extractSequences(String fastaPath, String outPath) throws IOException{
		this.sextractor.extractFromFasta(fastaPath, true, outPath);
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("annotationPath").hasArg().withDescription("path to the GTF file containing the coordinates to extract").create(Thunder.OPT_PATH_ANNOTATION));
		options.addOption(OptionBuilder.withArgName("fastaPath").hasArg().withDescription("path to the FASTA file containing the sequence from which to extract the coordinates").create(Thunder.OPT_PATH_INPUT));
		options.addOption(OptionBuilder.withArgName("transcriptWhiteList").hasArg().withDescription("text file containing transcript IDs to include").create("f"));
		options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("path to which to output the extracted FASTA sequences").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("featuresToUse").hasArg().withDescription("comma separated list of features to extract from the GTF (e.g. exon,CDS) [default: exon]").create("k"));
		options.addOption(OptionBuilder.withArgName("suppressNs").hasArg().withDescription("true/false remove transcripts with N bases in the sequence [default: true]").create("N"));
		options.addOption(OptionBuilder.withDescription("for EMBOSS translation, if specified will convert stop codons (specified by *) to unknown residues (specified by X)").create("clean"));
		options.addOption(OptionBuilder.withArgName("additionalAttributes").hasArg().withDescription("comma separated list of additional attributes to append to the transcript ID in the output fasta (e.g. transcript_type,transcript_name)").create("t"));
		//		options.addOption(new Option(Thunder.OPT_CHOICE_DB_FORCE_REFRESH, "reads and adds the annotation to the database, regardless of whether the data already exists"));
		return options;
	}



	public static void main(String[] args) throws Exception {

		//args = new String[]{"ExtractSequencesFromGTF","-a","/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v18.annotation.gtf","-i","/Users/robk/WORK/YALE_offline/ANNOTATIONS/hg19.fa","-o","/Users/robk/WORK/YALE_offline/ANNOTATIONS/TEST.fa"};

		//args = new String[]{"ExtractSequencesFromGTF","-a","/Users/robk/WORK/YALE_offline/My Proteomics/HEK/RNA-seq/merged.gtf","-i","/Users/robk/WORK/YALE_offline/ANNOTATIONS/hg19.fa","-o","/Users/robk/WORK/YALE_offline/ANNOTATIONS/TEST.fa"};

		/*args = new String[]{"ExtractSequencesFromGTF",
				"-a","/Users/robk/Box Sync/Work for other people/BaseSpaceApp/TestingOutput/control_vs_comparison.merged.gtf",
				"-i","/Users/robk/WORK/YALE_offline/ANNOTATIONS/hg38.fa",
				"-o","/Users/robk/Box Sync/Work for other people/BaseSpaceApp/TestingOutput/control_vs_comparison.merged.gtf.fasta",
				//"-ID","/Users/robk/Box Sync/Work for other people/BaseSpaceApp/RefSeq2SwissProt/HUMAN_9606_idmapping_selected.filtered.forApp.txt",
				"-ID","/Users/robk/Box Sync/Work for other people/BaseSpaceApp/RefSeq2SwissProt/HUMAN_9606_idmapping_selected.forApp.txt",
				"-NM2NP","/Users/robk/Box Sync/Work for other people/BaseSpaceApp/RefSeq2SwissProt/HomoSapien_NM_NP_conversion_10_16_2014.csv"};
		*/

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption(Thunder.OPT_PATH_ANNOTATION) && cmdArgs.hasOption(Thunder.OPT_PATH_INPUT) && cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT)){

			ExtractSequencesFromGTF thisInstance = new ExtractSequencesFromGTF();
			boolean suppressNs = true;

			if(cmdArgs.hasOption("f")){
				thisInstance.readIncludeList(cmdArgs.getOptionValue("f"));
			}
			if(cmdArgs.hasOption("N")){
				if(cmdArgs.getOptionValue("N").equalsIgnoreCase("false"))
					suppressNs = false;
			}
			if(cmdArgs.hasOption("k")){
				thisInstance.defineFeatureList(cmdArgs.getOptionValue("k"));
			}else{
				thisInstance.defineFeatureList("exon");
			}
			if(cmdArgs.hasOption("t")){
				thisInstance.defineAdditionalAttributes(cmdArgs.getOptionValue("t"));
			}

			thisInstance.readGTF(cmdArgs.getOptionValue(Thunder.OPT_PATH_ANNOTATION), suppressNs);
			thisInstance.extractSequences(cmdArgs.getOptionValue(Thunder.OPT_PATH_INPUT), cmdArgs.getOptionValue(Thunder.OPT_PATH_OUTPUT));

		
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ExtractSequencesFromGTF", getCmdLineOptions());
			System.out.println();
		}

	}

}
