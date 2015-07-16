package transcriptome;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;
import java.io.Writer;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Iterator;

import main.Thunder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class DEPRECIATED_RibosomeFootprintEM {
	public static boolean verbose = false;

	public String recordSummary(SAMRecord record){
		return record.getReadName()+"\t"+record.getReferenceName()+"\t"+record.getIntegerAttribute("NH")+"\t"+record.getIntegerAttribute("HI");
	}

	public void readSAM(File input, File output) throws Exception{
		System.out.println(this.getTime()+" Reading alignments: "+input.getAbsolutePath());
		SAMFileReader inputSam = new SAMFileReader(input);
		SAMRecord thisRecord;
		SAMRecordIterator it = inputSam.iterator();
		//int count = 0;

		Writer out = new BufferedWriter(new FileWriter(output));
		out.write("readID\ttranscriptID\n");

		//System.out.println(inputSam.getFileHeader().getSortOrder().toString());
		//if(inputSam.getFileHeader().getSortOrder().compareTo(SortOrder.queryname) == 0){
		//	System.out.println("Sorted correctly!");
		//}
		//it.assertSorted(SortOrder.queryname);
		StringBuffer sb;
		String thisGeneID = "";
		String thisReadID = "";
		boolean multimapper = false;
		int n = 0;
		ArrayList<String> tmp_transcripts;

		while(it.hasNext()){
			thisRecord = it.next();
			sb = new StringBuffer();
			tmp_transcripts = new ArrayList<String>();

			if(thisRecord.getIntegerAttribute("HI") == 1){
				thisReadID = thisRecord.getReadName();

				sb.append(thisRecord.getReadName()+"\t"+thisRecord.getReferenceName()+"\n");

				tmp_transcripts.add(thisRecord.getReferenceName());
				multimapper = false;

				n = thisRecord.getIntegerAttribute("NH");

				if(transcriptID_2_geneID.containsKey(thisRecord.getReferenceName())){
					thisGeneID = transcriptID_2_geneID.get(thisRecord.getReferenceName());

					//System.out.println("i=1: "+recordSummary(thisRecord)+"\t"+thisGeneID);

					for(int i=2;i<=n;i++){
						thisRecord = it.next();
						//System.out.print("i="+i+": "+recordSummary(thisRecord));

						if(thisRecord.getIntegerAttribute("HI") == i  &  thisRecord.getReadName().equals(thisReadID)){

							if(!tmp_transcripts.contains(thisRecord.getReferenceName())){
								tmp_transcripts.add(thisRecord.getReferenceName());
								sb.append(thisRecord.getReadName()+"\t"+thisRecord.getReferenceName()+"\n");
							}
							//System.out.println("\t"+transcriptID_2_geneID.get(thisRecord.getReferenceName()));

							if(!thisGeneID.equals(transcriptID_2_geneID.get(thisRecord.getReferenceName())))
								multimapper = true;
						}else{
							System.err.println("Error 2: This BAM needs to be sorted by readID!");
						}
					}
				}
			}else{
				System.err.println("Error 1: This BAM needs to be sorted by readID!");
			}

			//out.write(thisRecord.getReadName()+"\t"+thisRecord.getReferenceName()+"\n");
			//System.out.println("\tmultimapper="+multimapper);

			if(!multimapper)
				out.write(sb.toString());

			//if(thisRecord.getIntegerAttribute("NH").equals(1)){
			//if(verbose) System.out.print(thisRecord.getReferenceName()+":"+thisRecord.getAlignmentStart()+"-"+thisRecord.getAlignmentEnd()+"   "+thisRecord.getReadString());
			//if(verbose) System.out.println("   NM:"+thisRecord.getIntegerAttribute("NM")+"   NH:"+thisRecord.getIntegerAttribute("NH")+"   "+thisRecord.getCigarString());
			//}
			//			count++;
			//			if(count > 20)
			//				break;
		}

		out.flush();
		out.close();
		inputSam.close();
		System.out.println(this.getTime()+" Done.");
	}






	public void readGTF(File input) throws Exception{
		System.out.println(this.getTime()+" Reading GTF: "+input.getAbsolutePath());
		BufferedReader in = new BufferedReader(new FileReader(input));
		String line = "";

		while((line=in.readLine()) != null){
			if(!line.startsWith("##")){
				parseLine(line);
			}
		}
		in.close();
		System.out.println(this.getTime()+" Done.");
	}

	private HashMap<String, ArrayList<String>> gene2Transcripts = new HashMap<String, ArrayList<String>>();
	private HashMap<String, String> transcriptID_2_geneID = new HashMap<String, String>();
	private HashMap<String, String> transcriptID_2_transcriptType = new HashMap<String, String>();
	private HashMap<String, Integer> transcriptID_2_transcriptLength = new HashMap<String, Integer>();
	private HashMap<String, Integer> transcriptID_2_CDSLength = new HashMap<String, Integer>();

	private void parseLine(String line){
		String[] bits = line.split(" |\t");;

		if(bits[2].equals("exon")){
			String transcriptID = trimAttribute(bits[11].trim());
			addTranscript2Gene(transcriptID, trimAttribute(bits[9].trim()));
			addTranscript2Type(transcriptID, trimAttribute(bits[19].trim()));
			addExonLength2Transcript(transcriptID, (Integer.valueOf(bits[4])-Integer.valueOf(bits[3])));
		}else if(bits[2].equals("CDS")){
			String transcriptID = trimAttribute(bits[11].trim());
			addTranscript2Gene(transcriptID, trimAttribute(bits[9].trim()));
			addTranscript2Type(transcriptID, trimAttribute(bits[19].trim()));
			addCDSLength2Transcript(transcriptID, (Integer.valueOf(bits[4])-Integer.valueOf(bits[3])));
		}
	}
	public static String trimAttribute(String in){
		return in.replaceAll("^\"|\";$", "");
	}
	private void addTranscript2Gene(String transcriptID, String geneID){
		if(!transcriptID_2_geneID.containsKey(transcriptID))
			transcriptID_2_geneID.put(transcriptID, geneID);

		if(!gene2Transcripts.containsKey(geneID))
			gene2Transcripts.put(geneID, new ArrayList<String>());

		if(!gene2Transcripts.get(geneID).contains(transcriptID))
			gene2Transcripts.get(geneID).add(transcriptID);
	}
	private void addTranscript2Type(String transcriptID, String transcriptType){
		if(!transcriptID_2_transcriptType.containsKey(transcriptID))
			transcriptID_2_transcriptType.put(transcriptID, transcriptType);
	}
	private void addExonLength2Transcript(String transcriptID, int length){
		if(!transcriptID_2_transcriptLength.containsKey(transcriptID)){
			transcriptID_2_transcriptLength.put(transcriptID, length);
		}else{
			transcriptID_2_transcriptLength.put(transcriptID, transcriptID_2_transcriptLength.get(transcriptID) + length);
		}
	}
	private void addCDSLength2Transcript(String transcriptID, int length){
		if(!transcriptID_2_CDSLength.containsKey(transcriptID)){
			transcriptID_2_CDSLength.put(transcriptID, length);
		}else{
			transcriptID_2_CDSLength.put(transcriptID, transcriptID_2_CDSLength.get(transcriptID) + length);
		}
	}



	private HashMap<String, Double> transcriptQuants = new HashMap<String, Double>();

	/**
	 * 
	 * @param inputFile
	 * @throws Exception
	 */
	public void readTranscriptExpressions(File inputFile) throws IOException{
		System.out.println(this.getTime()+" Reading eXpress expression data: "+inputFile.getAbsolutePath());
		BufferedReader in = new BufferedReader(new FileReader(inputFile));
		String line = "";
		String[] bits;

		// burn the first line
		in.readLine();

		while((line=in.readLine()) != null){
			bits = line.split("\t");

			transcriptQuants.put(bits[1].trim(), Double.valueOf(bits[14].trim()).doubleValue());
		}
		in.close();
		System.out.println(this.getTime()+" Done.");
	}


	/**
	 * 
	 */
	public void readTranscriptExpressions(){
		System.out.println(this.getTime()+" Initialising null priors...");

		Iterator<String> iterator = transcriptID_2_geneID.keySet().iterator();
		String transcriptID = "";

		while(iterator.hasNext()){
			transcriptID = iterator.next();
			transcriptQuants.put(transcriptID, 1.0/(gene2Transcripts.get(transcriptID_2_geneID.get(transcriptID)).size()+0.0));
		}

		System.out.println(this.getTime()+" Done.");
	}

	/**
	 * 
	 * @param outputFile
	 * @throws IOException
	 */
	public void transcriptExpression2prior(File outputFile) throws IOException{
		Writer out = new BufferedWriter(new FileWriter(outputFile));
		out.write("geneID\ttranscriptID\tlength\tprior\n");

		Iterator<String> iterator_gene = gene2Transcripts.keySet().iterator();
		ArrayList<String> transcripts;
		Iterator<String> iterator_transcript;
		double sum = 0.0;
		double frac = 0.0;
		String geneID = "";
		String transcriptID = "";
		int length = 0;
		while(iterator_gene.hasNext()){
			geneID = iterator_gene.next();


			transcripts = gene2Transcripts.get(geneID);
			iterator_transcript = transcripts.iterator(); 

			// Calculate the total TPM for this gene
			sum = 0.0;
			while(iterator_transcript.hasNext()){
				transcriptID = iterator_transcript.next();
				if(transcriptQuants.containsKey(transcriptID))
					sum += transcriptQuants.get(transcriptID);
			}

			// Convert the transcript TPM to a fraction of the gene-level totalTPM
			iterator_transcript = transcripts.iterator();
			while(iterator_transcript.hasNext()){
				transcriptID = iterator_transcript.next();

				length = 0;
				if(transcriptID_2_transcriptLength.containsKey(transcriptID))
					length = transcriptID_2_transcriptLength.get(transcriptID);

				frac = 0;
				if(transcriptQuants.containsKey(transcriptID)  &  sum > 0){
					frac = transcriptQuants.get(transcriptID) / sum;
				}

				if(sum == 0.0){
					frac = 1.0/(transcripts.size())+0.0;
				}

				out.write(geneID+"\t"+transcriptID+"\t"+length+"\t"+frac+"\n");	
			}
		}

		out.flush();
		out.close();
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("annotationPath").hasArg().withDescription("Path to the GTF file containing the coordinates to extract").create(Thunder.OPT_PATH_ANNOTATION));
		options.addOption(OptionBuilder.withArgName("results.xprs").hasArg().withDescription("[optional] Text file containing RNA-seq transcript expression quantifications from eXpress").create("e"));
		options.addOption(OptionBuilder.withArgName("footprint alignments").hasArg().withDescription("SAM/BAM file containing alignments").create("f"));
		options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("Path to which to output results").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("maxIterations_EM").hasArg().withDescription("Maximum number of EM iterations [default: 500]").create("N"));
		options.addOption(OptionBuilder.withDescription("Keep intermediate files").create("DEV"));
		return options;
	}


	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {

		//		args = new String[]{"ParseFootprintAlignments", 
		//				"-f", "/Users/robk/Box Sync/Work/HEK293_RNAseq/STAR_eXpress/D1-Footprint_1x50_gencode.v18.annotation_mappedReads.bam", 
		//				"-a", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v18.annotation.filtered.gtf",
		//				//"-e", "/Users/robk/Box Sync/Work/HEK293_RNAseq/STAR_eXpress/A1-Total_2x75/results.xprs",
		//				"-o", "/Users/robk/Box Sync/Work/HEK293_RNAseq/EM_isoformQuant/firstRound_fullFootprintEM_fromJing/FootPrint/data/TEST_NEW"};
		//		
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption(Thunder.OPT_PATH_ANNOTATION) && cmdArgs.hasOption("f") && cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT)){
			System.out.println();
			DEPRECIATED_RibosomeFootprintEM engine = new DEPRECIATED_RibosomeFootprintEM();

			String output_path = cmdArgs.getOptionValue("o");
			String tmp_path = output_path+"/tmp";
			File output_priors = new File(output_path+"/priors.txt");
			File output_alignments = new File(tmp_path+"/alignments.txt");


			// Create the output directory and temp dir if they don't already exist
			engine.makeOutputDirectories(output_path);

			// Read the transcripts in the GTF
			engine.readGTF(new File(cmdArgs.getOptionValue(Thunder.OPT_PATH_ANNOTATION)));

			// Read the eXpress transcript quantifications or choose flat priors
			if(cmdArgs.hasOption("e")){
				engine.readTranscriptExpressions(new File(cmdArgs.getOptionValue("e")));
			}else{
				engine.readTranscriptExpressions();
			}

			// Convert the transcript expressions to priors and output
			engine.transcriptExpression2prior(output_priors);

			// Parse the footprint alignments
			engine.readSAM(new File(cmdArgs.getOptionValue("f")), output_alignments);

			// Make temporary copies of the perl EM executables:
			engine.makePerlExecutables(tmp_path);

			// set the number of iterations to cap the EM
			int maxIterations = 500;
			if(cmdArgs.hasOption("N"))
				maxIterations = Integer.valueOf(cmdArgs.getOptionValue("N")).intValue();
			
			// Run the EM!
			System.out.println(engine.getTime()+" Running EM (maxIterations: "+maxIterations+") -- this might take a few minutes...");
			engine.runEM("perl -I "+tmp_path+" "+tmp_path+"/main.foot.print.pl "+maxIterations+" "+output_priors+" "+output_alignments+" "+output_path);

			// Clean up
			if(!cmdArgs.hasOption("DEV"))
				engine.removeFile(new File(tmp_path));

			System.out.println(engine.getTime()+" All Done!");
			
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ParseFootprintAlignments", getCmdLineOptions());
			System.out.println();
		}


	}
	
	public String getTime(){
		return((new SimpleDateFormat("yyyy/MM/dd HH:mm:ss")).format(Calendar.getInstance().getTime()));
	}
	

	public void runEM(String perlCommand) throws IOException, InterruptedException{
		Process p = Runtime.getRuntime().exec(perlCommand);
		p.waitFor();
		BufferedReader reader_stdout = new BufferedReader(new InputStreamReader(p.getInputStream()));
		BufferedReader reader_stderr = new BufferedReader(new InputStreamReader(p.getErrorStream()));

		String line = "";	
		StringBuffer sb_stdout = new StringBuffer();
		StringBuffer sb_stderr = new StringBuffer();
		while ((line = reader_stdout.readLine())!= null) 
			sb_stdout.append(line + "\n");
		while ((line = reader_stderr.readLine())!= null) 
			sb_stderr.append(line + "\n");

		System.err.println(sb_stderr.toString());
		//System.out.println(sb_stdout.toString());
	}

	public void makePerlExecutables(String output_path) throws IOException{
		makePerlExecutable(output_path, "main.foot.print.pl");
		makePerlExecutable(output_path, "footPrint_EM.pm");
		makePerlExecutable(output_path, "proprocessing_footPrint_EM.pm");
	}

	private void makePerlExecutable(String output_path, String exeName) throws IOException{
		InputStream is = this.getClass().getResourceAsStream("/EM/"+exeName);
		if (null != is) {
			BufferedReader reader = new BufferedReader(new InputStreamReader(is, "UTF-8"));

			Writer out = new BufferedWriter(new FileWriter(output_path+"/"+exeName));
			String line = "";

			while((line=reader.readLine()) != null)
				out.write(line+"\n");

			out.flush();
			out.close();
			reader.close();
		}
	}



	public void makeOutputDirectories(String output_path){
		File tmp = new File(output_path);
		if(!tmp.exists())
			tmp.mkdir();
		tmp = new File(output_path+"/tmp");
		if(!tmp.exists())
			tmp.mkdir();
	}

	public void removeFile(File file){
		if(file.isDirectory()){
			File[] children = file.listFiles();
			for(int i=0;i<children.length;i++)
				removeFile(children[i]);
		}
		file.delete();
		//System.out.println("Deleted: "+file.getAbsolutePath());
	}

}
