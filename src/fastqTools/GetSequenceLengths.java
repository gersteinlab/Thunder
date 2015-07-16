package fastqTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

import main.Thunder;
import objects.FastX;
import objects.FastX_Record;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class GetSequenceLengths {

	private int maxExpectedSequenceSize = 1000000000;
	private boolean countCollapsedSequences = false;
	
	private HashMap<Integer,Integer> readLengths = new HashMap<Integer,Integer>();


	
	/**
	 * 
	 * @param inputFile
	 * @param countCollapsedSequences
	 * @return
	 * @throws IOException
	 */
	public int[] getReadLengths(File inputFile) throws IOException{

		
		/** Read the lengths of the sequences in the fasta/fastq file **/
		FastX fastx = new FastX(inputFile);
		FastX_Record thisRecord;
		while((thisRecord = fastx.readNext()) != null){
			addReadLength(thisRecord.getID(), thisRecord.getSequence());
		}
		fastx.close();
		

		/** Find the maximum sequence-length **/
		Iterator<Integer> it = this.readLengths.keySet().iterator();
		int maxLength = 0;
		int thisLength;
		while(it.hasNext()){
			thisLength = it.next();
			if(thisLength > maxLength){ maxLength = thisLength; }
		}

		/** Convert HashMap of sequence-lengths to an integer array **/
		int[] readLengthArray = new int[maxLength+1];
		for(int i=0;i<=maxLength;i++){
			if(this.readLengths.containsKey(i)){
				readLengthArray[i] = this.readLengths.get(i); 
			}else{
				readLengthArray[i] = 0;
			}
		}

		return readLengthArray;
	}


	/**
	 * 
	 * @param fastaID
	 * @param sequence
	 */
	private void addReadLength(String fastaID, String sequence){
		int thisSequenceLength = 0;
		int sequenceCount = 1;

		thisSequenceLength = sequence.length();

		if(countCollapsedSequences){
			sequenceCount = Integer.valueOf(fastaID.split("\\s+")[1].trim()).intValue();
		}

		if(thisSequenceLength > maxExpectedSequenceSize){
			System.err.println("sequence too long ("+thisSequenceLength+") -- "+fastaID+" -- "+sequence);
			System.exit(0);
		}

		if(!this.readLengths.containsKey(thisSequenceLength)){
			this.readLengths.put(thisSequenceLength, 0);
		}

		this.readLengths.put(thisSequenceLength, this.readLengths.get(thisSequenceLength)+sequenceCount);
	}



	/**
	 * 
	 * @param readLengths
	 * @return
	 */
	public static String[] readLengthsToString(int[] readLengths){
		String[] readLengthStrings = new String[]{"",""};
		for(int i=0;i<readLengths.length;i++){
			readLengthStrings[0] += i+"\t";
			readLengthStrings[1] += readLengths[i]+"\t";
		}
		return(readLengthStrings);
	}



	/**
	 * 
	 * @param readLengths
	 * @param filePath
	 */
	public static void writeReadLengths(int[] readLengths, String filePath){	
		try{
			String[] readLengthStrings = readLengthsToString(readLengths);
			BufferedWriter out = new BufferedWriter(new FileWriter(new File(filePath)));
			out.write(readLengthStrings[0]+"\n"+readLengthStrings[1]+"\n");
			out.flush();
			out.close();
		} catch (IOException e1) {
			System.err.println("ERROR: Failed to write to specified output file");
			e1.printStackTrace();
		}
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(new Option("countCollapsedSequences", "if FASTA/Q IDs are of the form 'ID 33' (where 33 is the number of identical collapsed sequences) also count these replicate sequences"));
		options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("output sequence lengths to a file [if not specified, lengths are printed to stdout]").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("maxExpectedLength").hasArg().withDescription("print sequences exceeding this maximum expected size (for debugging + QC)").create("m"));
		return options;
	}


	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{

		//args = new String[]{"GetSequenceLengths","-m","51", "/Users/robk/Downloads/mm.circRNA.fa"};
		//args = new String[]{"GetSequenceLengths","-m","75", "/Users/robk/Downloads/test.fq"};
		
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		@SuppressWarnings("unchecked")
		Iterator<String> it = cmdArgs.getArgList().iterator();

		//		System.out.println("cmdArgs.getArgList().size()= "+cmdArgs.getArgList().size());

		if(cmdArgs.getArgList().size() == 2){
			GetSequenceLengths reader = new GetSequenceLengths();

			if(cmdArgs.hasOption("m"))
				reader.maxExpectedSequenceSize = Integer.valueOf(cmdArgs.getOptionValue("m")).intValue();

			reader.countCollapsedSequences = cmdArgs.hasOption("countCollapsedSequences");

			int[] readLengths = new int[0];
			String thisArg = "";
			int count = 0;
			while(it.hasNext()){
				thisArg = it.next();
				if(count == 1){
					//					System.out.println(thisArg);
					readLengths = reader.getReadLengths(new File(thisArg));

					if(cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT)){
						writeReadLengths(readLengths, thisArg);
					}else{
						// If no output file specified, print to stdout:
						String[] tmp = readLengthsToString(readLengths);
						System.out.println(tmp[0]);
						System.out.println(tmp[1]);
					}
				}
				count++;
			}

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" GetSequenceLengths [options] <sequenceFile>", "", getCmdLineOptions(), "");
			System.out.println();
		}
	}
}	