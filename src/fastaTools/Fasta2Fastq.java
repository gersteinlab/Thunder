package fastaTools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.util.Iterator;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import fastqTools.SequenceReader;
import fastqTools.SequenceRecord;

public class Fasta2Fastq {


	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("output fastq sequences to this file [if not specified, sequences are printed to stdout]").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("insertQchar").hasArg().withDescription("dummy quality value (for Phred+33 best to use \'I\', else 'g' for Phred+64)) [default: \'I\' ]").create("q"));
		return options;
	}




	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{

		//args = new String[]{"GetSequenceLengths","-m","100", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/test.fa"};

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		@SuppressWarnings("unchecked")
		Iterator<String> it = cmdArgs.getArgList().iterator();

		if(cmdArgs.getArgList().size() == 2){

			// Override default max sequence length if set
			// NB: The longest protein in SwissProt is TITIN_MOUSE at 35213 amino-acids
			String qualValue = "I"; // this is the maximum value that X!Tandem will accept before segfaulting
			if(cmdArgs.hasOption("q"))
				qualValue = cmdArgs.getOptionValue("q").substring(0, 1);

			String thisArg = "";
			int count = 0;
			while(it.hasNext()){
				thisArg = it.next();
				if(count == 1){
					SequenceReader reader = new SequenceReader(thisArg);
					//SequenceReader reader = new SequenceReader(System.in);
					SequenceRecord tmp;

					if(cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT)){
						BufferedWriter bw = new BufferedWriter(new FileWriter(cmdArgs.getOptionValue(Thunder.OPT_PATH_OUTPUT)));
						while((tmp=reader.readNextRecord()) != null){
							tmp.setQuality(new String(new char[tmp.getSequenceLength()]).replace("\0", qualValue));
							bw.write(tmp.toString()+"\n");
						}
						bw.flush();
						bw.close();
					}else{
						while((tmp=reader.readNextRecord()) != null){
							tmp.setQuality(new String(new char[tmp.getSequenceLength()]).replace("\0", qualValue));
							System.out.println(tmp.toString());
						}
					}
				}
				count++;
			}
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" Fasta2Fastq [options] <sequenceFile>", "", getCmdLineOptions(), "");
			System.out.println();
		}
	}
}	
