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

public class Fasta_U2T {


	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("output fastq sequences to this file [if not specified, sequences are printed to stdout]").create(Thunder.OPT_PATH_OUTPUT));
		return options;
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		@SuppressWarnings("unchecked")
		Iterator<String> it = cmdArgs.getArgList().iterator();

		if(cmdArgs.getArgList().size() == 2){

			String thisArg = "";
			int count = 0;
			while(it.hasNext()){
				thisArg = it.next();
				if(count == 1){
					SequenceReader reader = new SequenceReader(thisArg);
					//SequenceReader reader = new SequenceReader(System.in);
					SequenceRecord tmp, tmp2;

					if(cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT)){
						BufferedWriter bw = new BufferedWriter(new FileWriter(cmdArgs.getOptionValue(Thunder.OPT_PATH_OUTPUT)));
						
						while((tmp=reader.readNextRecord()) != null){
							tmp2 = new SequenceRecord(tmp.getSequenceID());
							tmp2.addSequenceString(tmp.getSequence().replace("U", "T"));
							bw.write(tmp2.toString()+"\n");
						}
						bw.flush();
						bw.close();
					}else{
						while((tmp=reader.readNextRecord()) != null){
							tmp2 = new SequenceRecord(tmp.getSequenceID());
							tmp2.addSequenceString(tmp.getSequence().replace("U", "T"));
							System.out.println(tmp2.toString());
						}
					}
				}
				count++;
			}
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" Fasta_U2T [options] <fastaFile>", "", getCmdLineOptions(), "");
			System.out.println();
		}
	}
}	
