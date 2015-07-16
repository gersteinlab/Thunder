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

public class FastaHeaderGrep {


	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("output sequences shorter than the maximum length to this file [if not specified, sequences are printed to stdout]").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("searchPattern").hasArg().withDescription("pattern to use as a filter in the fasta headers").create("s"));
		return options;
	}




	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{

		//args = new String[]{"FastaHeaderGrep","-s","hsa-", "/Users/robk/Downloads/hairpin.fa"};

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		@SuppressWarnings("unchecked")
		Iterator<String> it = cmdArgs.getArgList().iterator();

		if(cmdArgs.getArgList().size() == 2){

			String searchString = "*"; 
			if(cmdArgs.hasOption("s"))
				searchString = cmdArgs.getOptionValue("s");

			String thisArg = "";
			it.next();
			while(it.hasNext()){
				thisArg = it.next();
				//System.out.println(cmdArgs.getArgList().size());
				
				SequenceReader reader = new SequenceReader(thisArg);
				//SequenceReader reader = new SequenceReader(System.in);
				SequenceRecord tmp;

				if(cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT)){
					BufferedWriter bw = new BufferedWriter(new FileWriter(cmdArgs.getOptionValue(Thunder.OPT_PATH_OUTPUT)));
					while((tmp=reader.readNextRecord()) != null){
						if(tmp.getSequenceID().matches(searchString))
							bw.write(tmp.toString()+"\n");
					}
					bw.flush();
					bw.close();
				}else{
					while((tmp=reader.readNextRecord()) != null){
						//System.out.println(tmp.getSequenceID());
						//if(tmp.getSequenceID().matches(searchString))
						if(tmp.getSequenceID().contains(searchString))
							System.out.println(tmp.toString());
					}
				}
			}
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" FastaHeaderGrep [options] <fastaFile>", "", getCmdLineOptions(), "");
			System.out.println();
		}
	}
}	
