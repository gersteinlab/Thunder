package fastaTools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Iterator;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import utils.IO_utils;
import fastqTools.SequenceReader;
import fastqTools.SequenceRecord;

public class FastaRemoveDuplicates {

	private ArrayList<String> _headers = new ArrayList<String>();
	private ArrayList<String> _sequences = new ArrayList<String>();
	
	/**
	 * 
	 * @param inputFasta
	 * @param outputFasta
	 * @throws IOException
	 */
	public FastaRemoveDuplicates(String inputFasta, String outputFasta, boolean alsoFilterDupSequences) throws IOException{
		SequenceReader reader = new SequenceReader(inputFasta);
		SequenceRecord tmp;
		int counter = 0;
		int removedCount = 0;
		
		if(outputFasta != null){
			BufferedWriter bw = new BufferedWriter(new FileWriter(outputFasta));
			while((tmp=reader.readNextRecord()) != null){
				counter ++;
				if(!_headers.contains(tmp.getSequenceID())){
					
					if(!alsoFilterDupSequences){
						// if we only want to filter by fasta header, output this record
						_headers.add(tmp.getSequenceID());
						bw.write(tmp.toString()+"\n");
					}else{
						// if we also want to remove duplicate sequences, check now
						if(!_sequences.contains(tmp.getSequence())){
							_sequences.add(tmp.getSequence());
							_headers.add(tmp.getSequenceID());
							bw.write(tmp.toString()+"\n");
						}else{
							removedCount ++;	
						}
					}
				}else{
					removedCount ++;
				}
			}
			bw.flush();
			bw.close();
		}else{
			while((tmp=reader.readNextRecord()) != null){
				counter ++;
				if(!_headers.contains(tmp.getSequenceID())){
					if(!alsoFilterDupSequences){
						// if we only want to filter by fasta header, output this record
						_headers.add(tmp.getSequenceID());
						System.out.println(tmp.toString());
					}else{
						// if we also want to remove duplicate sequences, check now
						if(!_sequences.contains(tmp.getSequence())){
							_sequences.add(tmp.getSequence());
							_headers.add(tmp.getSequenceID());
							System.out.println(tmp.toString());
						}else{
							removedCount ++;	
						}
					}
				}else{
					removedCount ++;
				}
			}
		}
		
		IO_utils.printLineErr("Done.  Read "+counter+" fasta entries, removed "+removedCount+" duplicates, output "+_headers.size()+" unique entries.");
	}
	
	
	
	
	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("output unique fasta entries to this file [if not specified, sequences are printed to stdout]").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("boolean").withDescription("also remove duplicate sequences (without this only duplicate headers will be removed)").create("s"));
		return options;
	}


	
	 


	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{

		//args = new String[]{"FastaRemoveDuplicates", "/Users/robk/Downloads/hairpin.fa"};

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		@SuppressWarnings("unchecked")
		Iterator<String> it = cmdArgs.getArgList().iterator();

		if(cmdArgs.getArgList().size() == 2){

			boolean alsoFilterDupSequences = false;
			if(cmdArgs.hasOption("s")){
				alsoFilterDupSequences = true;
			}
			
			String thisArg = "";
			it.next();
			while(it.hasNext()){
				thisArg = it.next();

				if(cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT))
					new FastaRemoveDuplicates(thisArg, cmdArgs.getOptionValue(Thunder.OPT_PATH_OUTPUT), alsoFilterDupSequences);
				else
					new FastaRemoveDuplicates(thisArg, null, alsoFilterDupSequences);
			}
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" FastaRemoveDuplicates [options] <fastaFile>", "", getCmdLineOptions(), "");
			System.out.println();
		}
	}
}	
