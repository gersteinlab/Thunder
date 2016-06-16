package fastqTools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import utils.IO_utils;

public class FilterFastxByHeaderList {

	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		//options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("output sequences shorter than the maximum length to this file [if not specified, sequences are printed to stdout]").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("IDListPath").hasArg().withDescription("file containing sequence IDs to use as a filter in the fasta headers (one ID per line)").create("IDs"));
		options.addOption(OptionBuilder.withArgName("whitelist").withDescription("[default] ID list is a whitelist of fasta headers to INCLUDE").create("w"));
		options.addOption(OptionBuilder.withArgName("blacklist").withDescription("[optional] ID list is a blacklist of fasta headers to EXCLUDE.  If both -w and -b are specified, -b takes preference").create("b"));
		options.addOption(OptionBuilder.withArgName("allowIDPrefixes").withDescription("[optional] ID can be a PREFIX of the fullID in the fasta/q header").create("p"));
		options.addOption(OptionBuilder.withArgName("allowIDSuffixes").withDescription("[optional] ID can be a SUFFIX of the fullID in the fasta/q header").create("s"));
		return options;
	}


	/**
	 * Read the ID list from a text file (one ID per line)
	 * @param listPath
	 * @return
	 * @throws IOException
	 */
	public static HashSet<String> readList(File listPath) throws IOException{
		IO_utils.printLineErr("Reading header list");
		HashSet<String> idList = new HashSet<String>();
		BufferedReader in = new BufferedReader(new FileReader(listPath));
		String line = "";
		while((line=in.readLine()) != null){
			idList.add(line.trim());
		}
		in.close();
		IO_utils.printLineErr("Read "+idList.size()+" headers");
		return idList;
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{

		//args = new String[]{"FilterFastxByIDList","-w","-p","-IDs","/Users/robk/Downloads/readsMappedToLibs.txt","/Users/robk/Downloads/endogenousUnaligned_ungapped.fq"};

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		@SuppressWarnings("unchecked")
		Iterator<String> it = cmdArgs.getArgList().iterator();

		if(cmdArgs.getArgList().size() >= 2){

			// set whitelist or blacklist
			boolean whitelist = true; 
			if(cmdArgs.hasOption("b"))
				whitelist = false;

			// read the ID list
			HashSet<String> idList = new HashSet<String>();
			if(cmdArgs.hasOption("IDs"))
				idList = readList(new File(cmdArgs.getOptionValue("IDs")));
			
			
			String thisArg = "";
			it.next();
			while(it.hasNext()){
				int countKept = 0;
				int countExcluded = 0;
				thisArg = it.next();
				IO_utils.printLineErr("Reading Fasta/q file:"+thisArg);
				//System.out.println(cmdArgs.getArgList().size());

				SequenceReader reader = new SequenceReader(thisArg);
				SequenceRecord tmp;
				while((tmp=reader.readNextRecord()) != null){

					// decide whether this fasta/q ID is in the list
					boolean matchesList = false;
					if(cmdArgs.hasOption("p")){
						matchesList = idList.contains(tmp.getSequenceID().split(" ")[0]);
					}else if(cmdArgs.hasOption("s")){
						String[] tmpHeader = tmp.getSequenceID().split(" ");
						matchesList = idList.contains(tmpHeader[tmpHeader.length-1]);
					}else{
						matchesList = idList.contains(tmp.getSequenceID());
					}

					if(matchesList){
						// this fasta/q header is in the list
						if(whitelist){
							System.out.println(tmp.toString());
							countKept ++;
						}else{
							countExcluded ++;
						}
					}else{
						// this fasta/q header is not in the list
						if(!whitelist){
							System.out.println(tmp.toString());
							countKept ++;
						}else{
							countExcluded ++;
						}
					}
				}
				IO_utils.printLineErr("\tReads kept: "+countKept);
				IO_utils.printLineErr("\tReads excluded: "+countExcluded);
			}
			IO_utils.printLineErr("Done.");
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" FilterFastxByHeaderList [options] <fasta/q File> [fasta/q File]", "", getCmdLineOptions(), "");
			System.out.println();
		}
	}
}	
