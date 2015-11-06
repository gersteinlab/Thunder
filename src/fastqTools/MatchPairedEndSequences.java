package fastqTools;

import java.io.BufferedOutputStream;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.OutputStreamWriter;
import java.util.HashMap;
import java.util.zip.GZIPOutputStream;

import main.Thunder;
import objects.FastX;
import objects.FastX_Record;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import utils.IO_utils;

public class MatchPairedEndSequences {
	
	private HashMap<String,FastX_Record> read1 = new HashMap<String,FastX_Record>();
	
	private boolean _verbose = false;
	public MatchPairedEndSequences(boolean verbose){
		_verbose = verbose;
	}
	
	/**
	 * 
	 * @param inputFile
	 * @param countCollapsedSequences
	 * @return
	 * @throws IOException
	 */
	public void matchFiles(File inputFile_read1, File inputFile_read2, String outputPrefix, int indexID) throws IOException{
		
		//BufferedWriter out1 = new BufferedWriter(new FileWriter(new File(outputPrefix+"_R1.fastq")));
		//BufferedWriter out2 = new BufferedWriter(new FileWriter(new File(outputPrefix+"_R2.fastq")));
		OutputStreamWriter out1 = new OutputStreamWriter(new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(outputPrefix+"_R1.fq.gz"))));
		OutputStreamWriter out2 = new OutputStreamWriter(new GZIPOutputStream(new BufferedOutputStream(new FileOutputStream(outputPrefix+"_R2.fq.gz"))));
		
		
		int totalReads_r1 = 0, totalReads_r2 = 0, matchedReads = 0;
		
		IO_utils.printLineOut("Reading mate1 reads.");
		FastX fastx = new FastX(inputFile_read1);
		FastX_Record thisRecord;
		while((thisRecord = fastx.readNext()) != null){
			if(_verbose){ 
				System.err.println(thisRecord.getID());
				System.err.println(thisRecord.getID().split(" ")[indexID]); 
			}
			read1.put(thisRecord.getID().split(" ")[indexID], thisRecord);
		}
		fastx.close();
		totalReads_r1 = read1.size();
		//Thunder.printLineOut("Done.");
		
		
		// Read second mates and output matches
		String tmpID;
		fastx = new FastX(inputFile_read2);
		IO_utils.printLineOut("Reading mate2 reads.");
		while((thisRecord = fastx.readNext()) != null){
			tmpID = thisRecord.getID().split(" ")[indexID];
			totalReads_r2 ++;
			if(read1.containsKey(tmpID)){
				// write matched sequences
				out1.write(read1.get(tmpID).toString()+"\n");
				out2.write(thisRecord.toString()+"\n");
				read1.remove(tmpID);		
				matchedReads ++;
			}
		}
		//Thunder.printLineOut("Done ("+read1.size()+" reads).");
		//Thunder.printLineOut("Done.");
		fastx.close();
		
		out1.flush();
		out2.flush();
		out1.close();
		out2.close();
		
		IO_utils.printLineOut("All done.  Matched "+matchedReads+" reads of "+totalReads_r1+" (R1) and "+totalReads_r2+" (R2)");
	}




	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("inputPath_Read1").hasArg().withDescription("File containing the first mate of each pair").create("r1"));
		options.addOption(OptionBuilder.withArgName("inputPath_Read2").hasArg().withDescription("File containing the second mate of each pair").create("r2"));
		options.addOption(OptionBuilder.withArgName("outputPrefix").hasArg().withDescription("Path prefix for output files").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("IDindex").hasArg().withDescription("[Optional] Zero based index fragment (space separated) to use as the ID of each read [default is 0]").create("index"));
		//options.addOption(OptionBuilder.withDescription("Print verbose messages (debugging only)").create("v"));
		return options;
	}


	
	
	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{
	
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		if(cmdArgs.hasOption("r1")  &&  cmdArgs.hasOption("r2")  &&  cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT)){
		
			int index = 0;
			if(cmdArgs.hasOption("index"))
				index = Integer.valueOf(cmdArgs.getOptionValue("index")).intValue();
			
			boolean verbose = false;
			if(cmdArgs.hasOption("v"))
				verbose = true;
			
			MatchPairedEndSequences reader = new MatchPairedEndSequences(verbose);
			reader.matchFiles(new File(cmdArgs.getOptionValue("r1")), new File(cmdArgs.getOptionValue("r2")), cmdArgs.getOptionValue(Thunder.OPT_PATH_OUTPUT), index);

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" MatchPairedEndSequences -r1 <sequenceFile_mate1> -r2 <sequenceFile_mate2> -o <outputPrefix>", "", getCmdLineOptions(), "");
			System.out.println();
		}
	}
}	