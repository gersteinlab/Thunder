package fastqTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;

import main.Thunder;
import objects.FastX;
import objects.FastX_Record;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.lang3.StringUtils;

public class RemoveHomopolymers {

	
	public static void processFile(String inputPath, String outputPath, double maxFrac, boolean rejected2stderr) throws IOException{
		
		BufferedWriter out = new BufferedWriter(new FileWriter(new File(outputPath)));
		int a,c,g,t;
		int sequencesRemoved = 0;
		
		FastX fastx = new FastX(new File(inputPath));
		FastX_Record thisRecord;
		while((thisRecord = fastx.readNext()) != null){
			int sum = 0;
			int max = 0;
			
			a = StringUtils.countMatches(thisRecord.getSequence().toUpperCase(), "A");
			sum += a;
			max = a;
			
			c = StringUtils.countMatches(thisRecord.getSequence().toUpperCase(), "C");
			sum += c;
			if(c > max)
				max = c;
			
			g = StringUtils.countMatches(thisRecord.getSequence().toUpperCase(), "G");
			sum += g;
			if(g > max)
				max = g;
			
			t = StringUtils.countMatches(thisRecord.getSequence().toUpperCase(), "T");
			sum += t;
			if(t > max)
				max = t;
			
			//System.out.println("max="+max+"\tsum="+sum);
			
			if((max+0.0)/(sum+0.0) < maxFrac)
				out.write(thisRecord.toString()+"\n");
			else{
				sequencesRemoved ++;
				if(rejected2stderr)
					System.err.println(thisRecord.toString());
			}
			
		}
		
		System.out.println("\nDone.  Sequences removed="+sequencesRemoved);
		fastx.close();
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
		options.addOption(new Option("verbose", "write the rejected sequences to stderr"));
		options.addOption(OptionBuilder.withArgName("inputPath").hasArg().withDescription("input fasta/fastq").create(Thunder.OPT_PATH_INPUT));
		options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("output sequence lengths to a file [if not specified, lengths are printed to stdout]").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("maxAcceptableNucleotideFrequency").hasArg().withDescription("output only sequences with fewer than this fraction of single nucleotides [default: 0.7]").create("m"));
		return options;
	}
	
	
	public static void main(String[] args) throws Exception {
		//args = new String[]{"RemoveHomopolymerRepeats","-m","0.7", "-"+Thunder.OPT_PATH_INPUT,"/Users/robk/Downloads/test.fq", "-"+Thunder.OPT_PATH_OUTPUT,"/Users/robk/Downloads/test.fq.out"};
		//args = new String[]{"RemoveHomopolymerRepeats","-m","0.7", "-"+Thunder.OPT_PATH_INPUT,"/Users/robk/Downloads/test.fq", "-"+Thunder.OPT_PATH_OUTPUT,"/Users/robk/Downloads/test.fq.out"};
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		//Iterator<String> it = cmdArgs.getArgList().iterator();
		
		
		
		if(cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT)  &&  cmdArgs.hasOption(Thunder.OPT_PATH_INPUT)){
			
			double maxFrac = 0.7;
			if(cmdArgs.hasOption("m"))
				maxFrac = Double.valueOf(cmdArgs.getOptionValue("m")).doubleValue();
			
			boolean verbose = false;
			if(cmdArgs.hasOption("verbose"))
				verbose = true;
			
			processFile(cmdArgs.getOptionValue(Thunder.OPT_PATH_INPUT), cmdArgs.getOptionValue(Thunder.OPT_PATH_OUTPUT), maxFrac, verbose);
			
		}else{
			HelpFormatter formatter = new HelpFormatter();
			//formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" RemoveHomopolymerRepeats", getCmdLineOptions());
			formatter.printHelp("RemoveHomopolymerRepeats", getCmdLineOptions());
			System.out.println();
		}


	}

}
