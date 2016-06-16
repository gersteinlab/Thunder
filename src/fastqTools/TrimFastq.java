package fastqTools;

import java.io.BufferedWriter;
import java.io.FileWriter;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class TrimFastq {



	/**
	 * 
	 * @param in
	 * @return
	 */
	public static SequenceRecord processRead(SequenceRecord in, int trim5p, int trim3p){
		//System.err.println("\n\n"+in.toString());
		SequenceRecord out = in;
		String seq = in.getSequence();
		String qual = in.getQuality();

		if(trim5p+trim3p >= seq.length()){
			out = null;
		}else{
			if(trim5p > 0){
				seq = seq.substring(trim5p);
				qual = qual.substring(trim5p);
			}
			if(trim3p > 0){
				seq = seq.substring(0, seq.length()-trim3p);
				qual = qual.substring(0, qual.length()-trim3p);
			}
			out = new SequenceRecord(in.getSequenceID());
			out.addSequenceString(seq);
			out.addQualityString(qual);
		}

		//System.err.println(out.toString());
		return out;
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("output sequences shorter than the maximum length to this file [if not specified, sequences are printed to stdout]").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("trim [int] bases from the 5' end of each read [default:0]").create("5p"));
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("trim [int] bases from the 3' end of each read [default:0]").create("3p"));
		options.addOption(OptionBuilder.withArgName("path").hasArg().withDescription("read input fastq at this path [default:stdin]").create("i"));
		options.addOption(OptionBuilder.withArgName("path").hasArg().withDescription("write input fastq to this path [default:stdout]").create("o"));
		return options;
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{

		//args = new String[]{"TrimFastq","-5p","3", "-3p","1", "-i","/Users/robk/Downloads/test.noRand.fq"};

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption("5p")  ||  cmdArgs.hasOption("3p")){

			// What to trim:
			int trim5p=0, trim3p=0;
			if(cmdArgs.hasOption("5p"))
				trim5p = Integer.valueOf(cmdArgs.getOptionValue("5p")).intValue();
			if(cmdArgs.hasOption("3p"))
				trim3p = Integer.valueOf(cmdArgs.getOptionValue("3p")).intValue();


			// Where to read from:
			SequenceReader reader;
			if(cmdArgs.hasOption("i"))
				reader = new SequenceReader(cmdArgs.getOptionValue("i"));
			else
				reader = new SequenceReader(System.in);


			// Continue?
			SequenceRecord tmpIn, tmpOut;
			//if(trim5p > 0  ||  trim3p > 0){
			// Where to write to:
			if(cmdArgs.hasOption("o")){
				BufferedWriter bw = new BufferedWriter(new FileWriter(cmdArgs.getOptionValue("o")));
				while((tmpIn=reader.readNextRecord()) != null){
					tmpOut = processRead(tmpIn, trim5p, trim3p);
					if(tmpOut != null)
						bw.write(tmpOut.toString()+"\n");
				}
				bw.flush();
				bw.close();
			}else{
				while((tmpIn=reader.readNextRecord()) != null){
					tmpOut = processRead(tmpIn, trim5p, trim3p);
					if(tmpOut != null)
						System.out.println(tmpOut.toString());
				}
			}
			reader.close();
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" TrimFastq [options] <sequenceFile>", "", getCmdLineOptions(), "");
			System.out.println();
		}

		
	}
}	
