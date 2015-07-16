package transcriptome;

import java.io.File;
import java.util.Iterator;

import main.Thunder;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class CIGAR_2_PWM {

	public static void readSAM(File input, int maxRecords) throws Exception{
		System.err.println("Reading alignments: "+input.getAbsolutePath());
		SAMFileReader inputSam = new SAMFileReader(input);
		SAMRecord thisRecord;
		SAMRecordIterator samIterator = inputSam.iterator();
		
		int records = 0;
		int count = 0;

		int[] m = new int[75];
		int[] s = new int[75];
		int[] i = new int[75];
		
		while(samIterator.hasNext()){
			
			records ++;
			count = 0;
			thisRecord = samIterator.next();
			
			Iterator<CigarElement> cigarElements = thisRecord.getCigar().getCigarElements().iterator();
			CigarElement tmp;
			CigarOperator tmp_op;
			
			while(cigarElements.hasNext()){
				tmp = cigarElements.next();
				tmp_op = tmp.getOperator();
				
				if(!tmp_op.equals(CigarOperator.N)  &  !tmp_op.equals(CigarOperator.D)){
					// ignore deletions to the reference (introns?)
					
					if(tmp_op.equals(CigarOperator.M)){
						for(int x=0;x<tmp.getLength();x++)
							m[count+x] ++;
						count += tmp.getLength();
					}else if(tmp_op.equals(CigarOperator.S)){
						for(int x=0;x<tmp.getLength();x++)
							s[count+x] ++;
						count += tmp.getLength();
					}else if(tmp_op.equals(CigarOperator.I)){
						for(int x=0;x<tmp.getLength();x++)
							i[count+x] ++;
						count += tmp.getLength();
					}
					else{
						System.err.println("Unknown CIGAR operator: "+tmp_op.toString()+" in "+thisRecord.getCigarString());
					}
				}
			}
			
			if(maxRecords > 0  &  records >= maxRecords)
				break;
			
		}
			
		System.err.println("SAM records read = "+records);
		
		System.out.print("M\t");
		for(int x=0;x<m.length-1;x++)
			System.out.print((m[x]+0.0)/records+"\t");
		System.out.println((m[m.length-1]+0.0)/records);
		
		System.out.print("S\t");
		for(int x=0;x<s.length-1;x++)
			System.out.print((s[x]+0.0)/records+"\t");
		System.out.println((s[s.length-1]+0.0)/records);
		
		System.out.print("I\t");
		for(int x=0;x<i.length-1;x++)
			System.out.print((i[x]+0.0)/records+"\t");
		System.out.println((i[i.length-1]+0.0)/records);
		
		inputSam.close();
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		//options.addOption(OptionBuilder.withArgName("annotationPath").hasArg().withDescription("Path to the GTF file containing the coordinates to extract").create(Thunder.OPT_PATH_ANNOTATION));
		//options.addOption(OptionBuilder.withArgName("results.xprs").hasArg().withDescription("[optional] Text file containing RNA-seq transcript expression quantifications from eXpress").create("e"));
		options.addOption(OptionBuilder.withArgName("alignments").hasArg().withDescription("SAM/BAM file containing alignments").create("f"));
		//options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("Path to which to output results").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("maxRecords").hasArg().withDescription("Maximum number of SAM records to read [default: -1]").create("N"));
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

		if(cmdArgs.hasOption("f")){

			int limitRecords = -1;
			//
			if(cmdArgs.hasOption("N"))
				limitRecords = Integer.valueOf(cmdArgs.getOptionValue("N")).intValue();
			// Parse the alignments
			readSAM(new File(cmdArgs.getOptionValue("f")), limitRecords);


			
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ParseFootprintAlignments", getCmdLineOptions());
			System.out.println();
		}


	}

}
