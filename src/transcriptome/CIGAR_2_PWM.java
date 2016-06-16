package transcriptome;

import java.io.File;
import java.util.ArrayList;
import java.util.Iterator;
import java.util.List;

import main.Thunder;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import objects.SAMRecordReduced;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import samTools.SAMReader;
import utils.IO_utils;
import utils.Object_utils;

public class CIGAR_2_PWM {


	public static void processReads(File input, int maxReads, boolean verbose){
		SAMReader engine = new SAMReader(input);

		//if(verbose)
		//	System.err.println("BAM sorted by readID: "+engine.isSortedByReadID());
		/*if( ! engine.isSortedByReadID()){
			System.err.println("BAM must be sorted by readID ('queryname')\nThis BAM is sorted by: "+engine.getSortOrder().toString());
			System.exit(0);
		}*/


		int reads = 0;
		int count = 0;

		int maxReadLength_guessed = 10000;

		int[] m = Object_utils.initArray_int(maxReadLength_guessed);;
		int[] s = Object_utils.initArray_int(maxReadLength_guessed);;
		int[] i = Object_utils.initArray_int(maxReadLength_guessed);;

		// re-set max read length to allow us to learn the actual max length from the reads 
		int maxReadLength_actual = 0;

		while(true){
			count = 0;
			ArrayList<SAMRecordReduced> alignmentsForThisRead = engine.getAlignmentsForNextRead();

			// if we've run out of reads
			if(alignmentsForThisRead == null  ||  alignmentsForThisRead.size() == 0)
				break;

			// find best alignment
			SAMRecordReduced bestAlignment = null;
			for(SAMRecordReduced thisRead: alignmentsForThisRead){
				if(thisRead.isPrimaryAlignment())
					bestAlignment = thisRead;
			}

			if(verbose)
				IO_utils.printLineErr("best alignment for read "+bestAlignment.getReadName()+" is "+bestAlignment.getReferenceName()+" (CIGAR:"+bestAlignment.getCigar()+")");

			// If readLength of best alignmnt is larger than the current max read length:
			if(bestAlignment.getReadSequence().length() > maxReadLength_actual){
				maxReadLength_actual = bestAlignment.getReadSequence().length();

				if(maxReadLength_actual > maxReadLength_guessed){
					// TODO, throw a warning!
				}
			}

			// Grab the CIGAR for this read, reverse if necessary
			List<CigarElement> tmpCigarElements = bestAlignment.getCigar().getCigarElements();
			List<CigarElement> cigarElements = new ArrayList<CigarElement>();

			if(bestAlignment.getReadNegativeStrandFlag()){						// if the read is negative strand, reverse the CIGAR elements 
				for(int index=(tmpCigarElements.size()-1);index>=0;index--){
					cigarElements.add(tmpCigarElements.get(index));
				}
			}else{																// if the read is positive strand, leave alone
				cigarElements = tmpCigarElements;
			}

			if(verbose)
				IO_utils.printLineErr("CIGAR size = "+cigarElements.size()+"");

			/*
			 * Loop through CIGAR string and add to the overall stats
			 */
			Iterator<CigarElement> cigarIterator = cigarElements.iterator();
			CigarElement tmp;
			CigarOperator tmp_op;
			while(cigarIterator.hasNext()){
				tmp = cigarIterator.next();
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
						System.err.println("Unknown CIGAR operator: "+tmp_op.toString()+" in "+bestAlignment.getCigar().toString());
					}
				}
			}

			reads ++;
			if(maxReads > 0  &  reads >= maxReads)
				break;
		}

		/*
		 * Print the result to stdout
		 */
		IO_utils.printLineErr("SAM records read = "+reads);
		IO_utils.printLineErr("max observed read length = "+maxReadLength_actual);

		if(reads > 0  &&  maxReadLength_actual > 0){

			System.out.print("M\t");
			//for(int x=0;x<m.length-1;x++)
			for(int x=0;x<(maxReadLength_actual-1);x++)
				System.out.print((m[x]+0.0)/reads+"\t");
			//System.out.println((m[m.length-1]+0.0)/reads);
			System.out.println((m[maxReadLength_actual-1]+0.0)/reads);

			System.out.print("S\t");
			//for(int x=0;x<s.length-1;x++)
			for(int x=0;x<(maxReadLength_actual-1);x++)
				System.out.print((s[x]+0.0)/reads+"\t");
			//System.out.println((s[s.length-1]+0.0)/reads);
			System.out.println((s[maxReadLength_actual-1]+0.0)/reads);

			System.out.print("I\t");
			//for(int x=0;x<i.length-1;x++)
			for(int x=0;x<(maxReadLength_actual-1);x++)
				System.out.print((i[x]+0.0)/reads+"\t");
			//System.out.println((i[i.length-1]+0.0)/reads);
			System.out.println((i[maxReadLength_actual-1]+0.0)/reads);

		}
		engine.close();
	}






	/*
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
	 */


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
		options.addOption(OptionBuilder.withDescription("Print *really* verbose output. Debugging only. Not recommended with large input BAM!").create("v"));
		return options;
	}


	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {

		/*args = new String[]{"CIGAR_2_PWM", 
				"-f", "/Users/robk/Box Sync/Work/HEK293_RNAseq/STAR_eXpress/D1-Footprint_1x50_gencode.v18.annotation_mappedReads.bam",
				"-N","10",
				"-v"
				//"-a", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v18.annotation.filtered.gtf",
				//"-e", "/Users/robk/Box Sync/Work/HEK293_RNAseq/STAR_eXpress/A1-Total_2x75/results.xprs",
				//"-o", "/Users/robk/Box Sync/Work/HEK293_RNAseq/EM_isoformQuant/firstRound_fullFootprintEM_fromJing/FootPrint/data/TEST_NEW"};
		};*/

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption("f")){

			int limitRecords = -1;
			//
			if(cmdArgs.hasOption("N"))
				limitRecords = Integer.valueOf(cmdArgs.getOptionValue("N")).intValue();

			boolean verbose = false;
			if(cmdArgs.hasOption("v"))
				verbose = true;

			// Parse the alignments
			//readSAM(new File(cmdArgs.getOptionValue("f")), limitRecords);
			processReads(new File(cmdArgs.getOptionValue("f")), limitRecords, verbose);



		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" CIGAR_2_PWM", getCmdLineOptions());
			System.out.println();
		}


	}

}
