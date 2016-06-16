package proteome;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import main.Thunder;
import objects.FastX;
import objects.FastX_Record;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

import utils.IO_utils;

public class ProcessPepDigest {

	public ProcessPepDigest(){

	}

	//private HashMap<String, HashMap<Integer,Integer>> isoformPeptideCounts = new HashMap<String, HashMap<Integer,Integer>>();
	//public HashMap<String, HashMap<Integer,Integer>> getIsoformPeptideCounts(){ return isoformPeptideCounts; } 

	public static HashMap<String, HashMap<Integer,Integer>> processEMBOSSoutput(File input) throws IOException{

		HashMap<String, HashMap<Integer,Integer>> isoformPeptideCounts = new HashMap<String, HashMap<Integer,Integer>>();

		BufferedReader in = new BufferedReader(new FileReader(input));
		String line = "";
		String[] bits;
		int peptideLength = 0;

		String thisIsoform = "";
		int thisFrame = 0;



		while((line=in.readLine()) != null){
			if(line.startsWith("# Sequence:")){
				bits = line.split(" ")[2].split("_");
				thisIsoform = bits[0];
				thisFrame = Integer.valueOf(bits[1]).intValue();
				//System.out.println(thisIsoform+"\t"+thisFrame);

				if(!isoformPeptideCounts.containsKey(thisIsoform))
					isoformPeptideCounts.put(thisIsoform, new HashMap<Integer,Integer>());

				isoformPeptideCounts.get(thisIsoform).put(thisFrame, 0);


				// Keep reading lines until we see some useful data!
				while((line=in.readLine()) != null){
					if(line.trim().startsWith("Start")){
						line=in.readLine();					// burn the header line
						break;
					}
				}
			}



			if(!line.startsWith("#")){
				bits = line.split(" ");
				if(bits.length > 10){
					//System.out.println(line);
					//System.out.println(bits.length);
					//System.out.println(bits[bits.length-1]);

					peptideLength = bits[bits.length-1].length();
					//System.out.println(peptideLength);
					if(peptideLength >= 6  &&  peptideLength <= 30){
						isoformPeptideCounts.get(thisIsoform).put(thisFrame, isoformPeptideCounts.get(thisIsoform).get(thisFrame) + 1);
					}
				}
			}
		}
		in.close();

		return(isoformPeptideCounts);
	}



	public static HashMap<String, Integer> digestProteinsFromFasta_withStats(File inputFasta, ArrayList<String> enzymes) throws IOException{
		HashMap<String, Integer> recordID_2_peptideCount = new HashMap<String, Integer>();
		HashMap<String, Integer> peptideSeq_2_peptideCount = new HashMap<String, Integer>();
		HashMap<String, Integer> peptideSeq_2_peptideLength = new HashMap<String, Integer>();

		String[] bits;
		FastX fastx = new FastX(inputFasta);
		FastX_Record thisRecord;
		while((thisRecord = fastx.readNext()) != null){

			// Separate out based on STOP codons
			String result = "";
			bits = thisRecord.getSequence().split("\\*");
			for(int i=0;i<bits.length;i++)
				result = result.concat(bits[i]+" ");

			if(enzymes.contains(ENZYME_TRYPSIN))
				result = digest_Trypsin(result);
			if(enzymes.contains(ENZYME_LYSC))
				result = digest_LysC(result);

			//System.out.println(result);

			recordID_2_peptideCount.put(thisRecord.getID(), countPeptides(result, 6, 30));

			String[] bits2 = result.split(" ");
			for(int i=0;i<bits2.length;i++){
				if(isEligiblePeptide(bits2[i], 6, 30)){
					if(!peptideSeq_2_peptideCount.containsKey(bits2[i])){
						peptideSeq_2_peptideCount.put(bits2[i], 0);
						peptideSeq_2_peptideLength.put(bits2[i], bits2[i].length());
					}
					peptideSeq_2_peptideCount.put(bits2[i], peptideSeq_2_peptideCount.get(bits2[i])+1);
				}
			}
		}

		int totalAAmass = 0;
		Iterator<String> it = peptideSeq_2_peptideLength.keySet().iterator();
		while(it.hasNext())
			totalAAmass += peptideSeq_2_peptideLength.get(it.next());

		System.out.println("# unique peptide sequences: "+peptideSeq_2_peptideCount.size());
		System.out.println("total peptide mass: "+totalAAmass);


		return(recordID_2_peptideCount);
	}

	public static HashMap<String, Integer> digestProteinsFromFasta(File inputFasta, ArrayList<String> enzymes) throws IOException{
		HashMap<String, Integer> results = new HashMap<String, Integer>();
		String[] bits;
		FastX fastx = new FastX(inputFasta);
		FastX_Record thisRecord;
		while((thisRecord = fastx.readNext()) != null){

			// Separate out based on STOP codons
			String result = "";
			bits = thisRecord.getSequence().split("\\*");
			for(int i=0;i<bits.length;i++)
				result = result.concat(bits[i]+" ");

			if(enzymes.contains(ENZYME_TRYPSIN))
				result = digest_Trypsin(result);
			if(enzymes.contains(ENZYME_LYSC))
				result = digest_LysC(result);

			//System.out.println(result);

			results.put(thisRecord.getID(), countPeptides(result, 6, 30));

		}

		return(results);
	}



	/**
	 * Count the observable peptide products of this digest
	 * @param seq
	 * @param minSize
	 * @param maxSize
	 * @return
	 */
	public static int countPeptides(String seq, int minSize, int maxSize){
		int peptideCount = 0;
		String[] bits = seq.split(" ");
		for(int i=0;i<bits.length;i++)
			//if(bits[i].length() >= minSize  &&  bits[i].length() <= maxSize)
			if(isEligiblePeptide(bits[i], minSize, maxSize))
				peptideCount ++;
		return(peptideCount);
	}


	/**
	 * Is this peptide is within an acceptible size range (based on what we can observe on the mass-spec)
	 * @param seq
	 * @param minSize
	 * @param maxSize
	 * @return
	 */
	public static boolean isEligiblePeptide(String seq, int minSize, int maxSize){
		if(seq.length() >= minSize  &&  seq.length() <= maxSize)
			return true;
		else
			return false;
	}


	public static String digest_Trypsin(String inputSeq){
		String result = "";
		String tmpK = "";
		String tmpR = "";
		String[] bits;

		inputSeq = inputSeq.replaceAll("KP", "-");
		inputSeq = inputSeq.replaceAll("RP", "=");

		bits = inputSeq.split("K");
		for(int i=0;i<bits.length;i++)
			tmpK = tmpK.concat(bits[i]+"K ");
		//System.out.println(tmpK);

		bits = tmpK.split("R");
		for(int i=0;i<bits.length;i++)
			tmpR = tmpR.concat(bits[i]+"R ");
		//System.out.println(tmpR);

		result = tmpR.replaceAll("-", "KP");
		result = result.replaceAll("=", "KR");
		return(result);	
	}

	public static String digest_LysC(String inputSeq){
		String result = "";
		String[] bits = inputSeq.split("K");
		for(int i=0;i<bits.length;i++)
			result = result.concat(bits[i]+"K ");
		//System.out.println(result);
		return(result);

	}

	public static final String ENZYME_TRYPSIN = "Trypsin";
	public static final String ENZYME_LYSC = "LysC";

	
	
	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName(".fa").hasArg().withDescription("Path fasta file containing amino acid sequences to digest").create("fasta"));
		options.addOption(OptionBuilder.withArgName("integer").hasArg().withDescription("[optional] minimum peptide length").create("min"));
		options.addOption(OptionBuilder.withArgName("integer").hasArg().withDescription("[optional] maximum peptide length").create("max"));
		return options;
	}
	
	public static void main(String[] args) throws IOException, ParseException{
		String fastaPath = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation.protein.fa";
		//String fastaPath = "/Users/robk/Downloads/test.protein.fa";
		//args = new String[]{"--fasta",fastaPath, "-min","6", "-max","30"};

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		if(cmdArgs.hasOption("fasta")){
			
			IO_utils.printLineErr("reading file: "+cmdArgs.getOptionValue("fasta"));
			
			ArrayList<String> enzymes = new ArrayList<String>();
			enzymes.add(ENZYME_TRYPSIN);
			enzymes.add(ENZYME_LYSC);
			//digestProteinsFromFasta(new File(fastaPath), enzymes);
			digestProteinsFromFasta_withStats(new File(cmdArgs.getOptionValue("fasta")), enzymes);
			
			
			IO_utils.printLineErr("Done");

			//HashMap<String, HashMap<Integer,Integer>> isoformPeptideCounts = ProcessPepDigest.processEMBOSSoutput(new File("/Users/robk/Desktop/EM_TEST/ProteoData/gencode.v18.annotation.protein.digest.head"));
			//HashMap<String, HashMap<Integer,Integer>> isoformPeptideCounts = ProcessPepDigest.processEMBOSSoutput(new File("/Users/robk/Desktop/EM_TEST/ProteoData/gencode.v18.annotation.protein.digest"));
			//System.err.println("Done!");

			/*Iterator<String> it = isoformPeptideCounts.keySet().iterator();
		while(it.hasNext()){
			String transcriptID = it.next();
			System.out.print(transcriptID+":\t");
			System.out.print(isoformPeptideCounts.get(transcriptID).get(1)+"\t");
			System.out.print(isoformPeptideCounts.get(transcriptID).get(2)+"\t");
			System.out.print(isoformPeptideCounts.get(transcriptID).get(3)+"\n");
		}*/
		}
		else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" PeptideDigest", getCmdLineOptions());
			System.err.println();
		}
	}

}
