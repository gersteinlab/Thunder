package fastqTools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import sequenceTools.Aligner_SmithWaterman;
import utils.IO_utils;

public class FindAdapter {

	private boolean _verbose = false;
	private SequenceReader _readReader, _adapterReader;

	public FindAdapter(String readsPath, String actualAdaptersPath){
		try {
			_readReader = new SequenceReader(readsPath);
		} catch (IOException e) {
			System.err.println("ERROR: could not find reads at path "+readsPath);
			e.printStackTrace();
		}

		try{
			_adapterReader = new SequenceReader(actualAdaptersPath);
		} catch (IOException e) {
			System.err.println("ERROR: could not find adapters at path "+actualAdaptersPath);
			e.printStackTrace();
		}

		readAdapters();
		IO_utils.printLineErr("Searching for best of "+_adapterSequences.size()+" potential adapter sequences");
		IO_utils.printLineErr("");
	}

	private void setVerbose(){
		_verbose=true;
	}

	private HashMap<String, String> _adapterSequences = new HashMap<String, String>(); 
	//private HashMap<String, Double> _adapterCurrentScore = new HashMap<String, Double>();
	private HashMap<String, Double> _adapterCumulativeScores = new HashMap<String, Double>();
	private int _nReadsProcessed = 0;
	private double _runningTotalOfScore = 0.0;
	/**
	 * Reads adapters from the given input file to a hashmap, 
	 * - also instantiates an alignment counter for each adapter which is used later to decide on the most likely adapter 
	 */
	private void readAdapters(){
		SequenceRecord thisSequence;
		try {
			while((thisSequence = _adapterReader.readNextRecord()) != null){
				if(_adapterSequences.containsKey(thisSequence)){ 
					System.err.println("WARNING: Ignoring duplicate sequence ID: "+thisSequence.getSequenceID());
				}else{
					_adapterSequences.put(thisSequence.getSequenceID(), thisSequence.getSequence());
					_adapterCumulativeScores.put(thisSequence.getSequenceID(), 0.0);
				}
			}
		} catch (IOException e) {
			System.err.println("ERROR: Failed to read adapter sequences");
			e.printStackTrace();
		}
	}


	//private int _rawReadLength = -1;

	/**
	 * 
	 * @param minimumNumberOfReadsToUse
	 * @param maximumNumberOfReadsToUse
	 * @param ratioOfMostLikelyAdapterToSecondBest
	 * @param minAvgAdapterMatches
	 * @throws IOException
	 */
	public void guessAdapterUsingReads(int minimumNumberOfReadsToUse, int maximumNumberOfReadsToUse, double ratioOfMostLikelyAdapterToSecondBest, double minAvgAdapterMatches) throws IOException{
		String winningAdapter = null;
		Aligner_SmithWaterman sw;
		SequenceRecord thisSequence;

		// For all reads
		while((thisSequence = _readReader.readNextRecord()) != null){
			_nReadsProcessed++;
			Iterator<String> adapters = _adapterSequences.keySet().iterator();
			String thisAdapterID;
			double thisReadMaxScore = 1.0;

			// check read length is consistent
			/*if(_rawReadLength < 0)
				_rawReadLength = thisSequence.getSequence().length();
			else if(_rawReadLength != thisSequence.getSequence().length()){
				IO_utils.printLineErr("ERROR: read length ("+thisSequence.getSequence().length()+") differs from previous length ("+_rawReadLength+") for read "+thisSequence.getSequenceID());
				System.exit(1);
			}*/


			// For all potential adapters
			while(adapters.hasNext()){
				thisAdapterID = adapters.next();

				// do the SW alignment
				sw = new Aligner_SmithWaterman(thisSequence.getSequence(), _adapterSequences.get(thisAdapterID));
				sw.findBestAlignments(false);

				// get alignment score
				//double thisScore = sw.getAlignmentScore()-((sw.getAlignmentStart_reference()-1)*2.0); // subtract ref alignment start to penalise 5' bases of adapter not matching
				/*double thisScore = (sw.getNumberOfMatches()+0.0)*Math.pow(sw.getMatchFractionOfOverlap(), 2);
				if(thisScore < 0  ||  sw.getAlignmentStart_query() < sw.getAlignmentStart_reference())
					thisScore = 0.0;*/
				double thisScore = sw.getWeightedScore();

				// add this score to the global running total
				_runningTotalOfScore += thisScore;

				// add this score to the local max?
				if(thisScore > thisReadMaxScore)
					thisReadMaxScore = thisScore;

				// add this score to this adapter's running total
				_adapterCumulativeScores.put(thisAdapterID, _adapterCumulativeScores.get(thisAdapterID) + thisScore);
				//_adapterCurrentScore.put(thisAdapterID, sw.getAlignmentScore());


				if(_verbose)
					System.err.println(sw.getAlignmentInfo()+"\t"+thisAdapterID+"\tfinalScore: "+thisScore);
				

				//sw.findBestAlignments(true);
				//sw.printAlignmentInfo();
				////sw.printDPMatrix();
				//System.out.println();
			}
			
			if(_verbose)
				System.err.println();
			

			if(_nReadsProcessed >= minimumNumberOfReadsToUse){
				// do we have a consensus on the most likely adapter?
				winningAdapter = getWinningAdapter(ratioOfMostLikelyAdapterToSecondBest);
				if(winningAdapter != null  &&  (_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed) >= minAvgAdapterMatches){
					// we can stop!
					IO_utils.printLineErr("FOUND! Most likely adapter (after "+_nReadsProcessed+" reads): "+winningAdapter+": "+_adapterSequences.get(winningAdapter));
					IO_utils.printLineErr("       with an average alignment score of "+(Math.round(100.0*_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed)/100.0)+" ("+getWinningAdapterScoreOverNextBest()+"x higher than the next best adapter: "+getSecondBestAdapter()+": "+_adapterSequences.get(getSecondBestAdapter())+")");
					break;
				}
			}

			if(_nReadsProcessed >= maximumNumberOfReadsToUse)
				break;
		}

		boolean alignmentAverageScoreOK = true;

		// if we have used all the reads, do we have a consensus on the most likely adapter?
		if(winningAdapter != null){
			if((_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed) < minAvgAdapterMatches){
				// if all reads have been processed but the alignments are all too bad, print the adapters
				IO_utils.printLineErr("WARNING: Unable to determine the most likely adapter (after using ALL "+_nReadsProcessed+" reads and requiring an avarage adapter alignment rate at >= "+minAvgAdapterMatches+"nt)");
				IO_utils.printLineErr("This average adapter alignment rate is too low: "+_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed+"nt");
				IO_utils.printLineErr("");
				//winningAdapter = null;
				alignmentAverageScoreOK = false;
			}
		}


		if(winningAdapter == null){
			winningAdapter = getWinningAdapter(ratioOfMostLikelyAdapterToSecondBest);

			if(winningAdapter != null){
				if((_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed) >= minAvgAdapterMatches){
					// we can stop!
					//IO_utils.printLineErr("FOUND! Most likely adapter (after "+_nReadsProcessed+" reads, which is all that are available) is "+winningAdapter+" with an average alignment score of "+(_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed)+" ("+getWinningAdapterScoreOverNextBest()+"x higher than "+getSecondBestAdapter()+", the next best adapter)");
					IO_utils.printLineErr("FOUND! Most likely adapter (after "+_nReadsProcessed+" reads, which is all that are available): "+winningAdapter+": "+_adapterSequences.get(winningAdapter));
					IO_utils.printLineErr("       with an average alignment score of "+(Math.round(100.0*_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed)/100.0)+" ("+getWinningAdapterScoreOverNextBest()+"x higher than the next best adapter: "+getSecondBestAdapter()+": "+_adapterSequences.get(getSecondBestAdapter())+")");
				}else{
					// if all reads have been processed but the alignments are all too bad, print the adapters
					IO_utils.printLineErr("WARNING: Unable to determine the most likely adapter (after "+_nReadsProcessed+" reads and requiring an avarage adapter alignment rate at >= "+minAvgAdapterMatches+"nt)");
					IO_utils.printLineErr("This average adapter alignment rate: "+_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed+"nt");
					IO_utils.printLineErr("");
					winningAdapter = null;
					alignmentAverageScoreOK = false;
				}
			}

			if(winningAdapter == null){

				if(alignmentAverageScoreOK){
					// if all reads have been processed and there is no consensus on the most likely adapter, print the adapters
					IO_utils.printLineErr("WARNING: Unable to determine the most likely adapter (after "+_nReadsProcessed+" reads and requiring the best adapter to be "+ratioOfMostLikelyAdapterToSecondBest+"x more likely than the next best)");
					IO_utils.printLineErr("");
				}

				IO_utils.printLineErr("The scores of all adapters are given below for debugging:");
				IO_utils.printLineErr("");
				IO_utils.printLineErr("AdapterID\tAvgAlignmentScore\tFractionOfTotalScore\tBestScoreOverThisScore");

				TreeMap<String,Double> adapterCumulativeScores_sorted = new TreeMap<String,Double>(new ValueComparator_Double(_adapterCumulativeScores));
				adapterCumulativeScores_sorted.putAll(_adapterCumulativeScores);
				Iterator<String> it = adapterCumulativeScores_sorted.keySet().iterator();
				int count = 0;
				double bestScore = 0.0;
				while(it.hasNext()){
					String thisAdapter = it.next();

					if(count == 0)
						bestScore = _adapterCumulativeScores.get(thisAdapter);

					IO_utils.printLineErr(thisAdapter+"\t"+(_adapterCumulativeScores.get(thisAdapter)/_nReadsProcessed)+"\t"+(_adapterCumulativeScores.get(thisAdapter) / _runningTotalOfScore)+"\t"+(bestScore/_adapterCumulativeScores.get(thisAdapter)));
					count++;
				}
			}
		}

		// Finally, print the winning adapter sequence to stdout
		if(winningAdapter != null  &&  (_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed) >= minAvgAdapterMatches){
			System.out.println(_adapterSequences.get(winningAdapter));
		}else{
			System.out.print("");
		}

	}

	private void outputStats(String path) throws IOException{
		BufferedWriter bw = new BufferedWriter(new FileWriter(path));

		bw.write("AdapterID\tAvgAlignmentScore\tFractionOfTotalScore\tBestScoreOverThisScore\n");

		TreeMap<String,Double> adapterCumulativeScores_sorted = new TreeMap<String,Double>(new ValueComparator_Double(_adapterCumulativeScores));
		adapterCumulativeScores_sorted.putAll(_adapterCumulativeScores);
		Iterator<String> it = adapterCumulativeScores_sorted.keySet().iterator();
		int count = 0;
		double bestScore = 0.0;
		while(it.hasNext()){
			String thisAdapter = it.next();

			if(count == 0)
				bestScore = _adapterCumulativeScores.get(thisAdapter);

			bw.write(thisAdapter+"\t"+(_adapterCumulativeScores.get(thisAdapter)/_nReadsProcessed)+"\t"+(_adapterCumulativeScores.get(thisAdapter) / _runningTotalOfScore)+"\t"+(bestScore/_adapterCumulativeScores.get(thisAdapter))+"\n");
			count++;
		}

		bw.flush();
		bw.close();
	}


	private Double getWinningAdapterScoreOverNextBest(){
		TreeMap<String,Double> adapterCumulativeScores_sorted = new TreeMap<String,Double>(new ValueComparator_Double(_adapterCumulativeScores));
		adapterCumulativeScores_sorted.putAll(_adapterCumulativeScores);
		Iterator<String> it = adapterCumulativeScores_sorted.keySet().iterator();
		String bestAdapter = it.next();
		String secondBestAdapter = it.next();
		return Math.round(10.0*_adapterCumulativeScores.get(bestAdapter) / _adapterCumulativeScores.get(secondBestAdapter))/10.0;
	}


	/**
	 * If there is a single adapter that has a total cumulative score greater than the
	 * threshold, return it's ID.  Else return null 
	 * @param totalScoreFraction
	 * @return
	 */
	private String getWinningAdapter(double secondBestScoreOverbestScore){
		TreeMap<String,Double> adapterCumulativeScores_sorted = new TreeMap<String,Double>(new ValueComparator_Double(_adapterCumulativeScores));
		adapterCumulativeScores_sorted.putAll(_adapterCumulativeScores);
		Iterator<String> it = adapterCumulativeScores_sorted.keySet().iterator();
		String winningAdapterID = null;

		String bestAdapter = it.next();
		String secondBestAdapter = it.next();
		if(_adapterCumulativeScores.get(bestAdapter) / _adapterCumulativeScores.get(secondBestAdapter) >= secondBestScoreOverbestScore){
			winningAdapterID = bestAdapter;
		}else{
			// if this winning adapter is Illumina_1.5_smallRNA_3p and the next best is Illumina_1.0_smallRNA_3p, pick the former!
			winningAdapterID = "Illumina_1.5_smallRNA_3p";
		}

		return winningAdapterID;
	}

	private String getSecondBestAdapter(){
		TreeMap<String,Double> adapterCumulativeScores_sorted = new TreeMap<String,Double>(new ValueComparator_Double(_adapterCumulativeScores));
		adapterCumulativeScores_sorted.putAll(_adapterCumulativeScores);
		Iterator<String> it = adapterCumulativeScores_sorted.keySet().iterator();

		//String bestAdapter = it.next();
		it.next();
		String secondBestAdapter = it.next();

		return secondBestAdapter;
	}







	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		//options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("output sequences shorter than the maximum length to this file [if not specified, sequences are printed to stdout]").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("Minimum number of reads to use for adapter identification [default: 10000]").create("n"));
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("Maximum number of reads to use for adapter identification [default: 1000000]").create("m"));
		options.addOption(OptionBuilder.withArgName("double").hasArg().withDescription("Minimum score ratio of the most likely adapter to the second best [default: 1.1]").create("r"));
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("Minimum average number of aligned read-adapter bases [default: 7]").create("s"));
		options.addOption(OptionBuilder.withArgName("path").hasArg().withDescription("Path to a fasta format file containing the adapter sequences to use in the search").create("a"));
		options.addOption(OptionBuilder.withArgName("path").hasArg().withDescription("Write a verbose log of the adapter ID process to the given file").create("summary"));
		options.addOption(OptionBuilder.withArgName("").withDescription("print verbose summary of each read to stderr").create("v"));
		return options;
	}

	public static void main(String[] args) throws Exception {

		//args = new String[]{"FindAdapter", "-n","10", "-v", "-r","1.3", "-s","7", "-a","/Users/robk/Downloads/adapters.fa", "/Users/robk/Downloads/tmp_reads_s.fq"};
		//args = new String[]{"FindAdapter", "-v", "-n","10", "-r","1.3", "-s","7", "-a","/Users/robk/Downloads/adapters.fa", "--summary","/Users/robk/Downloads/tmp_reads_s.adapterLog", "/Users/robk/Downloads/tmp_reads_s.fq"};
		//args = new String[]{"FindAdapter", "-n","100", "-r","1.3", "-s","7", "-v", "-a","/Users/robk/Downloads/adapters.fa", "--summary","/Users/robk/Downloads/miR-Pool-5pmoles-TruSeqRecJ-4N-2_S8_L001_R1_001.adapterLog", "/Users/robk/Downloads/miR-Pool-5pmoles-TruSeqRecJ-4N-2_S8_L001_R1_001.fastq"};
		
		//args = new String[]{"FindAdapter", "-n","1000", "-v", "-m","100000", "-s","7", "-a","/Users/robk/Downloads/adapters.fa", "/Users/robk/Downloads/MixA-a_S12_L001_R1_001.fastq", "-summary","/Users/robk/Downloads/ADAPTER_SUMMARY.txt"};
		//args = new String[]{"FindAdapter", "-n","1000", "-v", "-m","100000", "-s","7", "-a","/Users/robk/Downloads/adapters.fa", "/Users/robk/Downloads/ERCCF-Zymo-50-c20-1_S4_L001_R1_001.fastq", "-summary","/Users/robk/Downloads/ADAPTER_SUMMARY.txt"};
		
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		if(cmdArgs.getArgList().size() >= 2  &&  cmdArgs.hasOption("a")){

			// override default minimum number of reads?
			int minNumberOfreads = 10000;
			if(cmdArgs.hasOption("n"))
				minNumberOfreads = Integer.valueOf(cmdArgs.getOptionValue("n")).intValue();

			// override default minimum number of reads?
			int maxNumberOfreads = 1000000;
			if(cmdArgs.hasOption("m"))
				maxNumberOfreads = Integer.valueOf(cmdArgs.getOptionValue("m")).intValue();

			// override default best to second-best adapter score ratio?
			double ratioOfMostLikelyAdapterToSecondBest = 1.1;
			if(cmdArgs.hasOption("r"))
				ratioOfMostLikelyAdapterToSecondBest = Double.valueOf(cmdArgs.getOptionValue("r")).doubleValue();

			// override default minimum alignment rate?
			double minAvgAdapterMatches = 7.0;
			if(cmdArgs.hasOption("s"))
				minAvgAdapterMatches = Double.valueOf(cmdArgs.getOptionValue("s")).doubleValue();

			//
			//Iterator<String> it = cmdArgs.getArgList().iterator();
			//while(it.hasNext()){
			FindAdapter engine = new FindAdapter(cmdArgs.getArgList().get(1).toString(), cmdArgs.getOptionValue("a"));
			
			if(cmdArgs.hasOption("v"))
				engine.setVerbose();
			
			engine.guessAdapterUsingReads(minNumberOfreads, maxNumberOfreads, ratioOfMostLikelyAdapterToSecondBest, minAvgAdapterMatches);
			//}

			if(cmdArgs.hasOption("summary"))
				engine.outputStats(cmdArgs.getOptionValue("summary"));

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" FindAdapter [options] <sequenceFile (or '-' for stdin)>", "", getCmdLineOptions(), "");
			System.out.println();
		}

	}

}



class ValueComparator_Double implements Comparator<String> {

	Map<String, Double> base;
	public ValueComparator_Double(Map<String, Double> base) {
		this.base = base;
	}

	// Note: this comparator imposes orderings that are inconsistent with equals.    
	public int compare(String a, String b) {
		if (base.get(a) >= base.get(b)) {
			return -1;
		} else {
			return 1;
		} // returning 0 would merge keys
	}
}
