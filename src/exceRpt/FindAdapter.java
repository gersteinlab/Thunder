package exceRpt;

import java.io.IOException;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import main.Thunder;
import sequenceTools.Aligner_SmithWaterman;
import fastqTools.SequenceReader;
import fastqTools.SequenceRecord;

public class FindAdapter {

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
		Thunder.printLineErr("Searching for best of "+_adapterSequences.size()+" potential adapter sequences");
		Thunder.printLineErr("");
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



	public void guessAdapterUsingReads(int minimumNumberOfReadsToUse, int maximumNumberOfReadsToUse, double ratioOfMostLikelyAdapterToSecondBest, double minAvgAdapterMatches) throws IOException{
		String winningAdapter = null;
		Aligner_SmithWaterman sw;
		SequenceRecord thisSequence;
		while((thisSequence = _readReader.readNextRecord()) != null){
			_nReadsProcessed++;
			Iterator<String> adapters = _adapterSequences.keySet().iterator();
			String thisAdapterID;
			double thisReadMaxScore = 1.0;
			while(adapters.hasNext()){
				thisAdapterID = adapters.next();

				// do the SW alignment
				sw = new Aligner_SmithWaterman(thisSequence.getSequence(), _adapterSequences.get(thisAdapterID));
				sw.findBestAlignments(false);
				
				// print alignment info for debug?
				//System.out.println(thisAdapterID);
				//System.out.println("sw.getAlignmentStart_reference(): "+sw.getAlignmentStart_reference());
				//System.out.println("sw.getAlignmentStart_query(): "+sw.getAlignmentStart_query());
				//sw.printAlignmentInfo();
				//System.out.println();
				
				
				// get alignment score
				double thisScore = sw.getAlignmentScore()-((sw.getAlignmentStart_reference()-1)*2.0); // subtract ref alignment start to penalise 5' bases of adapter not matching
				if(thisScore < 0)
					thisScore = 0.0;
				
				
				// add this score to the global running total
				_runningTotalOfScore += thisScore;

				// add this score to the local max?
				if(thisScore > thisReadMaxScore)
					thisReadMaxScore = thisScore;  
						
				// add this score to this adapter's running total
				_adapterCumulativeScores.put(thisAdapterID, _adapterCumulativeScores.get(thisAdapterID) + thisScore);
				//_adapterCurrentScore.put(thisAdapterID, sw.getAlignmentScore());


				//sw.findBestAlignments(true);
				//sw.printAlignmentInfo();
				////sw.printDPMatrix();
				//System.out.println();
			}

			if(_nReadsProcessed >= minimumNumberOfReadsToUse){
				// do we have a consensus on the most likely adapter?
				winningAdapter = getWinningAdapter(ratioOfMostLikelyAdapterToSecondBest);
				if(winningAdapter != null){
					// we can stop!
					Thunder.printLineErr("FOUND! Most likely adapter (after "+_nReadsProcessed+" reads) is "+winningAdapter+" with an average alignment score of "+(_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed)+" ("+getWinningAdapterScoreOverNextBest()+"x higher than "+getSecondBestAdapter()+", the next best adapter)");
					break;
				}
			}
			
			if(_nReadsProcessed >= maximumNumberOfReadsToUse)
				break;
		}

		boolean alignmentAverageScoreOK = true;

		// if we have used all the reads, do we have a consensus on the most likely adapter?
		if(winningAdapter == null){
			winningAdapter = getWinningAdapter(ratioOfMostLikelyAdapterToSecondBest);

			if(winningAdapter != null){
				if(_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed >= minAvgAdapterMatches){
					// we can stop!
					Thunder.printLineErr("FOUND! Most likely adapter (after "+_nReadsProcessed+" reads, which is all that are available) is "+winningAdapter+" with an average alignment score of "+(_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed)+" ("+getWinningAdapterScoreOverNextBest()+"x higher than "+getSecondBestAdapter()+", the next best adapter)");
				}else{
					// if all reads have been processed but the alignments are all too bad, print the adapters
					Thunder.printLineErr("WARNING: Unable to determine the most likely adapter (after "+_nReadsProcessed+" reads and requiring an avarage adapter alignment rate at >= "+minAvgAdapterMatches+"nt)");
					Thunder.printLineErr("This average adapter alignment rate: "+_adapterCumulativeScores.get(winningAdapter)/_nReadsProcessed+"nt");
					Thunder.printLineErr("");
					winningAdapter = null;
					alignmentAverageScoreOK = false;
				}
			}

			if(winningAdapter == null){

				if(alignmentAverageScoreOK){
					// if all reads have been processed and there is no consensus on the most likely adapter, print the adapters
					Thunder.printLineErr("WARNING: Unable to determine the most likely adapter (after "+_nReadsProcessed+" reads and requiring the best adapter to be "+ratioOfMostLikelyAdapterToSecondBest+"x more likely than the next best)");
					Thunder.printLineErr("");
				}

				Thunder.printLineErr("The scores of all adapters are given below for debugging:");
				Thunder.printLineErr("");
				Thunder.printLineErr("AdapterID\tAvgAlignmentScore\tFractionOfTotalScore\tBestScoreOverThisScore");

				TreeMap<String,Double> adapterCumulativeScores_sorted = new TreeMap<String,Double>(new ValueComparator_Double(_adapterCumulativeScores));
				adapterCumulativeScores_sorted.putAll(_adapterCumulativeScores);
				Iterator<String> it = adapterCumulativeScores_sorted.keySet().iterator();
				int count = 0;
				double bestScore = 0.0;
				while(it.hasNext()){
					String thisAdapter = it.next();

					if(count == 0)
						bestScore = _adapterCumulativeScores.get(thisAdapter);

					Thunder.printLineErr(thisAdapter+"\t"+(_adapterCumulativeScores.get(thisAdapter)/_nReadsProcessed)+"\t"+(_adapterCumulativeScores.get(thisAdapter) / _runningTotalOfScore)+"\t"+(bestScore/_adapterCumulativeScores.get(thisAdapter)));
					count++;
				}
			}
		}

		// Finally, print the winning adapter sequence to stdout
		if(winningAdapter != null){
			System.out.println(_adapterSequences.get(winningAdapter));
		}else{
			System.out.print("");
		}

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
		}

		return winningAdapterID;
	}
	
	private String getSecondBestAdapter(){
		TreeMap<String,Double> adapterCumulativeScores_sorted = new TreeMap<String,Double>(new ValueComparator_Double(_adapterCumulativeScores));
		adapterCumulativeScores_sorted.putAll(_adapterCumulativeScores);
		Iterator<String> it = adapterCumulativeScores_sorted.keySet().iterator();
		
		String bestAdapter = it.next();
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
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("Minimum number of reads to use for adapter identification [default: 1000]").create("n"));
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("Maximum number of reads to use for adapter identification [default: 100000]").create("m"));
		options.addOption(OptionBuilder.withArgName("double").hasArg().withDescription("Minimum score ratio of the most likely adapter to the second best [default: 1.1]").create("r"));
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("Minimum average number of aligned read-adapter bases [default: 7]").create("s"));
		options.addOption(OptionBuilder.withArgName("path").hasArg().withDescription("Path to a fasta format file containing the adapter sequences to use in the search").create("a"));
		return options;
	}

	public static void main(String[] args) throws Exception {

		//args = new String[]{"FindAdapter", "-n","1000", "-r","1.3", "-s","7", "-a","/Users/robk/Downloads/adapters.fa", "/Users/robk/Downloads/tmp_reads_s.fq"};
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		if(cmdArgs.getArgList().size() >= 2  &&  cmdArgs.hasOption("a")){

			// override default minimum number of reads?
			int minNumberOfreads = 1000;
			if(cmdArgs.hasOption("n"))
				minNumberOfreads = Integer.valueOf(cmdArgs.getOptionValue("n")).intValue();

			// override default minimum number of reads?
			int maxNumberOfreads = 100000;
			if(cmdArgs.hasOption("m"))
				maxNumberOfreads = Integer.valueOf(cmdArgs.getOptionValue("m")).intValue();

			// override default best to second-best adapter score ratio?
			double ratioOfMostLikelyAdapterToSecondBest = 1.1;
			if(cmdArgs.hasOption("r"))
				ratioOfMostLikelyAdapterToSecondBest = Double.valueOf(cmdArgs.getOptionValue("r")).doubleValue();

			// override default minimum alignment rate?
			double minAvgAdapterMatches = 7;
			if(cmdArgs.hasOption("s"))
				minAvgAdapterMatches = Double.valueOf(cmdArgs.getOptionValue("s")).doubleValue();

			//
			//Iterator<String> it = cmdArgs.getArgList().iterator();
			//while(it.hasNext()){
			FindAdapter engine = new FindAdapter(cmdArgs.getArgList().get(1).toString(), cmdArgs.getOptionValue("a"));
			engine.guessAdapterUsingReads(minNumberOfreads, maxNumberOfreads, ratioOfMostLikelyAdapterToSecondBest, minAvgAdapterMatches);
			//}

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
