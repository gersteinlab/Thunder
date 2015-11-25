package footprintAlignments;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;

import main.Thunder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import objects.Transcript;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import utils.IO_utils;
import annotation.ReadGTF;
import annotation.TranscriptAnnotation;


/*protected static class TranscriptSummary {
	String transcript;
	String type;
	int numberOfReads;
	int maxInFrame;  
	int midInFrame;  
	int minInFrame;  

	double alignmentPValue;
	boolean significance;

	TranscriptSummary(String transcriptID, String transcriptType, int totalReads, int maxAligned, int midAligned, int minAligned, double adjustedPValue) {
		transcript = transcriptID;
		type = transcriptType;
		numberOfReads = totalReads;
		maxInFrame = maxAligned;
		midInFrame = midAligned;
		minInFrame = minAligned;
		alignmentPValue = this.getAlignmentPValue();
		significance = alignmentPValue <= adjustedPValue;
	}

	double getAlignmentPValue() {
		// find the absolute number of aligned reads by multiplying the total numberOfReads with the percentage that are in frame
		BigDecimal sum = new BigDecimal("0");  // in the end, we subtract `sum` from 1 to get a probability

		for (int i = 0; i < maxInFrame; i++) {
			for (int j = 0; j < maxInFrame; j++) {
				for (int m = 0; m < maxInFrame; m++) {
					if (i + j + m == numberOfReads) { 
						//                            denominatorDecimal = factorialN.divide(denominatorDecimal, 200, RoundingMode.HALF_UP);
						sum = sum.add(binPatterns(numberOfReads, i, j, m));
					} 
				}
			}
		}

		BigDecimal totalPossibilities = BigDecimal.valueOf(3).pow(numberOfReads);
		BigDecimal matches = totalPossibilities.subtract(sum);

		return matches.divide(totalPossibilities, 200, RoundingMode.HALF_UP).doubleValue(); 

	}

	private BigDecimal binPatterns(int n, int j, int m, int t) {
		return new BigDecimal(getFactorial(n)).divide(new BigDecimal(getFactorial(j).multiply(getFactorial(m)).multiply(getFactorial(t))), 200, RoundingMode.HALF_UP);
	}

	private BigInteger getFactorial(int n) {
		if (factorialCache.containsKey(n)) {
			return factorialCache.get(n);
		} else {
			BigInteger factorialN = factorial(n);
			factorialCache.put(n, factorialN);
			return factorialN;
		}
	}


	public void writeSummaryToFile(BufferedWriter outputFile) throws IOException {
		outputFile.write(transcript + "\t");
		outputFile.write(type + "\t");
		outputFile.write(numberOfReads + "\t");
		outputFile.write(String.valueOf(maxInFrame + "\t"));
		outputFile.write(String.valueOf(midInFrame + "\t"));
		outputFile.write(String.valueOf(minInFrame + "\t"));
		outputFile.write(String.valueOf(alignmentPValue + "\t"));
		outputFile.write(String.valueOf(significance));
		outputFile.newLine();
	}


}*/

//
//class IntegerComparator implements Comparator {
//    private Map base;
//    public IntegerComparator(Map base) { this.base = base; }
//
//    // Note: this comparator imposes orderings that are inconsistent with equals.
//    public int compare(Object a, Object b) {
//        if ((Integer) a >= (Integer) b) {
//            return -1;
//        } else {
//            return 1;
//        } // returning 0 would merge keys
//    }
//}

public class FootprintFrameCalculator {
	private TranscriptAnnotation _annotation;
	private int _minFootprints = 1;

	int numberOfTests;  // the number of transcripts (hypotheses). This will be used to determine our new p-value adjusted multiple testing correction. 
	double adjustedPValue;
	int minimumNumberOfReads;
	//private static final HashMap<Integer, BigInteger> factorialCache = new HashMap<Integer, BigInteger>();  // don't redo computations


	// TODO: perhaps have some other interfaces if clients already have the gencode data
	public FootprintFrameCalculator(String annotation, int minFootprints) throws Exception {
		_annotation = ReadGTF.readGTF(annotation);
		_minFootprints = minFootprints;
		makeReadOffset_guessed();
	}







	// stores transcriptome alignmnets
	private Map<String, ArrayList<SAMRecord>> _transcript2readAlignments = new HashMap<String, ArrayList<SAMRecord>>();

	/**
	 * Reads transcriptome alignments in the given file
	 * @param alignmentFile
	 */
	private void readAlignments(File alignmentFile, int maxReads){
		int count_allReads = 0;
		int count_uniqueReads = 0;
		int count_allAlignments = 0;
		int count_uniqueAlignments = 0;

		SAMFileReader inputBam = new SAMFileReader(alignmentFile);

		// to check for multi-maps:
		String thisReadID="";
		String lastReadID="";
		String thisGeneID="";
		String lastGeneID="";
		boolean allOK = true;
		ArrayList<SAMRecord> tmpReads = new ArrayList<SAMRecord>();

		IO_utils.printLineErr("Reading alignments...");
		for (final SAMRecord thisRead : inputBam) {

			String transcriptID = thisRead.getReferenceName();
			thisReadID = thisRead.getReadName();

			if(!transcriptID.equals("*")){
				count_allAlignments ++;

				thisGeneID = _annotation.getGeneForTranscript(transcriptID);

				if(thisGeneID == null){
					IO_utils.printLineErr("WARNING: no geneID in annotation corresponding to transcript \'"+transcriptID+"\' for read \'"+thisReadID+"\'");
				}else{

					// for the very first read:
					if(count_allAlignments == 1){
						count_allReads ++;
						lastReadID = thisReadID;
						lastGeneID = thisGeneID;
					}

					// check if this is the same read:
					if(thisReadID.equals(lastReadID)){
						// if it is the same read and maps to the same gene, keep it
						if(thisGeneID.equals(lastGeneID)){
							// everything ok so far
							tmpReads.add(thisRead);
						}else{
							// read maps to +1 gene, suppress it
							allOK = false;
						}
					}else{
						// this is a new read
						count_allReads ++;

						if(allOK){
							count_uniqueReads ++;
							// add this read and transcripts to the global list...
							count_uniqueAlignments += tmpReads.size();
							Iterator<SAMRecord> it = tmpReads.iterator();
							while(it.hasNext()){
								SAMRecord tmpRead = it.next();
								String tmpTranscriptID = tmpRead.getReferenceName();
								// add read to transcript 
								if(!_transcript2readAlignments.containsKey(tmpTranscriptID)) {
									_transcript2readAlignments.put(tmpTranscriptID, new ArrayList<SAMRecord>());
								}
								_transcript2readAlignments.get(tmpTranscriptID).add(tmpRead);
							}
						}
						tmpReads = new ArrayList<SAMRecord>();
						tmpReads.add(thisRead);
						allOK = true;
						lastReadID = thisReadID;
						lastGeneID = thisGeneID;
					}
				}
			}
			if(count_allReads >= maxReads){
				IO_utils.printLineErr("Using the first "+maxReads+" reads.");
				break;
			}

		}
		if (inputBam != null) inputBam.close();

		IO_utils.printLineErr("Done - Read "+count_allAlignments+" alignments of "+count_allReads+" reads.");
		IO_utils.printLineErr("     - Discarded "+(count_allAlignments-count_uniqueAlignments)+" alignments ("+((Math.round((count_allAlignments-count_uniqueAlignments)*10000.0/count_allAlignments)/100.0))+"%) from "+(count_allReads-count_uniqueReads)+" reads ("+((Math.round((count_allReads-count_uniqueReads)*10000.0/count_allReads)/100.0))+"%) due to multi-gene ambiguity.");
	}




















	private HashMap<Integer, HashMap<Double, Integer>> _readLengths2offsets = new HashMap<Integer, HashMap<Double, Integer>>();

	/**
	 * 
	 */
	public void printReadLengths2offsets(HashMap<Integer,Double[]> globalOffsets){
		// get the min and max read length
		int minLength=0, maxLength=0;
		Iterator<Integer> lengths = globalOffsets.keySet().iterator();
		while(lengths.hasNext()){
			int thisLength = lengths.next();
			if(thisLength < minLength)
				minLength = thisLength;
			else if(thisLength > maxLength)
				maxLength = thisLength;
		}

		// print the result:
		for(int i=minLength;i<=maxLength;i++){
			if(globalOffsets.containsKey(i)){
				IO_utils.printLineErr("       Read-length:"+i+"  CDS_offset:"+globalOffsets.get(i)[0]+"  weight:"+(Math.round(100.0*globalOffsets.get(i)[1])/100.0));
			}
		}
	}

	/*

	private void printReadLengths2offsets(HashMap<Integer,Double[]> globalOffsets){
		Iterator<Integer> lengths = globalOffsets.keySet().iterator();
		while(lengths.hasNext()){
			int thisLength = lengths.next();
			IO_utils.printLineErr("       Read-length:"+thisLength+"  CDS offset:"+globalOffsets.get(thisLength)[0]+"  weight:"+(Math.round(100.0*globalOffsets.get(thisLength)[1])/100.0));			
		}
	}*/

	/**
	 * 
	 * @param outputPath
	 * @throws IOException 
	 */
	public void writeReadLengths2offsets(String outputPath) throws IOException{

		// get the min and max read length
		int minLength=0, maxLength=0;
		Iterator<Integer> lengths = _readLengths2offsets.keySet().iterator();
		while(lengths.hasNext()){
			int thisLength = lengths.next();
			if(thisLength < minLength)
				minLength = thisLength;
			else if(thisLength > maxLength)
				maxLength = thisLength;
		}

		BufferedWriter out = new BufferedWriter(new FileWriter(outputPath));
		out.write("ReadLength\tOffsetToCDS\tNumberOfReads\n");

		for(int i=minLength;i<=maxLength;i++){
			//int thisLength = i;
			int totalReads = 0;
			int bestReads = 0;
			double bestOffset = 0.0;

			if(_readLengths2offsets.containsKey(i)){
				Iterator<Double> offsets = _readLengths2offsets.get(i).keySet().iterator();
				while(offsets.hasNext()){
					Double thisOffset = offsets.next();
					int nReads = _readLengths2offsets.get(i).get(thisOffset);
					out.write(i+"\t"+thisOffset+"\t"+nReads+"\n");

					totalReads += nReads;
					if(nReads > bestReads){
						bestReads = nReads;
						bestOffset = thisOffset;
					}
				}

				double fractionReads = (bestReads+0.0)/(totalReads+0.0);
				//if(thisLength >= 26  &&  thisLength <= 33  &&  fractionReads > 0.5){
				//if(fractionReads > 0.5){
				if(i >= 20  &&  i <= 40  &&  fractionReads > 0.5){
					_readLengthWeights_empirical.put(i, new Double[]{bestOffset, fractionReads});
				}
			}
		}
		out.flush();
		out.close();
	}

	// stores the learned offsets for different read lengths
	// <readLength, [offset, fractionOfReadsOfThisLengthWithThisOffset]>
	private HashMap<Integer,Double[]> _readLengthWeights_empirical = new HashMap<Integer,Double[]>();
	private HashMap<Integer,Double[]> _readLengthWeights_guessed = new HashMap<Integer,Double[]>();


	/**
	 * Create read length vs CDS offsets based on observation that 3' variability more likely than 5'
	 * @param readLength
	 * @return
	 */
	private void makeReadOffset_guessed(){
		//double weight = 0.01;
		double weight = 1.0;
		_readLengthWeights_guessed.put(26, new Double[]{2.0,weight});
		_readLengthWeights_guessed.put(27, new Double[]{2.5,weight});
		_readLengthWeights_guessed.put(28, new Double[]{0.0,weight});
		_readLengthWeights_guessed.put(29, new Double[]{0.5,weight});
		_readLengthWeights_guessed.put(30, new Double[]{1.0,weight});
		_readLengthWeights_guessed.put(31, new Double[]{0.5,weight});
		_readLengthWeights_guessed.put(32, new Double[]{1.0,weight});
	}



	/**
	 * Calculate the distance in nucleotides of the mid-point of these reads to the annotated CDS translation frame
	 * @param thisTranscript
	 * @return
	 */
	private HashMap<SAMRecord, Double> calculateOffsetForReadsInThisTranscript(Transcript thisTranscript){
		ArrayList<SAMRecord> reads = _transcript2readAlignments.get(thisTranscript.getID());
		HashMap<SAMRecord, Double> readOffsets = new HashMap<SAMRecord, Double>();

		// compute the CDS offset from the beginning/end of this transcript (to be able to pick CDS reads and to compute distance from start codon)
		int cdsOffset_start = thisTranscript.getCDSStartFromTxStartInTxCoords();
		int cdsOffset_end = thisTranscript.getCDSStopFromTxStopInTxCoords();

		/*if(thisTranscript.getID().equals("ENST00000331710.7")){
			System.out.print(thisTranscript.getID());
			System.out.println("\t"+cdsOffset_start);
		}*/

		// loop through each read
		Iterator<SAMRecord> readsIterator = reads.iterator();
		while(readsIterator.hasNext()){
			SAMRecord thisRead = readsIterator.next();

			// read strand must be +ve w.r.t the transcript!
			boolean negativeStrand = thisRead.getReadNegativeStrandFlag();

			// if the middle of this read is within the CDS
			double readMid = thisRead.getAlignmentStart()+0.0 + ((thisRead.getCigar().getReferenceLength()+0.0)/2.0);
			if(readMid >= cdsOffset_start  &&  readMid >= cdsOffset_end  &&  !negativeStrand){
				// record the offset of this read from the start of the CDS
				readOffsets.put(thisRead, (readMid - cdsOffset_start) % 3);
			}else{
				readOffsets.put(thisRead, -999.0);
			}
		}
		return readOffsets;
	}

	//private static final int _USE_OFFSET_NONE = 0;
	private static final int _USE_OFFSET_GUESSED = 1;
	private static final int _USE_OFFSET_EMPIRICAL = 2;


	/**
	 * @throws IOException 
	 * 
	 */
	public void computeReadMidOffset2CDS(String outputPath, int useOffsets, boolean useFrameToSelectBestTranscript) throws IOException{

		HashMap<Integer, Double[]> globalOffsets = _readLengthWeights_guessed;
		if(useOffsets == _USE_OFFSET_GUESSED){
			//IO_utils.printLineErr("Learning global read-length vs CDS offset parameters using guessed read-length frame offsets:");
			IO_utils.printLineErr("Learning global read-length vs CDS offset parameters...");
			//printReadLengths2offsets(globalOffsets);
			//IO_utils.printLineErr("");
		}else if(useOffsets == _USE_OFFSET_EMPIRICAL){
			IO_utils.printLineErr("Adjusting footprint frame predictions using learned read-length vs CDS offset parameters:");
			globalOffsets = _readLengthWeights_empirical;
			printReadLengths2offsets(globalOffsets);
			//IO_utils.printLineErr("");
		}

		BufferedWriter out = new BufferedWriter(new FileWriter(outputPath));
		out.write("GeneID\tTranscriptID\tTranscriptBiotype\tStrand\tChromosome\tTotalNumberOfReads\tReadsWrongSize\tReadsOutsideCDS\tEligibleReads\tReadsInCorrectFrame\tReadsInWrongFrame\tFractionReadsInCorrectFrame\tmoderatedReadsInCorrectFrame\tmoderatedReadsInWrongFrame\tFractionModeratedReadsInCorrectFrame\n");

		// for progress bar
		int total = _annotation.getMap_gene2transcript().size();
		int counter_gene = 0;
		int counter_transcripts = 0;

		int lastVal = 0;
		IO_utils.printProgressBar(counter_gene);


		// loop through all genes
		Iterator<String> genes = _annotation.getMap_gene2transcript().keySet().iterator();
		String thisGeneID;
		while(genes.hasNext()){
			thisGeneID = genes.next();

			// progress bar
			counter_gene ++;
			if((int)Math.round((counter_gene+0.0)*100.0/(total+0.0)) > lastVal){
				lastVal = (int)Math.round((counter_gene+0.0)*100.0/(total+0.0));
				IO_utils.printProgressBar(lastVal);
			}

			// set up a map to count the number of reads for each transcript
			HashMap<String, Integer> transcript2readCount = new HashMap<String, Integer>(); 
			int maxCount = 0;
			String bestTranscript = null;

			// set up a map to record the read offsets each transcript
			// <transcriptID, <read, offset>>
			HashMap<String, HashMap<SAMRecord, Double>> transcript2readOffset = new HashMap<String, HashMap<SAMRecord, Double>>();

			// loop through all transcripts of this gene
			Iterator<String> transcripts = _annotation.getMap_gene2transcript().get(thisGeneID).iterator();
			while(transcripts.hasNext()){
				counter_transcripts ++;

				Transcript thisTranscript = _annotation.getTranscript(transcripts.next());

				// if this is a protein coding transcript:
				if(thisTranscript.getTranscriptBiotype().equals("protein_coding")){

					// if this transcript has reads aligned to it:
					if(_transcript2readAlignments.containsKey(thisTranscript.getID())){

						// calculate the read-mid offsets from the annotated coding frame of this transcript and add to the list for this gene:
						HashMap<SAMRecord, Double> readOffsets = calculateOffsetForReadsInThisTranscript(thisTranscript);
						transcript2readOffset.put(thisTranscript.getID(), readOffsets);

						// record the number of valid CDS reads mapping to this transcript
						transcript2readCount.put(thisTranscript.getID(), readOffsets.size());
						if(readOffsets.size() > maxCount){
							maxCount = readOffsets.size();
							bestTranscript = thisTranscript.getID();
						}

						// write this transcript to the log
						int[] counts = getReadFrameCounts(thisTranscript, readOffsets, globalOffsets);
						double fracCorrect = -1.0;
						if(counts[0]+counts[1] > 0)
							fracCorrect = (counts[0]+0.0)/(counts[0]+counts[1]+0.0);
						double fracCorrect_moderated = -1.0;
						if(counts[4]+counts[5] > 0)
							fracCorrect_moderated = (counts[4]+0.0)/(counts[4]+counts[5]+0.0);
						out.write(thisGeneID+"\t"+thisTranscript.getID()+"\t"+thisTranscript.getTranscriptBiotype()+"\t"+thisTranscript.getStrand()+"\t"+thisTranscript.getChromosome()+"\t"+
								readOffsets.size()+"\t"+counts[2]+"\t"+counts[3]+"\t"+(counts[0]+counts[1])+"\t"+counts[0]+"\t"+counts[1]+"\t"+fracCorrect+"\t"+counts[4]+"\t"+counts[5]+"\t"+fracCorrect_moderated+"\n");

					}
				}
			}


			// use the stats from the best transcript
			if(bestTranscript != null  &&  maxCount >= _minFootprints){//  &&  _annotation.getTranscript(bestTranscript).getStrand().equals("+")){
				HashMap<SAMRecord, Double> tmp_stats = transcript2readOffset.get(bestTranscript);
				Iterator<SAMRecord> it = tmp_stats.keySet().iterator();

				// to store read frame guesses:
				int correctFrame = 0;
				int wrongFrame = 0;
				int ignored_wrongLength = 0;

				// loop through all reads
				while(it.hasNext()){
					SAMRecord thisRead = it.next();
					double thisReadOffset = tmp_stats.get(thisRead);
					// if this read is within the CDS add to the global tally
					if(thisReadOffset >= 0.0){
						int thisReadLength = thisRead.getCigar().getReferenceLength();

						// Store read length vs. offset in the global result list
						if(!_readLengths2offsets.containsKey(thisReadLength))
							_readLengths2offsets.put(thisReadLength, new HashMap<Double, Integer>());
						if(!_readLengths2offsets.get(thisReadLength).containsKey(thisReadOffset))
							_readLengths2offsets.get(thisReadLength).put(thisReadOffset, 0);
						_readLengths2offsets.get(thisReadLength).put(thisReadOffset, _readLengths2offsets.get(thisReadLength).get(thisReadOffset) + 1);

						// can these reads guess the correct frame for this transcript after compensating for read length?
						if(getReadOffset(thisReadLength, globalOffsets) > -9.0){
							if(thisReadOffset == getReadOffset(thisReadLength, globalOffsets))
								correctFrame ++;
							else
								wrongFrame ++;
						}else{
							ignored_wrongLength ++;
						}
					}
				}
				_transcriptFrameGuesses.put(bestTranscript, new Integer[]{correctFrame, wrongFrame, ignored_wrongLength});
			}

		}
		out.flush();
		out.close();
		IO_utils.printLineErr("Done - "+counter_gene+" genes, "+counter_transcripts+" transcripts");
	}

	//
	private HashMap<String, Integer[]> _transcriptFrameGuesses = new HashMap<String, Integer[]>();
	/**
	 * 
	 * @param path
	 * @throws IOException 
	 */
	public void writeTranscriptFrameGuesses(String outputPath) throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputPath));
		out.write("TranscriptID\tStrand\tReadsInCorrectFrame\tReadsInWrongFrame\tReadsWrongSize\n");
		Iterator<String> transcripts = _transcriptFrameGuesses.keySet().iterator();
		while(transcripts.hasNext()){
			String thisTranscript = transcripts.next();
			Integer[] counts = _transcriptFrameGuesses.get(thisTranscript);
			out.write(thisTranscript+"\t"+_annotation.getTranscript(thisTranscript).getStrand()+"\t"+counts[0]+"\t"+counts[1]+"\t"+counts[2]+"\n");
		}
		out.flush();
		out.close();
	}













	private int[] getReadFrameCounts(Transcript thisTranscript, HashMap<SAMRecord, Double> readOffsets, HashMap<Integer, Double[]> globalOffsets){
		// to store read frame guesses:
		int correctFrame = 0, wrongFrame = 0;
		int ignored_wrongLength = 0, ignored_notInCDS = 0;

		//
		int moderated_correctFrame = 0, moderated_wrongFrame = 0;
		double tmp = 0.0;

		Iterator<SAMRecord> it = readOffsets.keySet().iterator();
		while(it.hasNext()){
			SAMRecord thisRead = it.next();
			int readLength = thisRead.getCigar().getReferenceLength();
			if(readOffsets.get(thisRead) >= 0.0){
				// can these reads guess the correct frame for this transcript after compensating for read length?
				if(getReadOffset(readLength, globalOffsets) > -9.0){
					if(readOffsets.get(thisRead) == getReadOffset(readLength, globalOffsets)){
						correctFrame ++;
						tmp = globalOffsets.get(readLength)[1];
						if(tmp > -9.0)
							moderated_correctFrame += Math.round(100*tmp);
					}else{
						wrongFrame ++;
						tmp = globalOffsets.get(readLength)[1];
						if(tmp > -9.0)
							moderated_wrongFrame += Math.round(100*tmp);
					}
				}else{
					ignored_wrongLength ++;
				}
			}else{
				ignored_notInCDS ++;
			}
		}
		return new int[]{correctFrame, wrongFrame, ignored_wrongLength, ignored_notInCDS, moderated_correctFrame, moderated_wrongFrame};
	}


	/**
	 * 
	 * @param readLength
	 * @return
	 */
	public double getLengthAdjustedReadMidPoint(int readLength, HashMap<Integer, Double[]> offsets){
		double readMid = (readLength+0.0)/2.0;
		double offset = getReadOffset(readLength, offsets);
		if(offset > -9.0)
			return readMid-offset;
		else
			return -999.0;
	}

	/**
	 * 
	 * @param readLength
	 * @return
	 */
	public double getReadOffset(int readLength, HashMap<Integer, Double[]> offsets){
		if(offsets.containsKey(readLength))
			return offsets.get(readLength)[0];
		else
			return -999.0;
	}










	/**
	 * 
	 * @param outputSummaryFile
	 * @throws IOException
	 */
	/*	public void calculateFrameConsistency(File outputSummaryFile) throws IOException {
		IO_utils.printLineErr("Calculating frame consistency for footprints of each transcript...");

		BufferedWriter outputFile = null;

		try {

			outputFile = new BufferedWriter(new FileWriter(outputSummaryFile));

			// large array of individual transcript summaries that we're eventually going to write to a file
			ArrayList<TranscriptSummary> transcriptSummaries = new ArrayList<TranscriptSummary>();

			this.numberOfTests = _transcript2readAlignments.size();
			this.adjustedPValue = 0.05 / this.numberOfTests; 

			minimumNumberOfReads = (int) Math.ceil(Math.log(adjustedPValue) / Math.log(1/3d)) + 1;

			IO_utils.printLineErr("Based on "+this.numberOfTests+" tests, any transcript with fewer than "+minimumNumberOfReads+" reads cannot achieve a statistically significant frame consistency");


			Iterator<Entry<String, ArrayList<SAMRecord>>> it = _transcript2readAlignments.entrySet().iterator();
			int total = _transcript2readAlignments.size();
			int counter = 0;
			int lastVal = 0;
			IO_utils.printProgressBar(counter);


			String fileHeader = "\"transcript_id\"\t\"transcript_type\"\t\"total_reads\"\t\"max_aligned\"\t\"mid_aligned\"\t\"min_aligned\"\t\"alignment_score\"\t\"significant?\"";
			outputFile.write(fileHeader);
			outputFile.newLine();

			while (it.hasNext()) {

				if(Math.round((counter+0.0)*100.0/(total+0.0)) > lastVal){
					lastVal = (int)Math.round((counter+0.0)*100.0/(total+0.0));
					IO_utils.printProgressBar(lastVal);
				}

				Entry<String, ArrayList<SAMRecord>> thisTranscript = it.next();

				String transcriptID = thisTranscript.getKey();
				ArrayList<SAMRecord> reads = thisTranscript.getValue();
				int totalReads = reads.size();

				HashMap<Integer, Integer> readFrames = new HashMap<Integer, Integer>();
				readFrames.put(0, 0);
				readFrames.put(1, 0);
				readFrames.put(2, 0);
				int readPosition;
				for (int i = 0; i < totalReads; i++) {
					readPosition = reads.get(i).getAlignmentStart() % 3;
					readFrames.put(readPosition, readFrames.get(readPosition) + 1);
				}
				int bin0 = readFrames.get(0);
				int bin1 = readFrames.get(1);
				int bin2 = readFrames.get(2);
				int maxAligned = Math.max(Math.max(bin0, bin1), bin2);
				int minAligned = Math.min(Math.min(bin0, bin1), bin2);
				int midAligned = (bin0 + bin1 + bin2) - maxAligned - minAligned;

				//String transcriptType = transcriptTypeMap.get(transcriptID);
				String transcriptType = _annotation.getTranscript(transcriptID).getTranscriptBiotype();


				TranscriptSummary thisSummary = new TranscriptSummary(transcriptID, transcriptType, totalReads, maxAligned, midAligned, minAligned, adjustedPValue);
				transcriptSummaries.add(thisSummary);

				thisSummary.writeSummaryToFile(outputFile);

				counter ++;
			}

		} catch (IOException ex) {
			Logger.getLogger(FootprintFrameCalculator.class.getName()).log(Level.SEVERE, null, ex);
		} finally {
			if (outputFile != null) outputFile.close();

			IO_utils.printLineErr("All done.");
			//System.out.println("You finished, yay!!!!");
		}
	}


	public static BigInteger numberOfCombinations(int stringLength, int k) {

		try {
			if (k < 1 || k > stringLength) {
				throw new IllegalArgumentException();
			}
		} catch (IllegalArgumentException e) {
			System.out.println("Illegal argument exception: " + e.getMessage());
			System.out.println("You cannot find the combinations of k symbols when k is above the given string length or below 1");
			System.exit(1);
		}

		// numerator factorial
		BigInteger numPermutations = new BigInteger("1");

		// implement classic permutations with `i` number of slots left over
		for (int i = stringLength; i > stringLength - k; i--) {
			numPermutations = numPermutations.multiply(new BigInteger(i + ""));
		}
		BigInteger kPrime = new BigInteger("1");
		for (int i = k; i > 0; i--) {
			kPrime = kPrime.multiply(new BigInteger(i + "")); 
		}

		return numPermutations.divide(kPrime);

	}

	public static BigInteger factorial(int n) {
		BigInteger result = BigInteger.ONE;
		while (n != 0) {
			result = result.multiply(BigInteger.valueOf(n));
			n--;
		}
		return result;
	} 
	 */



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		//options.addOption(new Option("verbose", "write the rejected sequences to stderr"));
		options.addOption(OptionBuilder.withArgName("annotationPath").hasArg().withDescription("Path to the GTF file containing the transcript/gene relationship").create(Thunder.OPT_PATH_ANNOTATION));
		options.addOption(OptionBuilder.withArgName("footprint alignments").hasArg().withDescription("SAM/BAM file containing alignments against the transcriptome").create("f"));
		options.addOption(OptionBuilder.withLongOpt("outputPrefix").withArgName("outputPath").hasArg().withDescription("[optional] Output prefix for results [default: same as input file]").create("o"));
		options.addOption(OptionBuilder.withArgName("minFootprints").hasArg().withDescription("[optional] Transcript must have more than this minimum number of footprints to be included [default: 1]").create("N"));
		options.addOption(OptionBuilder.withArgName("maxReads").hasArg().withDescription("[optional] Only read the first R reads from the input bam [default: 10000000]").create("R"));
		return options;
	}



	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		/*args = new String[]{"-f","/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Footprinting/FP_original_1_Genome_Aligned.toTranscriptome.sorted.bam",
				"-a","/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation_noSelenocysteine.gtf",
				"-o","/Users/robk/Downloads/FootprintFrameAnalysis_TEST_",
				"-N","3"};
		 */

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		//if(cmdArgs.hasOption(Thunder.OPT_PATH_ANNOTATION)  &&  cmdArgs.hasOption("f")  &&  cmdArgs.hasOption("o")){
		if(cmdArgs.hasOption(Thunder.OPT_PATH_ANNOTATION)  &&  cmdArgs.hasOption("f")){
			System.err.println();

			String outputPrefix = (new File(cmdArgs.getOptionValue("f"))).getAbsolutePath();
			if(cmdArgs.hasOption("o")){
				outputPrefix = cmdArgs.getOptionValue("o");
			}

			// check override for min # of footprints
			int minFP = 1;
			if(cmdArgs.hasOption("N"))
				minFP = Integer.valueOf(cmdArgs.getOptionValue("N"));

			// start and read the annotation information
			FootprintFrameCalculator engine = new FootprintFrameCalculator(cmdArgs.getOptionValue(Thunder.OPT_PATH_ANNOTATION), minFP);

			// check override for min # of footprints
			int maxReads = 10000000;
			if(cmdArgs.hasOption("R"))
				maxReads = Integer.valueOf(cmdArgs.getOptionValue("R"));
			// Read the footprint alignments
			engine.readAlignments(new File(cmdArgs.getOptionValue("f")), maxReads);

			// Compare footprint midpoints to annotated codon positions
			engine.computeReadMidOffset2CDS(outputPrefix+"_FramesPerTranscript_guessed.txt", _USE_OFFSET_GUESSED, false);
			engine.writeReadLengths2offsets(outputPrefix+"_ReadLength_vs_codonOffset.txt");
			//engine.writeTranscriptFrameGuesses(outputPrefix+"_TranscriptFrames_guessed.txt");

			engine.computeReadMidOffset2CDS(outputPrefix+"_FramesPerTranscript_learned.txt", _USE_OFFSET_EMPIRICAL, false);
			//engine.writeTranscriptFrameGuesses(outputPrefix+"_TranscriptFrames_learned.txt");


			// Calculate frame-consistency of reads for each transcript (to be used in the EM) 
			//engine.calculateFrameConsistency(new File(cmdArgs.getOptionValue("o")+"FrameConsistencyPerTranscript.txt"));
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" FootprintFrameAnaysis", getCmdLineOptions());
			System.err.println();
		}
	}

}
