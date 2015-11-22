package footprintAlignments;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.math.BigDecimal;
import java.math.BigInteger;
import java.math.RoundingMode;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;

import main.Thunder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import objects.Transcript;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import utils.IO_utils;
import annotation.ReadGTF;
import annotation.TranscriptAnnotation;


/**
 *
 * @author J. David McPeek
 */
public class FootprintFrameCalculator {
	private TranscriptAnnotation _annotation;
	private int _minFootprints = 1;
	
	int numberOfTests;  // the number of transcripts (hypotheses). This will be used to determine our new p-value adjusted multiple testing correction. 
	double adjustedPValue;
	int minimumNumberOfReads;
	private static final HashMap<Integer, BigInteger> factorialCache = new HashMap<Integer, BigInteger>();  // don't redo computations


	// TODO: perhaps have some other interfaces if clients already have the gencode data
	public FootprintFrameCalculator(String annotation, int minFootprints) throws Exception {
		_annotation = ReadGTF.readGTF(annotation);
		_minFootprints = minFootprints;
		/*System.out.println("annotation.getTranscripts().keySet().size() = "+_annotation.getTranscripts().keySet().size());
		Iterator<String> it = _annotation.getTranscripts().keySet().iterator();
		while(it.hasNext()){
			String thisTranscript = it.next();
			System.out.println(_annotation.getTranscripts().get(thisTranscript).toString());
		}*/
	}

	protected static class TranscriptSummary {
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


	}

	// stores transcriptome alignmnets
	private Map<String, ArrayList<SAMRecord>> _transcript2readAlignments = new HashMap<String, ArrayList<SAMRecord>>();;

	/**
	 * Reads transcriptome alignments in the given file
	 * @param alignmentFile
	 */
	private void readAlignments(File alignmentFile){
		int count = 0;

		SAMFileReader inputBam = new SAMFileReader(alignmentFile);

		IO_utils.printLineErr("Reading alignments...");
		for (final SAMRecord thisRead : inputBam) {
			count ++;

			//if(count < 100){

			// we're only looking at reads that match 29 base pairs
			//if (thisRead.getCigar().getReferenceLength() != 29) {
			//	continue;
			//}
			String transcriptID = thisRead.getReferenceName();
			if(!_transcript2readAlignments.containsKey(transcriptID)) {
				_transcript2readAlignments.put(transcriptID, new ArrayList<SAMRecord>());
			}
			_transcript2readAlignments.get(transcriptID).add(thisRead);

			//}
		}
		if (inputBam != null) inputBam.close();

		IO_utils.printLineErr("Done- read "+count+" alignments.");
	}

	private HashMap<Double, HashMap<Double, Integer>> _readLengths2offsets = new HashMap<Double, HashMap<Double, Integer>>();


	/**
	 * 
	 */
	public void printReadLengths2offsets(){
		// print the result:
		Iterator<Double> lengths = _readLengths2offsets.keySet().iterator();
		while(lengths.hasNext()){
			Double thisLength = lengths.next();
			Iterator<Double> offsets = _readLengths2offsets.get(thisLength).keySet().iterator();
			while(offsets.hasNext()){
				Double thisOffset = offsets.next();
				System.out.println(thisLength+"\t"+thisOffset+"\t"+_readLengths2offsets.get(thisLength).get(thisOffset));
			}
		}
	}

	/**
	 * 
	 * @param outputPath
	 * @throws IOException 
	 */
	public void writeReadLengths2offsets(String outputPath) throws IOException{
		BufferedWriter out = new BufferedWriter(new FileWriter(outputPath));
		out.write("ReadLength\tOffsetToCDS\tNumberOfReads\n");

		Iterator<Double> lengths = _readLengths2offsets.keySet().iterator();
		while(lengths.hasNext()){
			Double thisLength = lengths.next();
			Iterator<Double> offsets = _readLengths2offsets.get(thisLength).keySet().iterator();
			while(offsets.hasNext()){
				Double thisOffset = offsets.next();
				out.write(thisLength+"\t"+thisOffset+"\t"+_readLengths2offsets.get(thisLength).get(thisOffset)+"\n");
			}
		}
		out.flush();
		out.close();
	}

	/**
	 * 
	 */
	public void compareFootprintsToAnnotatedCDS(){
		IO_utils.printLineErr("Calculating frame call accuracy...");


		// loop through all genes
		Iterator<String> genes = _annotation.getMap_gene2transcript().keySet().iterator();
		String thisGeneID;
		while(genes.hasNext()){
			thisGeneID = genes.next();

			// set up a map to count the number of reads for each transcript
			HashMap<String, Integer> transcript2readCount = new HashMap<String, Integer>(); 
			int maxCount = 0;
			String bestTranscript = null;

			// set up a map to record the read offsets each transcript
			HashMap<String, Double[][]> transcript2readOffset = new HashMap<String, Double[][]>();

			// loop through all transcripts of this gene
			Iterator<String> transcripts = _annotation.getMap_gene2transcript().get(thisGeneID).iterator();
			while(transcripts.hasNext()){
				Transcript thisTranscript = _annotation.getTranscript(transcripts.next());

				// if this is a protein coding transcript:
				if(thisTranscript.getTranscriptBiotype().equals("protein_coding")){

					// if this transcript has reads aligned to it:
					if(_transcript2readAlignments.containsKey(thisTranscript.getID())){
						ArrayList<SAMRecord> reads = _transcript2readAlignments.get(thisTranscript.getID());

						// reset the counter of valid (CDS) reads for this transcript
						int allReadCount = 0;
						int validReadCount = 0;
						transcript2readOffset.put(thisTranscript.getID(), new Double[reads.size()][2]);

						// compute the CDS offset from the beginning/end of this transcript (to be able to pick CDS reads and to compute distance from start codon)
						int cdsOffset_start = thisTranscript.getCDSStartFromTxStartInTxCoords();
						int cdsOffset_end = thisTranscript.getCDSStopFromTxStopInTxCoords();

						// loop through each read
						Iterator<SAMRecord> readsIterator = reads.iterator();
						while(readsIterator.hasNext()){

							SAMRecord thisRead = readsIterator.next();

							double readLength = thisRead.getCigar().getReferenceLength()+0.0;
							transcript2readOffset.get(thisTranscript.getID())[allReadCount][0] = readLength;

							// if the middle of this read is within the CDS
							double readMid = thisRead.getAlignmentStart()+0.0 + (readLength/2.0);
							if(readMid >= cdsOffset_start  &&  readMid >= cdsOffset_end){
								// add to the read count for this transcript
								validReadCount ++;

								// record the offset of this read from the start of the CDS
								transcript2readOffset.get(thisTranscript.getID())[allReadCount][1] = (readMid - cdsOffset_start) % 3;
							}else{
								transcript2readOffset.get(thisTranscript.getID())[allReadCount][1] = -999.0;
							}
							allReadCount ++;
						}

						// record the number of valid CDS reads mapping to this transcript
						transcript2readCount.put(thisTranscript.getID(), validReadCount);
						if(validReadCount > maxCount){
							maxCount = validReadCount;
							bestTranscript = thisTranscript.getID();
						}
					}
				}
			}

			// use the stats from the best transcript
			if(bestTranscript != null  &&  maxCount >= _minFootprints){
				Double[][] tmp_stats = transcript2readOffset.get(bestTranscript);
				// loop through all reads
				for(int i=0;i<tmp_stats.length;i++){
					// if this read is within the CDS add to the global tally
					if(tmp_stats[i][1] >= 0.0){
						if(!_readLengths2offsets.containsKey(tmp_stats[i][0]))
							_readLengths2offsets.put(tmp_stats[i][0], new HashMap<Double, Integer>());

						if(!_readLengths2offsets.get(tmp_stats[i][0]).containsKey(tmp_stats[i][1]))
							_readLengths2offsets.get(tmp_stats[i][0]).put(tmp_stats[i][1], 0);

						_readLengths2offsets.get(tmp_stats[i][0]).put(tmp_stats[i][1], _readLengths2offsets.get(tmp_stats[i][0]).get(tmp_stats[i][1]) + 1);
					}
				}
			}
		}
		IO_utils.printLineErr("Done.");
	}




	public void calculateFrameConsistency(File outputSummaryFile) throws IOException {
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



			String fileHeader = "\"transcript_id\"\t\"transcript_type\"\t\"total_reads\"\t\"max_aligned\"\t\"mid_aligned\"\t\"min_aligned\"\t\"alignment_score\"\t\"significant?\"";
			outputFile.write(fileHeader);
			outputFile.newLine();

			while (it.hasNext()) {

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




	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(new Option("verbose", "write the rejected sequences to stderr"));
		options.addOption(OptionBuilder.withArgName("annotationPath").hasArg().withDescription("Path to the GTF file containing the transcript/gene relationship").create(Thunder.OPT_PATH_ANNOTATION));
		options.addOption(OptionBuilder.withArgName("footprint alignments").hasArg().withDescription("SAM/BAM file containing alignments against the transcriptome").create("f"));
		options.addOption(OptionBuilder.withLongOpt("outputPrefix").withArgName("outputPath").hasArg().withDescription("Output prefix for results").create("o"));
		options.addOption(OptionBuilder.withArgName("minFootprints").hasArg().withDescription("Transcript must have more than this minimum number of footprints to be included [default: 1]").create("N"));
		return options;
	}



	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		/*args = new String[]{"-f","/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Footprinting/FP_original_1_Genome_Aligned.toTranscriptome.sorted.bam",
							"-a","/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation_noSelenocysteine.gtf",
							"-o","/Users/robk/Downloads/FootprintFrameAnalysis_TEST_",
							"-N","1"};
		*/
		
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption(Thunder.OPT_PATH_ANNOTATION)  &&  cmdArgs.hasOption("f")  &&  cmdArgs.hasOption("o")){
			System.err.println();

			// start and read the annotation information
			FootprintFrameCalculator engine = new FootprintFrameCalculator(cmdArgs.getOptionValue(Thunder.OPT_PATH_ANNOTATION), Integer.valueOf(cmdArgs.getOptionValue("N")));

			// Read the footprint alignments
			engine.readAlignments(new File(cmdArgs.getOptionValue("f")));

			// Compare footprint midpoints to annotated codon positions
			engine.compareFootprintsToAnnotatedCDS();
			engine.writeReadLengths2offsets(cmdArgs.getOptionValue("o")+"FrameConsistencySummary.txt");

			// Calculate frame-consistency of reads for each transcript (to be used in the EM) 
			engine.calculateFrameConsistency(new File(cmdArgs.getOptionValue("o")+"FrameConsistencyPerTranscript.txt"));
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" FootprintFrameAnaysis", getCmdLineOptions());
			System.err.println();
		}
	}

}
