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
import objects.SAMRecordReduced;
import objects.Transcript;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import samTools.SAMReader;
import utils.IO_utils;
import annotation.ReadGTF;
import annotation.TranscriptAnnotation;



public class FootprintFrameCalculator {

	private boolean _printProgress = true;
	public void suppressProgressBar(){ _printProgress = false; }
	public boolean printProgressBar(){ return _printProgress; }

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



	//public String _printVerboseForGeneName = "ENSG00000213380.11";



	// stores transcriptome alignmnets
	private Map<String, ArrayList<SAMRecordReduced>> _transcript2readAlignments = new HashMap<String, ArrayList<SAMRecordReduced>>();

	/**
	 * Reads *transcriptome* alignments in the given file
	 * @param alignmentFile
	 */
	private void readAlignments(File alignmentFile, long maxReads){
		
		IO_utils.printLineErr("Reading alignments...");
		
		int count_allReads = 0;
		int count_uniqueReads = 0;
		//int count_allAlignments = 0;
		int count_uniqueAlignments = 0;

		//SAMFileReader inputBam = new SAMFileReader(alignmentFile);

		SAMReader engine = new SAMReader(alignmentFile);
		//if(verbose)
		//	System.err.println("BAM sorted by readID: "+engine.isSortedByReadID());

		if( ! engine.isSortedByReadID()){
			System.err.println("BAM must be sorted by readID ('queryname')\nThis BAM is sorted by: "+engine.getSortOrder().toString());
			System.exit(0);
		}


		ArrayList<SAMRecordReduced> thisRead;
		// Loop through all reads
		while((thisRead = engine.getAlignmentsForNextRead()) != null){

			count_allReads ++;
			String thisTranscriptID = "";

			// if this read maps to a single genic locus:
			if(SAMReader.isReadUniquelyMappedToOneGene(thisRead, _annotation.getMap_transcript2gene())){

				count_uniqueReads ++;
				boolean firstAlignment = true;

				// loop again through the alignments and add to the gene list
				for(SAMRecordReduced thisAlignment: thisRead){

					count_uniqueAlignments ++;
					
					// get the transcriptID for this alignment:
					thisTranscriptID = thisAlignment.getReferenceName();
					
					// if this read is in the CDS of this transcript:
					if(_annotation.getTranscript(thisTranscriptID).hasCDS()){
					
						double readMid = thisAlignment.getAlignmentStart()+((thisAlignment.getAlignmentEnd()-thisAlignment.getAlignmentStart()+0.0)/2.0);
						
						// Only allow alignments to the ORF:
						if(_annotation.getTranscript(thisTranscriptID).isCoordInCDS(readMid)){
					
							// add read to transcript 
							if(!_transcript2readAlignments.containsKey(thisTranscriptID)) {
								_transcript2readAlignments.put(thisTranscriptID, new ArrayList<SAMRecordReduced>());
							}
							_transcript2readAlignments.get(thisTranscriptID).add(thisAlignment);

							// add read length to the list:
							if(firstAlignment){
								int readLength = thisAlignment.getReadSequence().length();
								if(!_readLengths2readCount.containsKey(readLength))
									_readLengths2readCount.put(readLength,0);
								_readLengths2readCount.put(readLength, _readLengths2readCount.get(readLength) + 1);
								firstAlignment = false;
							}
						}
					}
				}
			}
			
			if(maxReads > 0  &&  count_uniqueReads >= maxReads){
				IO_utils.printLineErr("** WARNING **  Using only the first "+maxReads+" reads (use -R -1 to use ALL reads).");
				break;
			}
		}
		
		engine.close();

		IO_utils.printLineErr("Done - Using "+count_uniqueAlignments+" alignments of "+count_allReads+" reads.");
		IO_utils.printLineErr("     - Discarded "+(count_allReads-count_uniqueReads)+" reads ("+((Math.round((count_allReads-count_uniqueReads)*10000.0/count_allReads)/100.0))+"%) due to multi-gene ambiguity.");

	}



	// HashMap<readLength, HashMap<CDSoffset, nReads>>
	private HashMap<Integer, HashMap<Double, Integer>> _readLength2offset2readCount = new HashMap<Integer, HashMap<Double, Integer>>();
	private HashMap<Integer, Integer> _readLengths2readCount = new HashMap<Integer, Integer>();

	//HashMap<readLength, HashMap<CDSoffset, weight>>
	private HashMap<Integer, HashMap<Double, Double>> _readLength2offset2weight = new HashMap<Integer, HashMap<Double, Double>>();
	public HashMap<Integer, HashMap<Double, Double>> getReadLengths2offsets2weights(){ return _readLength2offset2weight; }
	//public HashMap<Integer, HashMap<Double, Integer>> getReadLengths2offsets(){ return _readLength2offset2readCount; }


	/**
	 * @throws IOException 
	 * 
	 */
	public void printReadLengths2offsets(HashMap<Integer,Double[]> globalOffsets) throws IOException{
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

		// print and write the result:
		//BufferedWriter out = new BufferedWriter(new FileWriter(outputPath+"_ReadLength2Offset2Weight.txt"));
		//out.write("ReadLength\tReadCount\tOffsetToCDS\tWeight\n");
		for(int i=minLength;i<=maxLength;i++){
			int readCount = 0;
			if(_readLengths2readCount.containsKey(i))
				readCount = _readLengths2readCount.get(i);

			if(globalOffsets.containsKey(i)){
				double weight = (Math.round(100.0*globalOffsets.get(i)[1])/100.0);
				if(weight >= 0.5)
					IO_utils.printLineErr("       Read-length:"+i+"  read-count:"+readCount+"  CDS-offset:"+globalOffsets.get(i)[0]+"  weight:"+weight);
				//out.write(i+"\t"+readCount+"\t"+globalOffsets.get(i)[0]+"\t"+weight+"\n");
			}
		}
		//out.flush();
		//out.close();
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
		Iterator<Integer> lengths = _readLength2offset2readCount.keySet().iterator();
		while(lengths.hasNext()){
			int thisLength = lengths.next();
			if(thisLength < minLength)
				minLength = thisLength;
			else if(thisLength > maxLength)
				maxLength = thisLength;
		}

		BufferedWriter out = new BufferedWriter(new FileWriter(outputPath));
		out.write("ReadLength\tOffsetToCDS\tNumberOfReads\tFractionOfReadsThisLength\n");

		for(int i=minLength;i<=maxLength;i++){
			int totalReads = 0;
			int bestReads = 0;
			double bestOffset = 0.0;

			if(_readLength2offset2readCount.containsKey(i)){

				// find the total number of reads of this length:
				Iterator<Double> offsets = _readLength2offset2readCount.get(i).keySet().iterator();
				while(offsets.hasNext()){
					Double thisOffset = offsets.next();
					totalReads += _readLength2offset2readCount.get(i).get(thisOffset);
				}

				// print the length/offset/count/fraction and compute the best offset
				offsets = _readLength2offset2readCount.get(i).keySet().iterator();
				while(offsets.hasNext()){
					Double thisOffset = offsets.next();
					int nReads = _readLength2offset2readCount.get(i).get(thisOffset);

					_readLength2offset2weight.get(i).put(thisOffset, (nReads+0.0)/(totalReads+0.0));
					out.write(i+"\t"+thisOffset+"\t"+nReads+"\t"+((nReads+0.0)/(totalReads+0.0))+"\n");

					if(nReads > bestReads){
						bestReads = nReads;
						bestOffset = thisOffset;
					}
				}

				double fractionReads = (bestReads+0.0)/(totalReads+0.0);
				//if(thisLength >= 26  &&  thisLength <= 33  &&  fractionReads > 0.5){
				//if(fractionReads > 0.5){
				if(i >= 10  &&  i <= 45  &&  fractionReads > 0.33){
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
	private HashMap<SAMRecordReduced, Double> calculateOffsetForReadsInThisTranscript(Transcript thisTranscript){
		ArrayList<SAMRecordReduced> reads = _transcript2readAlignments.get(thisTranscript.getID());
		HashMap<SAMRecordReduced, Double> readOffsets = new HashMap<SAMRecordReduced, Double>();

		// compute the CDS offset from the beginning/end of this transcript (to be able to pick CDS reads and to compute distance from start codon)
		int cdsOffset_start = thisTranscript.getCDSStartFromTxStartInTxCoords();
		//int cdsOffset_end = thisTranscript.getCDSStopFromTxStopInTxCoords();
		//int	cdsLength = thisTranscript.getTotalCodingExonLength();
		
		//if(thisTranscript.getGeneID().equals(_printVerboseForGeneName))
		//	System.out.println(thisTranscript.getGeneID()+"\t"+thisTranscript.getID()+"\t"+cdsLength+"\t"+cdsOffset_start);

		/*if(thisTranscript.getID().equals("ENST00000331710.7")){
			System.out.print(thisTranscript.getID());
			System.out.println("\t"+cdsOffset_start);
		}*/

		// loop through each read
		Iterator<SAMRecordReduced> readsIterator = reads.iterator();
		while(readsIterator.hasNext()){
			SAMRecordReduced thisRead = readsIterator.next();

			// read strand must be +ve w.r.t the transcript!
			boolean negativeStrand = thisRead.getReadNegativeStrandFlag();

			// if the middle of this read is within the CDS
			double readMid = thisRead.getAlignmentStart()+0.0 + ((thisRead.getCigar().getReferenceLength()+0.0)/2.0);

			//if(thisTranscript.getGeneID().equals(_printVerboseForGeneName))
			//	System.out.println("\t"+thisRead.getReadName()+"\t"+negativeStrand+"\t"+readMid+"\t"+thisTranscript.isCoordInCDS(readMid));

			//if(readMid >= cdsOffset_start  &&  readMid >= cdsOffset_end  &&  !negativeStrand){
			if(!negativeStrand  &&  thisTranscript.isCoordInCDS(readMid)){
				// record the offset of this read from the start of the CDS
				readOffsets.put(thisRead, (readMid - cdsOffset_start) % 3);
			}else{
				readOffsets.put(thisRead, -999.0);
			}
		}
		return readOffsets;
	}


	private static final int _USE_OFFSET_GUESSED = 1;
	private static final int _USE_OFFSET_EMPIRICAL = 2;


	/**
	 * @throws IOException 
	 * 
	 */
	public void computeReadMidOffset2CDS(String outputPath, int useOffsets, boolean useFrameToSelectBestTranscript, boolean outputRead2TranscriptFrameMatrix) throws IOException{

		BufferedWriter out = null;
		HashMap<Integer, Double[]> globalOffsets = _readLengthWeights_guessed;
		if(useOffsets == _USE_OFFSET_GUESSED){
			//IO_utils.printLineErr("Learning global read-length vs CDS offset parameters using guessed read-length frame offsets:");
			IO_utils.printLineErr("Learning global read-length vs CDS offset parameters...");
			//printReadLengths2offsets(globalOffsets);
			//IO_utils.printLineErr("");
			out = new BufferedWriter(new FileWriter(outputPath+"_FramesPerTranscript_guessed.txt"));
		}else if(useOffsets == _USE_OFFSET_EMPIRICAL){
			IO_utils.printLineErr("Adjusting footprint frame predictions using learned read-length vs CDS offset parameters:");
			globalOffsets = _readLengthWeights_empirical;
			printReadLengths2offsets(globalOffsets);
			//IO_utils.printLineErr("");
			out = new BufferedWriter(new FileWriter(outputPath+"_FramesPerTranscript_learned.txt"));
		}

		out.write("GeneID\tTranscriptID\tTranscriptBiotype\tStrand\tChromosome\tTotalNumberOfReads\tReadsWrongSize\tReadsOutsideCDS\tEligibleReads\tReadsInCorrectFrame\tReadsInWrongFrame\tFractionReadsInCorrectFrame\tmoderatedReadsInCorrectFrame\tmoderatedReadsInWrongFrame\tFractionModeratedReadsInCorrectFrame\n");

		BufferedWriter readInFrameWithTranscript = null;
		if(outputRead2TranscriptFrameMatrix){
			readInFrameWithTranscript = new BufferedWriter(new FileWriter(outputPath+"_Read2TranscriptFrame.txt"));
			readInFrameWithTranscript.write("GeneID\tTranscriptID\tReadID\tReadLength\tReadSenseToTranscript\tReadMidFrameOffset\tReadInCDS\tReadLikelyInFrame\n");
		}

		// for progress bar
		int total = _annotation.getMap_gene2transcript().size();
		int counter_gene = 0;
		int counter_transcripts = 0;

		int lastVal = 0;
		if(printProgressBar())
			IO_utils.printProgressBar(counter_gene);


		//
		// loop through all genes
		//
		Iterator<String> genes = _annotation.getMap_gene2transcript().keySet().iterator();
		String thisGeneID;
		while(genes.hasNext()){
			thisGeneID = genes.next();

			// progress bar
			counter_gene ++;
			if(printProgressBar()){
				if((int)Math.round((counter_gene+0.0)*100.0/(total+0.0)) > lastVal){
					lastVal = (int)Math.round((counter_gene+0.0)*100.0/(total+0.0));
					IO_utils.printProgressBar(lastVal);
				}
			}

			// set up a map to count the number of reads for each transcript
			HashMap<String, Integer> transcript2readCount = new HashMap<String, Integer>(); 
			int maxCount = 0;
			String bestTranscript = null;

			// set up a map to record the read offsets each transcript
			// <transcriptID, <read, offset>>
			HashMap<String, HashMap<SAMRecordReduced, Double>> transcript2readOffset = new HashMap<String, HashMap<SAMRecordReduced, Double>>();

			// loop through all transcripts of this gene
			Iterator<String> transcripts = _annotation.getMap_gene2transcript().get(thisGeneID).iterator();
			while(transcripts.hasNext()){
				counter_transcripts ++;

				Transcript thisTranscript = _annotation.getTranscript(transcripts.next());

				//if(thisGeneID.equals(_printVerboseForGeneName))
				//	System.out.print(thisGeneID+"\t"+thisTranscript.getID());


				// if this is a protein coding transcript:
				if(thisTranscript.getTranscriptBiotype().equals("protein_coding")){

					// if this transcript has reads aligned to it:
					if(_transcript2readAlignments.containsKey(thisTranscript.getID())){

						// calculate the read-mid offsets from the annotated coding frame of this transcript and add to the list for this gene:
						HashMap<SAMRecordReduced, Double> readOffsets = calculateOffsetForReadsInThisTranscript(thisTranscript);
						transcript2readOffset.put(thisTranscript.getID(), readOffsets);

						// record the number of valid CDS reads mapping to this transcript
						transcript2readCount.put(thisTranscript.getID(), readOffsets.size());
						if(readOffsets.size() > maxCount){
							maxCount = readOffsets.size();
							bestTranscript = thisTranscript.getID();
						}

						// write this transcript to the log
						//int[] counts = getReadFrameCounts(thisTranscript, readOffsets, globalOffsets);
						int[] counts = getReadFrameCounts(thisTranscript, readOffsets, globalOffsets, readInFrameWithTranscript);

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
				HashMap<SAMRecordReduced, Double> tmp_stats = transcript2readOffset.get(bestTranscript);
				Iterator<SAMRecordReduced> it = tmp_stats.keySet().iterator();

				// to store read frame guesses:
				int correctFrame = 0;
				int wrongFrame = 0;
				int ignored_wrongLength = 0;

				// loop through all reads
				while(it.hasNext()){
					SAMRecordReduced thisRead = it.next();
					double thisReadOffset = tmp_stats.get(thisRead);
					// if this read is within the CDS add to the global tally
					if(thisReadOffset >= 0.0){
						int thisReadLength = thisRead.getCigar().getReferenceLength();

						// Store read length vs. offset in the global result list
						if(!_readLength2offset2readCount.containsKey(thisReadLength)){
							_readLength2offset2readCount.put(thisReadLength, new HashMap<Double, Integer>());
							_readLength2offset2weight.put(thisReadLength, new HashMap<Double, Double>());
						}

						if(!_readLength2offset2readCount.get(thisReadLength).containsKey(thisReadOffset)){
							_readLength2offset2readCount.get(thisReadLength).put(thisReadOffset, 0);
							_readLength2offset2weight.get(thisReadLength).put(thisReadOffset, 0.0);
						}
						_readLength2offset2readCount.get(thisReadLength).put(thisReadOffset, _readLength2offset2readCount.get(thisReadLength).get(thisReadOffset) + 1);

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


		if(outputRead2TranscriptFrameMatrix){
			readInFrameWithTranscript.flush();
			readInFrameWithTranscript.close();
		}



		out.flush();
		out.close();
		
		if(printProgressBar())
			IO_utils.printErr("\n");
		
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


	/**
	 * 
	 * @param thisTranscript
	 * @param readOffsets
	 * @param globalOffsets
	 * @return
	 * @throws IOException 
	 */
	private int[] getReadFrameCounts(Transcript thisTranscript, HashMap<SAMRecordReduced, Double> readOffsets, HashMap<Integer, Double[]> globalOffsets, BufferedWriter readInFrameWithTranscript) throws IOException{
		// to store read frame guesses:
		int correctFrame = 0, wrongFrame = 0;
		int ignored_wrongLength = 0, ignored_notInCDS = 0;

		//
		int moderated_correctFrame = 0, moderated_wrongFrame = 0;

		Double[] thisReadLength2Weight;
		double tmp = 0.0;

		Iterator<SAMRecordReduced> it = readOffsets.keySet().iterator();
		while(it.hasNext()){
			SAMRecordReduced thisRead = it.next();
			int readLength = thisRead.getCigar().getReadLength();
			double thisReadOffset = readOffsets.get(thisRead);
			boolean readSense = !thisRead.getReadNegativeStrandFlag();

			// write this read/transcript in frame?
			if(readInFrameWithTranscript != null)
				readInFrameWithTranscript.write(thisTranscript.getGeneID()+"\t"+thisTranscript.getID()+"\t"+thisRead.getReadName()+"\t"+readLength+"\t"+readSense);
			//readInFrameWithTranscript.write("GeneID\tTranscriptID\tReadID\tReadLength\tReadMidFrameOffset\tReadInCDS\tReadLikelyInFrame\n");

			if(thisReadOffset >= 0.0){

				// can these reads guess the correct frame for this transcript after compensating for read length?
				if(getReadOffset(readLength, globalOffsets) > -9.0){

					// write this read/transcript in frame?
					if(readInFrameWithTranscript != null)
						readInFrameWithTranscript.write("\t"+thisReadOffset);

					if(thisReadOffset == getReadOffset(readLength, globalOffsets)){
						correctFrame ++;
						thisReadLength2Weight = globalOffsets.get(readLength);
						tmp = thisReadLength2Weight[1];
						if(tmp > -9.0)
							moderated_correctFrame += Math.round(100*tmp);

						// write this read/transcript in frame?
						if(readInFrameWithTranscript != null)
							readInFrameWithTranscript.write("\t"+"true"+"\t"+"1"+"\n");

					}else{
						wrongFrame ++;
						tmp = globalOffsets.get(readLength)[1];
						if(tmp > -9.0)
							moderated_wrongFrame += Math.round(100*tmp);

						// write this read/transcript in frame?
						if(readInFrameWithTranscript != null)   
							readInFrameWithTranscript.write("\t"+"true"+"\t"+"0"+"\n");
					}
				}
				else{
					// this read is too short or too long
					if(readInFrameWithTranscript != null)
						readInFrameWithTranscript.write("\t"+"-1.0"+"\t"+"true"+"\t"+"0"+"\n");
				}
			}else{
				// this read is not in the coding sequence
				ignored_notInCDS ++;
				if(readInFrameWithTranscript != null)
					readInFrameWithTranscript.write("\t"+"0.0"+"\t"+"false"+"\t"+"0"+"\n");
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
		options.addOption(OptionBuilder.withArgName("maxReads").hasArg().withDescription("[optional] Only read the first R reads from the input bam [default: 1,000,000,000]").create("R"));
		options.addOption(OptionBuilder.withArgName("true/false").withDescription("[optional] Output read-level summary [default: false]").create("outputRead2Transcript"));
		options.addOption(OptionBuilder.withDescription("Do not print progress to stderr").create("noprog"));
		return options;
	}



	/**
	 * @param args the command line arguments
	 */
	public static void main(String[] args) throws Exception {

		/*args = new String[]{"-f","/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Footprinting/FP_original_2_Genome_Aligned.toTranscriptome.sorted.bam",
				"-a","/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation_noSelenocysteine.gtf",
				"-o","/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Footprinting/FootprintFrameAnalysis_TEST3_",
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

			// suppress printing of the progress bar
			if(cmdArgs.hasOption("noprog"))
				engine.suppressProgressBar();

			// check override for min # of footprints
			long maxReads = 1000000000;
			if(cmdArgs.hasOption("R"))
				maxReads = Integer.valueOf(cmdArgs.getOptionValue("R"));
			// Read the footprint alignments
			//engine.readAlignments(new File(cmdArgs.getOptionValue("f")), maxReads);

			boolean outputRead2TranscriptFrameMatrix = false;
			if(cmdArgs.hasOption("outputRead2Transcript"))
				outputRead2TranscriptFrameMatrix = true;

			engine.doAll(new File(cmdArgs.getOptionValue("f")), maxReads, outputPrefix, false, outputRead2TranscriptFrameMatrix);

			// Calculate frame-consistency of reads for each transcript (to be used in the EM) 
			//engine.calculateFrameConsistency(new File(cmdArgs.getOptionValue("o")+"FrameConsistencyPerTranscript.txt"));
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" FootprintFrameAnaysis", getCmdLineOptions());
			System.err.println();
		}
	}


	/**
	 * 
	 * @param readAlignments
	 * @param maxReads
	 * @param outputPrefix
	 * @param useFrameToSelectBestTranscript
	 * @param outputRead2TranscriptFrameMatrix
	 * @throws IOException
	 */
	public void doAll(File readAlignments, long maxReads, String outputPrefix, boolean useFrameToSelectBestTranscript, boolean outputRead2TranscriptFrameMatrix) throws IOException{
		// Read the footprint alignments
		readAlignments(readAlignments, maxReads);

		// Compare footprint midpoints to annotated codon positions
		computeReadMidOffset2CDS(outputPrefix, _USE_OFFSET_GUESSED, useFrameToSelectBestTranscript, false);
		writeReadLengths2offsets(outputPrefix+"_ReadLength_vs_codonOffset.txt");
		//engine.writeTranscriptFrameGuesses(outputPrefix+"_TranscriptFrames_guessed.txt");

		computeReadMidOffset2CDS(outputPrefix, _USE_OFFSET_EMPIRICAL, useFrameToSelectBestTranscript, outputRead2TranscriptFrameMatrix);
		//engine.writeTranscriptFrameGuesses(outputPrefix+"_TranscriptFrames_learned.txt");

	}

}
