package objects;

import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Comparator;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map;
import java.util.TreeMap;

import utils.IO_utils;


public class RandomBarcodeStats {

	// HashMap<insert sequence,   HashMap<barcode sequence,   number of copies of this barcode>>
	private HashMap<String, HashMap<String, Integer>> _readsAndBarcodes = new HashMap<String, HashMap<String, Integer>>();		// if sequences have BOTH and 5' and 3' barcode
	private HashMap<String, HashMap<String, Integer>> _readsAndBarcodes_5p = new HashMap<String, HashMap<String, Integer>>();	// if sequences have ONLY a 5' barcode
	private HashMap<String, HashMap<String, Integer>> _readsAndBarcodes_3p = new HashMap<String, HashMap<String, Integer>>();	// if sequences have ONLY a 3' barcode

	// HashMap<insert sequence,   number of copies of this insert>
	private HashMap<String, Integer> _readCounts = new HashMap<String, Integer>();
	private Comparator_HashInteger _bvc =  new Comparator_HashInteger(_readCounts);
	private TreeMap<String,Integer> _readCounts_sorted = new TreeMap<String,Integer>(_bvc);


	// Housekeeping variables
	private int _nBarcodeBases, _totalBarcodeBases;
	private boolean _look5p, _look3p;



	/**
	 * Instantiates the class
	 * @param barcodeLength Number of nucleotides in the random barcode on a SINGLE adapter
	 * @param is5p Is there a random barcode on the 5' adapter?
	 * @param is3p Is there a random barcode on the 3' adapter?
	 */
	public RandomBarcodeStats(int barcodeLength, boolean is5p, boolean is3p){
		_nBarcodeBases = barcodeLength;
		_look5p = is5p;
		_look3p = is3p;
		//_nBarcodeBases = barcodeLength;
		_totalBarcodeBases = barcodeLength;
		if(_look5p  &&  _look3p)
			_totalBarcodeBases = 2*barcodeLength;
	}



	/**
	 * Create running tally of the barcode - insert frequencies for computing the global barcode probabilities for the KL divergence
	 * 
	 *   Barcode sequence, insert sequences
	 */
	//private HashMap<String, ArrayList<String>> _barcodes2inserts_5p3p = new HashMap<String, ArrayList<String>>(); 
	private HashMap<String, ArrayList<String>> _barcodes2inserts_5p = new HashMap<String, ArrayList<String>>();
	private HashMap<String, ArrayList<String>> _barcodes2inserts_3p = new HashMap<String, ArrayList<String>>();

	/**
	 * Adds an insert sequence and its barcode(s) to the stats lists
	 * @param insertSequence
	 * @param barcode
	 */
	public void add(String insertSequence, String barcode){

		// Add this insert to the running counts
		if(!_readCounts.containsKey(insertSequence))
			_readCounts.put(insertSequence, 0);
		_readCounts.put(insertSequence, _readCounts.get(insertSequence) + 1);


		// If there is a barcode at each end of the insert
		if(_look5p  &&  _look3p){
			if(!_readsAndBarcodes.containsKey(insertSequence))
				_readsAndBarcodes.put(insertSequence, new HashMap<String, Integer>());

			if(!_readsAndBarcodes.get(insertSequence).containsKey(barcode))
				_readsAndBarcodes.get(insertSequence).put(barcode, 0);
			_readsAndBarcodes.get(insertSequence).put(barcode, _readsAndBarcodes.get(insertSequence).get(barcode) + 1);

			/*if(!_barcodes2inserts_5p3p.containsKey(barcode))
				_barcodes2inserts_5p3p.put(barcode, new ArrayList<String>());
			if(!_barcodes2inserts_5p3p.get(barcode).contains(insertSequence))
				_barcodes2inserts_5p3p.get(barcode).add(insertSequence);*/
		}

		// If there is a barcode at the 5' end of the insert
		if(_look5p){
			if(!_readsAndBarcodes_5p.containsKey(insertSequence))
				_readsAndBarcodes_5p.put(insertSequence, new HashMap<String, Integer>());

			String barcode_5p = barcode;
			if(_look5p  &&  _look3p)
				barcode_5p = barcode.substring(0, _nBarcodeBases);

			if(!_readsAndBarcodes_5p.get(insertSequence).containsKey(barcode_5p))
				_readsAndBarcodes_5p.get(insertSequence).put(barcode_5p, 0);
			_readsAndBarcodes_5p.get(insertSequence).put(barcode_5p, _readsAndBarcodes_5p.get(insertSequence).get(barcode_5p) + 1);

			if(!_barcodes2inserts_5p.containsKey(barcode_5p))
				_barcodes2inserts_5p.put(barcode_5p, new ArrayList<String>());
			if(!_barcodes2inserts_5p.get(barcode_5p).contains(insertSequence))
				_barcodes2inserts_5p.get(barcode_5p).add(insertSequence);
		}

		// If there is a barcode at the 3' end of the insert
		if(_look3p){
			if(!_readsAndBarcodes_3p.containsKey(insertSequence))
				_readsAndBarcodes_3p.put(insertSequence, new HashMap<String, Integer>());

			String barcode_3p = barcode;
			if(_look5p  &&  _look3p)
				barcode_3p = barcode.substring(_nBarcodeBases, _nBarcodeBases*2);

			if(!_readsAndBarcodes_3p.get(insertSequence).containsKey(barcode_3p))
				_readsAndBarcodes_3p.get(insertSequence).put(barcode_3p, 0);
			_readsAndBarcodes_3p.get(insertSequence).put(barcode_3p, _readsAndBarcodes_3p.get(insertSequence).get(barcode_3p) + 1);

			if(!_barcodes2inserts_3p.containsKey(barcode_3p))
				_barcodes2inserts_3p.put(barcode_3p, new ArrayList<String>());
			if(!_barcodes2inserts_3p.get(barcode_3p).contains(insertSequence))
				_barcodes2inserts_3p.get(barcode_3p).add(insertSequence);
		}
	}



	//private HashMap<String, Double> _globalBarcodeProbabilities_5p3p = new HashMap<String, Double>();
	private HashMap<String, Double> _globalBarcodeProbabilities_5p = new HashMap<String, Double>();
	private HashMap<String, Double> _globalBarcodeProbabilities_3p = new HashMap<String, Double>();

	/**
	 * Compute probability of each 4N barcode sequence given its global frequency 
	 */
	public void computeGlobalBarcodeProbabilities(){
		Iterator<String> it;
		if(_barcodes2inserts_5p.size() > 0){
			it = _barcodes2inserts_5p.keySet().iterator();
			int runningTotal = 0;
			while(it.hasNext())
				runningTotal += _barcodes2inserts_5p.get(it.next()).size();

			it = _barcodes2inserts_5p.keySet().iterator();
			while(it.hasNext()){
				String thisBarcode = it.next();
				_globalBarcodeProbabilities_5p.put(thisBarcode, new Double(_barcodes2inserts_5p.get(thisBarcode).size()/(runningTotal+0.0)));
				//System.out.println(thisBarcode+">"+"\t"+_barcodes2inserts_5p.get(thisBarcode).size()/(runningTotal+0.0));
			}
		}

		if(_barcodes2inserts_3p.size() > 0){
			it = _barcodes2inserts_3p.keySet().iterator();
			int runningTotal = 0;
			while(it.hasNext())
				runningTotal += _barcodes2inserts_3p.get(it.next()).size();

			it = _barcodes2inserts_3p.keySet().iterator();
			while(it.hasNext()){
				String thisBarcode = it.next();
				_globalBarcodeProbabilities_3p.put(thisBarcode, new Double(_barcodes2inserts_3p.get(thisBarcode).size()/(runningTotal+0.0)));
				//System.out.println("<"+thisBarcode+"\t"+_barcodes2inserts_3p.get(thisBarcode).size()/(runningTotal+0.0));
			}
		}
	}

	private HashMap<String, Integer> _barcodePlusLocalInsertCount_5p = new HashMap<String, Integer>();
	private HashMap<String, Integer> _barcodePlusLocalInsertCount_3p = new HashMap<String, Integer>();
	private HashMap<String, Integer> _barcode2barcodeCount = new HashMap<String, Integer>();
	private HashMap<String, Integer> _insert2insertCount = new HashMap<String, Integer>();
	
	public void computeBarcodeLigationLocalSequenceContext(String outputPath, int nBases) throws IOException{
		Iterator<String> it, it2;
		
		if(_barcodes2inserts_5p.size() > 0){
			it = _barcodes2inserts_5p.keySet().iterator();
			while(it.hasNext()){
				String thisBarcode = it.next();
				it2 = _barcodes2inserts_5p.get(thisBarcode).iterator();
				while(it2.hasNext()){
					String insertSeq = it2.next();
					String barcode2insert = thisBarcode.substring(thisBarcode.length()-nBases)+"|"+insertSeq.substring(0, nBases);
					if(!_barcodePlusLocalInsertCount_5p.containsKey(barcode2insert))
						_barcodePlusLocalInsertCount_5p.put(barcode2insert, 0);
					_barcodePlusLocalInsertCount_5p.put(barcode2insert, _barcodePlusLocalInsertCount_5p.get(barcode2insert)+1);
					
					String insert2insert = insertSeq.substring(0, nBases)+"|"+insertSeq.substring(insertSeq.length()-nBases);
					if(!_insert2insertCount.containsKey(insert2insert))
						_insert2insertCount.put(insert2insert, 0);
					_insert2insertCount.put(insert2insert, _insert2insertCount.get(insert2insert)+1);
				}
			}
			writeResult(outputPath+".global_ligationEff_5p", _barcodePlusLocalInsertCount_5p);
			writeResult(outputPath+".global_ligationEff_insert2insert", _insert2insertCount);
		}
		
		
		if(_barcodes2inserts_3p.size() > 0){
			it = _barcodes2inserts_3p.keySet().iterator();
			while(it.hasNext()){
				String thisBarcode = it.next();
				it2 = _barcodes2inserts_3p.get(thisBarcode).iterator();
				while(it2.hasNext()){
					String insertSeq = it2.next();
					String barcode2insert = insertSeq.substring(insertSeq.length()-nBases)+"|"+thisBarcode.substring(0, nBases);
					if(!_barcodePlusLocalInsertCount_3p.containsKey(barcode2insert))
						_barcodePlusLocalInsertCount_3p.put(barcode2insert, 0);
					_barcodePlusLocalInsertCount_3p.put(barcode2insert, _barcodePlusLocalInsertCount_3p.get(barcode2insert)+1);
					
					String insert2insert = insertSeq.substring(0, nBases)+"|"+insertSeq.substring(insertSeq.length()-nBases);
					if(!_insert2insertCount.containsKey(insert2insert))
						_insert2insertCount.put(insert2insert, 0);
					_insert2insertCount.put(insert2insert, _insert2insertCount.get(insert2insert)+1);
				}
			}
			writeResult(outputPath+".global_ligationEff_3p", _barcodePlusLocalInsertCount_3p);
			writeResult(outputPath+".global_ligationEff_insert2insert", _insert2insertCount);
		}
	}
	private void writeResult(String outputPath, HashMap<String, Integer> input) throws IOException{
		PrintWriter out = new PrintWriter(new FileWriter(outputPath));
		TreeMap<String,Integer> barcodePlusLocalInsertCount_5p_sorted = new TreeMap<String,Integer>(new Comparator_HashInteger(input));
		barcodePlusLocalInsertCount_5p_sorted.putAll(input);
		Iterator<String> it = barcodePlusLocalInsertCount_5p_sorted.keySet().iterator();
		while(it.hasNext()){
			String seq = it.next();
			out.write(seq+"\t"+input.get(seq)+"\n");
		}
		out.flush();
		out.close();
	}



	/**
	 * 
	 * @param outputPath
	 * @throws IOException
	 */
	public void processBarcodeStatistics(String outputPath) throws IOException{

		// compute global barcode frequencies to adjust the KL divergence computation
		IO_utils.printLineErr("Computing global barcode probabilities");
		computeGlobalBarcodeProbabilities();
		writeGlobalAdapterStats(outputPath);
		computeBarcodeLigationLocalSequenceContext(outputPath, 2);
		

		PrintWriter out = new PrintWriter(new FileWriter(outputPath));
		out.write("InsertSequence\tnInsertSequences");
		//\tKLdivergence\tchisq_testStatistic\tchisq_pValue\tPWM_A\tPWM_C\tPWM_G\tPWM_T\tmostAbundantBarcode_1\tmostAbundantBarcode_2\tmostAbundantBarcode_3\n");

		String pwmHeader = "\tPWM_A\tPWM_C\tPWM_G\tPWM_T\n";
		boolean printPWM = !_look5p&&_look3p; // only print PWM for 3' or 5' when there aren't both barcodes present

		if(_look5p){
			//out.write("\tnUniqueBarcodes_5p\tKLdivergence_5p\tchisq_testStatistic_5p\tchisq_pValue_5p\ttopBarcode_5p\tsecondBarcode_5p\tthirdBarcode_5p");
			out.write("\tnUniqueBarcodes_5p\tKLdivergence_5p\ttopBarcode_5p\tsecondBarcode_5p\tthirdBarcode_5p");
			if(printPWM)
				out.write(pwmHeader);
		}if(_look3p){
			//out.write("\tnUniqueBarcodes_3p\tKLdivergence_3p\tchisq_testStatistic_3p\tchisq_pValue_3p\ttopBarcode_3p\tsecondBarcode_3p\tthirdBarcode_3p");
			out.write("\tnUniqueBarcodes_3p\tKLdivergence_3p\ttopBarcode_3p\tsecondBarcode_3p\tthirdBarcode_3p");
			if(printPWM)
				out.write(pwmHeader);
		}if(_look5p  &&  _look3p){
			//out.write("\tnUniqueBarcodes_5p3p\tKLdivergence_5p3p\tchisq_testStatistic_5p3p\tchisq_pValue_5p3p\ttopBarcode_5p3p\tsecondBarcode_5p3p\tthirdBarcode_5p3p"+pwmHeader);
			out.write("\tnUniqueBarcodes_5p3p\tKLdivergence_5p3p\ttopBarcode_5p3p\tsecondBarcode_5p3p\tthirdBarcode_5p3p"+pwmHeader);
		}

		// sort inserts based on decreasing frequency
		_readCounts_sorted.putAll(_readCounts);
		Iterator<String> it_reads = _readCounts_sorted.keySet().iterator();

		String thisReadSeq;
		while(it_reads.hasNext()){
			thisReadSeq = it_reads.next();
			int totalCount = _readCounts.get(thisReadSeq);

			// write this insert sequence and the number of times this insert was sequenced
			out.write(thisReadSeq+"\t"+totalCount);

			// write the number of times each barcode was observed, broken down by barcode position
			if(_look5p){
				out.write("\t" + _readsAndBarcodes_5p.get(thisReadSeq).size());
				out.write("\t" + processBarcodeStatsForInsert(thisReadSeq, _readsAndBarcodes_5p.get(thisReadSeq), totalCount, _nBarcodeBases, BARCODE_POSITION_5P, printPWM));
				//out.write("\t" + processBarcodeStatsForInsert(thisReadSeq, _readsAndBarcodes_5p.get(thisReadSeq), totalCount, 4, printPWM));
			}if(_look3p){
				out.write("\t" + _readsAndBarcodes_3p.get(thisReadSeq).size());
				out.write("\t" + processBarcodeStatsForInsert(thisReadSeq, _readsAndBarcodes_3p.get(thisReadSeq), totalCount, _nBarcodeBases, BARCODE_POSITION_3P, printPWM));
			}
			if(_look5p  &&  _look3p){
				out.write("\t" + _readsAndBarcodes.get(thisReadSeq).size());
				out.write("\t" + processBarcodeStatsForInsert(thisReadSeq, _readsAndBarcodes.get(thisReadSeq), totalCount, 2*_nBarcodeBases, BARCODE_POSITION_5P3P, true));
			}
			out.write("\n");

			//insertCount++;
		}
		out.flush();
		out.close();
		
		writeResult(outputPath+".global_ligationEff_barcode2barcode", _barcode2barcodeCount);
	}



	public static final int BARCODE_POSITION_5P3P = 0;
	public static final int BARCODE_POSITION_5P = 1;
	public static final int BARCODE_POSITION_3P = 2;


	/**
	 * 
	 * @param insertSequence
	 * @param theseBarcodes
	 * @param totalCount
	 * @param barcodeBases
	 * @param doPWM
	 * @return
	 */
	private String processBarcodeStatsForInsert(String insertSequence, HashMap<String, Integer> theseBarcodes, int totalCount, int barcodeBases, int barcodePosition, boolean doPWM){
		//private String processBarcodeStatsForInsert(String insertSequence, HashMap<String, Integer> theseBarcodes, int totalCount, int barcodeBases, boolean doPWM){
		String outputString = "";
		String barcodeOutputString = "";

		TreeMap<String,Integer> barcodeCounts_sorted = new TreeMap<String,Integer>(new Comparator_HashInteger(theseBarcodes));
		barcodeCounts_sorted.putAll(theseBarcodes);
		Iterator<String> it_barcodes = barcodeCounts_sorted.keySet().iterator();

		Double entropy = 0.0;
		double[][] pwm = makeEmptyPWM();

		//KolmogorovSmirnovTest ks = new KolmogorovSmirnovTest();
		//ChiSquareTest chisq = new ChiSquareTest();
		//long[] observedCounts = new long[(int)Math.pow(4, barcodeBases)];
		double[] expectedCounts = new double[(int)Math.pow(4, barcodeBases)];
		double expectedFrac = 1.0/Math.pow(4, barcodeBases);
		for(int i=0;i<expectedCounts.length;i++){ expectedCounts[i] = expectedFrac; }

		int count = 0;
		while(it_barcodes.hasNext()){
			String thisBarcode = it_barcodes.next();
			int thisBarcodeCount = theseBarcodes.get(thisBarcode);

			// add 5p and/or 3p barcode counts to the global stats
			//updateGlobalAdapterStats(thisBarcode, thisBarcodeCount, barcodeBases, barcodePosition, insertSequence);
			//updateGlobalAdapterStats(thisBarcode, 1, barcodeBases, barcodePosition, insertSequence);
			
			// Add the end of the 5' and the beginning of the 3' barcode to the global count to compare against barcode-insert ligation stats
			if(barcodePosition == BARCODE_POSITION_5P3P){
				String barcodeBits = thisBarcode.substring((int)Math.round(barcodeBases/2.0)-2, (int)Math.round(barcodeBases/2.0)+2);
				barcodeBits = barcodeBits.substring(0,2)+"|"+barcodeBits.substring(2);
				//System.out.println(thisBarcode+"\t"+barcodeBits);
				if(!_barcode2barcodeCount.containsKey(barcodeBits))
					_barcode2barcodeCount.put(barcodeBits, 0);
				_barcode2barcodeCount.put(barcodeBits, _barcode2barcodeCount.get(barcodeBits)+1);
				
			}

			// add this barcode to the PWM
			pwm = addBarcodeToPWM(pwm, thisBarcode, theseBarcodes.get(thisBarcode), totalCount);

			// calculate the kulback-liebler (normalised entropy)
			double nPoss = totalCount; 
			if(nPoss > Math.pow(4, barcodeBases)){ nPoss = Math.pow(4, barcodeBases); }
			double barcodeProb = 1.0/(nPoss+0.0);
			if(barcodePosition == BARCODE_POSITION_5P)
				barcodeProb = this._globalBarcodeProbabilities_5p.get(thisBarcode);
			else if(barcodePosition == BARCODE_POSITION_3P)
				barcodeProb = this._globalBarcodeProbabilities_3p.get(thisBarcode);
			entropy -= calculateKL(barcodeProb, thisBarcodeCount, barcodeBases, theseBarcodes.size());

			// add the observed count for this barcode
			//observedCounts[count] = thisBarcodeCount;

			// format barcode and count for printing
			if(count < 3)
				barcodeOutputString = barcodeOutputString.concat("\t"+thisBarcode+" "+thisBarcodeCount);
			count ++;
		}

		// if we don't have three barcodes, pad the space
		for(int i=count;i<3;i++)
			barcodeOutputString = barcodeOutputString.concat("\t");


		//System.err.print("\t"+pwmToString(pwm));
		//System.err.println(barcodeOutputString);

		outputString = outputString.concat(""+Math.round(entropy*1000)/1000.0);
		//out.printf("\t%e\t%e",chisq.chiSquare(expectedCounts, observedCounts), chisq.chiSquareTest(expectedCounts, observedCounts));
		//outputString = outputString.concat("\t" + chisq.chiSquare(expectedCounts, observedCounts) +"\t"+ chisq.chiSquareTest(expectedCounts, observedCounts));
		outputString = outputString.concat(barcodeOutputString);
		if(doPWM)
			outputString = outputString.concat(pwmToString(pwm));

		return outputString;
	}




	/**
	 * 
	 * @param outputPath
	 * @throws IOException
	 */
	public void writeGlobalAdapterStats(String outputPath) throws IOException{
		if(_look5p){
			PrintWriter out = new PrintWriter(new FileWriter(outputPath+".global_5p"));
			TreeMap<String,Double> globalBarcodeProbabilities_5p_sorted = new TreeMap<String,Double>(new Comparator_HashDouble(_globalBarcodeProbabilities_5p));
			globalBarcodeProbabilities_5p_sorted.putAll(_globalBarcodeProbabilities_5p);
			Iterator<String> it = globalBarcodeProbabilities_5p_sorted.keySet().iterator();
			while(it.hasNext()){
				String barcode = it.next();
				out.write(barcode+"\t"+_barcodes2inserts_5p.get(barcode).size()+"\t"+_globalBarcodeProbabilities_5p.get(barcode)+"\n");
			}
			out.flush();
			out.close();
		}
		if(_look3p){
			PrintWriter out = new PrintWriter(new FileWriter(outputPath+".global_3p"));
			TreeMap<String,Double> globalBarcodeProbabilities_3p_sorted = new TreeMap<String,Double>(new Comparator_HashDouble(_globalBarcodeProbabilities_3p));
			globalBarcodeProbabilities_3p_sorted.putAll(_globalBarcodeProbabilities_3p);
			Iterator<String> it = globalBarcodeProbabilities_3p_sorted.keySet().iterator();
			while(it.hasNext()){
				String barcode = it.next();
				out.write(barcode+"\t"+_barcodes2inserts_3p.get(barcode).size()+"\t"+_globalBarcodeProbabilities_3p.get(barcode)+"\n");
			}
			out.flush();
			out.close();
		}

	}



	/*
	//private int _adapterIndex = 0;
	private void addToGlobalAdapterCounts(String tmpBarcodeString, String barcodePosition, int count){
		if(!_globalAdapterCounts.containsKey(barcodePosition))
			_globalAdapterCounts.put(barcodePosition, new HashMap<String, Integer>());

		if(!_globalAdapterCounts.get(barcodePosition).containsKey(tmpBarcodeString))
			_globalAdapterCounts.get(barcodePosition).put(tmpBarcodeString, 0);

		_globalAdapterCounts.get(barcodePosition).put(tmpBarcodeString, _globalAdapterCounts.get(barcodePosition).get(tmpBarcodeString)+count);
	}

	private void addToAdaptersVsInserts(String tmpBarcodeString, int count, String thisInsert){
		if(!_adaptersVsInserts.containsKey(tmpBarcodeString))
			_adaptersVsInserts.put(tmpBarcodeString, new HashMap<String,Integer>());
		_adaptersVsInserts.get(tmpBarcodeString).put(thisInsert, count);
	}



	private HashMap<String, HashMap<String,Integer>> _globalAdapterCounts = new HashMap<String, HashMap<String,Integer>>();
	private HashMap<String,Integer> _globalAdapterCounts_5p = new HashMap<String,Integer>();
	private HashMap<String,Integer> _globalAdapterCounts_3p = new HashMap<String,Integer>();
	private HashMap<String,HashMap<String,Integer>> _adaptersVsInserts = new HashMap<String,HashMap<String,Integer>>();




	 */

	/**
	 * 
	 * @param thisBarcode
	 * @param thisBarcodeCount
	 * @param thisInsert
	 */
	/*
	private void updateGlobalAdapterStats(String thisBarcode, int thisBarcodeCount, int barcodeBases, int barcodePosition, String thisInsert){
		if(barcodePosition == BARCODE_POSITION_5P3P){
			System.out.println(thisBarcode.substring(0, barcodeBases/2)+"<>"+thisBarcode.substring(barcodeBases/2)+"\t"+thisBarcodeCount);

			addToGlobalAdapterCounts(thisBarcode.substring(0, barcodeBases/2)+"<>"+thisBarcode.substring(barcodeBases/2), "5p3p", thisBarcodeCount);
			//addToGlobalAdapterCounts(thisBarcode.substring(0, barcodeBases/2)+"<", "5p", thisBarcodeCount);
			//addToGlobalAdapterCounts(">"+thisBarcode.substring(barcodeBases/2), "3p", thisBarcodeCount);
			//addToAdaptersVsInserts(thisBarcode.substring(0, barcodeBases)+"<", thisBarcodeCount, thisInsert);
			//addToAdaptersVsInserts(">"+thisBarcode.substring(barcodeBases, thisBarcode.length()-barcodeBases), thisBarcodeCount, thisInsert);
		}else if(barcodePosition == BARCODE_POSITION_5P){
			System.out.println(thisBarcode+"<");
			addToGlobalAdapterCounts(thisBarcode+"<", "5p", thisBarcodeCount);
			//addToAdaptersVsInserts(thisBarcode.substring(0, barcodeBases)+"<", thisBarcodeCount, thisInsert);
		}else if(barcodePosition == BARCODE_POSITION_3P){
			System.out.println(">"+thisBarcode);
			addToGlobalAdapterCounts(">"+thisBarcode, "3p", thisBarcodeCount);
			//addToAdaptersVsInserts(">"+thisBarcode.substring(0, _nBarcodeBases), thisBarcodeCount, thisInsert);
		}
	}
	 */






	/**
	 * Calculates the Kulbach-Liebler divergence for the given barcode counts
	 * @param numberOfObservationsOfALLBarcodeSequences
	 * @param numberOfObservationsOfThisBarcodeSequence
	 * @param nBarcodeBases
	 * @param nObservedBarcodes
	 * @return
	 */
	public static double calculateKL(double barcodeProbability, int numberOfObservationsOfThisBarcodeSequence, int nBarcodeBases, int nObservedBarcodes){

		double frequency = (double) numberOfObservationsOfThisBarcodeSequence * barcodeProbability;
		return frequency * (Math.log(frequency*nObservedBarcodes) / Math.log(2));
	}


	/**
	 * Functions for computing the PWM of bases in the random barcode
	 * @return
	 */
	private double[][] makeEmptyPWM(){
		double[][] pwm = new double[4][_totalBarcodeBases];
		for(int i=0;i<pwm.length;i++)
			for(int j=0;j<pwm[i].length;j++)
				pwm[i][j]= 0;
		return pwm;
	}
	private double[][] addBarcodeToPWM(double[][] pwm, String barcode, int nBarcodes, int totalCount){
		for(int x=0;x<barcode.length();x++){
			char thisChar = barcode.charAt(x);
			if(thisChar == 'A')
				pwm[0][x] += (nBarcodes+0.0)/(totalCount+0.0);
			else if(thisChar == 'C')
				pwm[1][x] += (nBarcodes+0.0)/(totalCount+0.0);
			else if(thisChar == 'G')
				pwm[2][x] += (nBarcodes+0.0)/(totalCount+0.0);
			else if(thisChar == 'T')
				pwm[3][x] += (nBarcodes+0.0)/(totalCount+0.0);
		}
		return pwm;
	}
	private String[] bases = new String[]{"A","C","G","T"};
	private String pwmToString(double[][] pwm){
		String result = "";
		for(int i=0;i<pwm.length;i++){
			result = result.concat("\t"+bases[i]+":{");
			for(int j=0;j<pwm[i].length-1;j++){
				result = result.concat(Math.round(pwm[i][j]*1000)/1000.0+" ");
			}
			result = result.concat(Math.round(pwm[i][pwm[i].length-1]*1000)/1000.0+"}");
		}
		return result ;
	}

}



class Comparator_HashInteger implements Comparator<String> {

	Map<String, Integer> base;
	public Comparator_HashInteger(Map<String, Integer> base) {
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

class Comparator_HashDouble implements Comparator<String> {

	Map<String, Double> base;
	public Comparator_HashDouble(Map<String, Double> base) {
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