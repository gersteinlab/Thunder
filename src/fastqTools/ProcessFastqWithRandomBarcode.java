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

public class ProcessFastqWithRandomBarcode {

	private SequenceReader _reader;
	private int _nBarcodeBases, _totalBarcodeBases;
	private boolean _look5p, _look3p, _calcBarcodeStats;
	
	public ProcessFastqWithRandomBarcode(String inputPath, int nBarcodeBases, boolean look5p, boolean look3p, boolean calcBarcodeStats) throws IOException{
		_reader = new SequenceReader(inputPath);
		_nBarcodeBases = nBarcodeBases;
		_look5p = look5p;
		_look3p = look3p;
		_calcBarcodeStats = calcBarcodeStats;
		
		_totalBarcodeBases = _nBarcodeBases; 
		if(_look5p && _look3p)
			_totalBarcodeBases = 2*_nBarcodeBases;
	}

	public void processFastqFile() throws IOException{
		SequenceRecord tmpIn, tmpOut;
		while((tmpIn= _reader.readNextRecord()) != null){
			tmpOut = processRead(tmpIn);
			if(tmpOut.getSequence().length() > 0)
				System.out.println(tmpOut.toString());
		}
	}

	public void processFastqFile(BufferedWriter bw) throws IOException{
		SequenceRecord tmpIn, tmpOut;
		while((tmpIn= _reader.readNextRecord()) != null){
			tmpOut = processRead(tmpIn);
			if(tmpOut.getSequence().length() > 0)
				bw.write(tmpOut.toString()+"\n");
		}
	}

	private HashMap<String, HashMap<String, Integer>> _readsAndBarcodes = new HashMap<String, HashMap<String, Integer>>(); 
	private HashMap<String, Integer> _readCounts = new HashMap<String, Integer>();
	private ValueComparator _bvc =  new ValueComparator(_readCounts);
    private TreeMap<String,Integer> _readCounts_sorted = new TreeMap<String,Integer>(_bvc);
    
	public SequenceRecord processRead(SequenceRecord in){
		String seq = in.getSequence();
		String qual = in.getQuality();
		//System.err.println(seq);
		//System.err.println(qual);

		String barcode="", trimmedReadSeq="", trimmedReadQual="";
		if(_look5p  &&  _look3p  &&  (seq.length()-1) > _nBarcodeBases*2){
			String barcode5p = seq.substring(0, _nBarcodeBases);
			String barcode3p = seq.substring(seq.length()-_nBarcodeBases, seq.length());
			barcode = barcode5p.concat(barcode3p); 
			trimmedReadSeq = seq.substring(_nBarcodeBases, seq.length()-_nBarcodeBases);
			trimmedReadQual = qual.substring(_nBarcodeBases, qual.length()-_nBarcodeBases);
		}else if(_look5p  &&  (seq.length()-1) > _nBarcodeBases){
			barcode = seq.substring(0, _nBarcodeBases);
			trimmedReadSeq = seq.substring(_nBarcodeBases, seq.length());
			trimmedReadQual = qual.substring(_nBarcodeBases, qual.length());
		}else if(_look3p  &&  (seq.length()-1) > _nBarcodeBases){
			barcode = seq.substring(seq.length()-_nBarcodeBases, seq.length());
			trimmedReadSeq = seq.substring(0, seq.length()-_nBarcodeBases);
			trimmedReadQual = qual.substring(0, qual.length()-_nBarcodeBases);
		}

		if(_calcBarcodeStats){
			if(trimmedReadSeq.length() > 0){
				if(!_readsAndBarcodes.containsKey(trimmedReadSeq)){
					_readsAndBarcodes.put(trimmedReadSeq, new HashMap<String, Integer>());
					_readCounts.put(trimmedReadSeq, 0);
				}
				if(!_readsAndBarcodes.get(trimmedReadSeq).containsKey(barcode)){
					_readsAndBarcodes.get(trimmedReadSeq).put(barcode, 0);
				}
				_readsAndBarcodes.get(trimmedReadSeq).put(barcode, _readsAndBarcodes.get(trimmedReadSeq).get(barcode) + 1);
				_readCounts.put(trimmedReadSeq, _readCounts.get(trimmedReadSeq) + 1);
			}
		}

		SequenceRecord toReturn = new SequenceRecord(in.getSequenceID());
		toReturn.addSequenceString(trimmedReadSeq);
		toReturn.addQualityString(trimmedReadQual);

		return toReturn;
	}

	
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
			result = result.concat(bases[i]+":{");
			for(int j=0;j<pwm[i].length-1;j++){
				result = result.concat(Math.round(pwm[i][j]*1000)/1000.0+" ");
			}
			result = result.concat(Math.round(pwm[i][pwm[i].length-1]*1000)/1000.0+"}\t");
		}
		return result ;
	}

	public void processBarcodeStatistics(){
		//KolmogorovSmirnovTest ks = new KolmogorovSmirnovTest();
		
		_readCounts_sorted.putAll(_readCounts);
		Iterator<String> it_reads = _readCounts_sorted.keySet().iterator();
		
		String thisReadSeq;
		while(it_reads.hasNext()){
			thisReadSeq = it_reads.next();
			int totalCount = _readCounts.get(thisReadSeq);

			double[][] pwm = makeEmptyPWM();
			
			System.err.print(thisReadSeq+"\t"+totalCount+"\t"+_readsAndBarcodes.get(thisReadSeq).size());
			String barcodeOutputString = "";
			
			TreeMap<String,Integer> barcodeCounts_sorted = new TreeMap<String,Integer>(new ValueComparator(_readsAndBarcodes.get(thisReadSeq)));
			barcodeCounts_sorted.putAll(_readsAndBarcodes.get(thisReadSeq));
			Iterator<String> it_barcodes = barcodeCounts_sorted.keySet().iterator();
			Double entropy = 0.0;
			int count = 0;
			while(it_barcodes.hasNext()){
				String thisBarcode = it_barcodes.next();

				// add this barcode to the PWM
				pwm = addBarcodeToPWM(pwm, thisBarcode, _readsAndBarcodes.get(thisReadSeq).get(thisBarcode), totalCount);
				
				// calculate the kulback-liebler (normalised entropy)
				double nPoss = totalCount; 
				if(nPoss > Math.pow(4, 8))
					nPoss = Math.pow(4, 8);
				double frequency = (double) _readsAndBarcodes.get(thisReadSeq).get(thisBarcode) / totalCount;
				entropy -= frequency * (Math.log(frequency*_readsAndBarcodes.get(thisReadSeq).size()) / Math.log(2));
				
				// format barcode and count for printing
				if(count < 3)
					barcodeOutputString = barcodeOutputString.concat("\t"+thisBarcode+" "+_readsAndBarcodes.get(thisReadSeq).get(thisBarcode));
				count ++;
			}

			if(totalCount > 0){
				//System.err.print("\tpVal_KS="+ks.kolmogorovSmirnovTest(barcodeCounts_db, uniformCounts));
				System.err.print("\t"+Math.round(entropy*1000)/1000.0);
			}//else
				//System.err.print("\tentropy=NA");
			
			System.err.print("\t"+pwmToString(pwm));
			System.err.println(barcodeOutputString);
		}
	}


	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("output sequences shorter than the maximum length to this file [if not specified, sequences are printed to stdout]").create(Thunder.OPT_PATH_OUTPUT));
		options.addOption(OptionBuilder.withArgName("randomBarcodeBaseCount").hasArg().withDescription("Number of degenerate bases in the random barcode [default: 4]").create("n"));
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("Maximum insert size (read length including 3' adapter)").create("max"));
		options.addOption(OptionBuilder.withArgName("barcode3p").withDescription("Look for barcode at 3' end of the read").create("3p"));
		options.addOption(OptionBuilder.withArgName("barcode5p").withDescription("Look for barcode at 5' end of the read").create("5p"));
		options.addOption(OptionBuilder.withArgName("boolean").withDescription("Calculate statistics on the random barcodes").create("stats"));
		return options;
	}


	/**
	 * @param args
	 */
	public static void main(String[] args) throws Exception{

		//args = new String[]{"ProcessFastqWithRandomBarcode", "-n","4", "--3p", "--5p", "-"+Thunder.OPT_PATH_OUTPUT,"/Users/robk/Downloads/randomBarcodeTest.out.fq", "/Users/robk/Downloads/randomBarcodeTest.fq"};
		//args = new String[]{"ProcessFastqWithRandomBarcode", "-n","4", "--3p", "--5p", "/Users/robk/Downloads/randomBarcodeTest.fq"};
		//args = new String[]{"ProcessFastqWithRandomBarcode", "-n", "4", "--5p", "--stats", "/Users/robk/Downloads/randomBarcodeTest.fq"};
		//args = new String[]{"ProcessFastqWithRandomBarcode", "-n", "4", "--3p", "--stats", "/Users/robk/Downloads/randomBarcodeTest.fq"};
		//args = new String[]{"ProcessFastqWithRandomBarcode", "-n", "4", "--stats", "/Users/robk/Downloads/randomBarcodeTest.fq"};

		//args = new String[]{"ProcessFastqWithRandomBarcode", "-n", "4", "--3p", "--5p", "--stats", "/Users/robk/Downloads/randomBarcodeTest.fq"};

		//"/Users/robk/Downloads/miR-Pool-5pmoles-TruSeqRecJ-4N-2_S8_L001_R1_001_TEST_RandomBarcode_miR-Pool-5pmoles-TruSeqRecJ-4N-2_S8_L001_R1_001.clipped.fastq"
		//args = new String[]{"ProcessFastqWithRandomBarcode", "-n", "4", "--3p", "--5p", "--stats", "/Users/robk/Downloads/test.fq"};

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		//@SuppressWarnings("unchecked")
		//Iterator<String> it = cmdArgs.getArgList().iterator();

		if(cmdArgs.getArgList().size() == 2){
			//System.out.println((String) cmdArgs.getArgList().get(1));

			int nBarcodeBases = 4;
			if(cmdArgs.hasOption("n"))
				nBarcodeBases = Integer.valueOf(cmdArgs.getOptionValue("n")).intValue();
			//System.out.println("nBarcodeBases="+nBarcodeBases);

			ProcessFastqWithRandomBarcode engine = new ProcessFastqWithRandomBarcode((String) cmdArgs.getArgList().get(1), nBarcodeBases, cmdArgs.hasOption("5p"), cmdArgs.hasOption("3p"), cmdArgs.hasOption("stats"));

			if(cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT)){
				BufferedWriter bw = new BufferedWriter(new FileWriter(cmdArgs.getOptionValue(Thunder.OPT_PATH_OUTPUT)));
				engine.processFastqFile(bw);
				bw.flush();
				bw.close();
			}else{
				engine.processFastqFile();
			}

			if(cmdArgs.hasOption("stats"))
				engine.processBarcodeStatistics();


		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" ProcessFastqWithRandomBarcode [options] <sequenceFile>", "", getCmdLineOptions(), "");
			System.out.println();
		}
	}
}	




class ValueComparator implements Comparator<String> {

    Map<String, Integer> base;
    public ValueComparator(Map<String, Integer> base) {
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
