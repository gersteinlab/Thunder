package fastqTools;

import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;

import main.Thunder;
import objects.RandomBarcodeStats;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import utils.IO_utils;

public class ProcessFastqWithRandomBarcode {

	private SequenceReader _reader;
	private int _nBarcodeBases;
	private boolean _look5p, _look3p, _calcBarcodeStats;
	private String _barcodesLocatedAt;
	private RandomBarcodeStats _stats;


	/**
	 * 
	 * @param inputPath
	 * @param nBarcodeBases
	 * @param look5p
	 * @param look3p
	 * @param calcBarcodeStats
	 * @throws IOException
	 */
	public ProcessFastqWithRandomBarcode(String inputPath, int nBarcodeBases, boolean look5p, boolean look3p, boolean calcBarcodeStats) throws IOException{
		_reader = new SequenceReader(inputPath);
		_nBarcodeBases = nBarcodeBases;
		_look5p = look5p;
		_look3p = look3p;
		_calcBarcodeStats = calcBarcodeStats;

		if(_look5p && _look3p)
			_barcodesLocatedAt = "5' and 3'";
		else if(_look5p)
			_barcodesLocatedAt = "5'";
		else if(_look3p)
			_barcodesLocatedAt = "3'";

		_stats = new RandomBarcodeStats(nBarcodeBases, _look5p, _look3p);
	}



	/**
	 * 
	 * @throws IOException
	 */
	public void processFastqFile() throws IOException{
		IO_utils.printLineErr("Removing "+_nBarcodeBases+"N random barcodes "+_barcodesLocatedAt+" to insert");
		SequenceRecord tmpIn, tmpOut;
		while((tmpIn= _reader.readNextRecord()) != null){
			tmpOut = processRead(tmpIn);
			if(tmpOut != null)
				System.out.println(tmpOut.toString());
		}
		System.err.println();
		IO_utils.printLineErr("\tInput "+_nReads_in+" reads");
		IO_utils.printLineErr("\tRemoved "+_nReads_suppressed_dimer+" barcode/adapter dimers");
		IO_utils.printLineErr("\tRemoved "+_nReads_suppressed_tooShort+" inserts that were shorter than "+_minInsertLength+" nt");
		IO_utils.printLineErr("\tOutput "+_nReads_out+" reads"); 
	}


	/**
	 * 
	 * @param bw
	 * @throws IOException
	 */
	public void processFastqFile(BufferedWriter bw) throws IOException{
		IO_utils.printLineErr("Removing "+_nBarcodeBases+"N random barcodes "+_barcodesLocatedAt+" to insert");
		SequenceRecord tmpIn, tmpOut;
		while((tmpIn= _reader.readNextRecord()) != null){
			tmpOut = processRead(tmpIn);
			if(tmpOut != null)
				bw.write(tmpOut.toString()+"\n");
		}
		IO_utils.printLineErr("\n\tInput "+_nReads_in+" reads");
		IO_utils.printLineErr("\tRemoved "+_nReads_suppressed_dimer+" barcode/adapter dimers");
		IO_utils.printLineErr("\tRemoved "+_nReads_suppressed_tooShort+" inserts that were shorter than "+_minInsertLength+" nt");
		IO_utils.printLineErr("\tOutput "+_nReads_out+" reads");
	}



	private int _nReads_in=0, _nReads_out=0, _nReads_suppressed_dimer=0, _nReads_suppressed_tooShort=0;
	private double _tmpSecs = 0.0;
	/**
	 * 
	 * @param in
	 * @return
	 */
	public SequenceRecord processRead(SequenceRecord in){
		String seq = in.getSequence();
		String qual = in.getQuality();
		_nReads_in ++;

		int n = 100000;
		if(_nReads_in % n == 0){
			System.err.print("\r"+IO_utils.getTime()+" Processed "+_nReads_in+" reads ("+Math.round(n/((System.currentTimeMillis()/1000.0)-_tmpSecs))+" reads/sec)     ");
			_tmpSecs = System.currentTimeMillis()/1000.0;
		}

		String barcode="", trimmedReadSeq="", trimmedReadQual="", newFastqHeader="";
		if(_look5p  &&  _look3p  &&  (seq.length()-1) > _nBarcodeBases*2){
			String barcode5p = seq.substring(0, _nBarcodeBases);
			String barcode3p = seq.substring(seq.length()-_nBarcodeBases, seq.length());
			barcode = barcode5p.concat(barcode3p); 
			trimmedReadSeq = seq.substring(_nBarcodeBases, seq.length()-_nBarcodeBases);
			trimmedReadQual = qual.substring(_nBarcodeBases, qual.length()-_nBarcodeBases);
			newFastqHeader = in.getSequenceID().concat(" {5p"+barcode5p+"}{3p"+barcode3p+"}");
		}
		else if(_look5p  &&  (seq.length()-1) > _nBarcodeBases){
			barcode = seq.substring(0, _nBarcodeBases);
			trimmedReadSeq = seq.substring(_nBarcodeBases, seq.length());
			trimmedReadQual = qual.substring(_nBarcodeBases, qual.length());
			newFastqHeader = in.getSequenceID().concat(" {5p"+barcode+"}");
		}
		else if(_look3p  &&  (seq.length()-1) > _nBarcodeBases){
			barcode = seq.substring(seq.length()-_nBarcodeBases, seq.length());
			trimmedReadSeq = seq.substring(0, seq.length()-_nBarcodeBases);
			trimmedReadQual = qual.substring(0, qual.length()-_nBarcodeBases);
			newFastqHeader = in.getSequenceID().concat(" {3p"+barcode+"}");
		}

		if(trimmedReadSeq.length() > 0){
			if(trimmedReadSeq.length() >= _minInsertLength){
				_nReads_out ++;

				if(_calcBarcodeStats)
					_stats.add(trimmedReadSeq, barcode);

				SequenceRecord toReturn = new SequenceRecord(newFastqHeader.replaceAll(" ", "_"));
				toReturn.addSequenceString(trimmedReadSeq);
				toReturn.addQualityString(trimmedReadQual);
				return toReturn;
			}else{
				_nReads_suppressed_tooShort ++;
				return null;
			}
		}else{
			_nReads_suppressed_dimer ++;
			return null;
		}
	}

	private int _minInsertLength = 16;
	public void setMinInsertLength(int minLength){ _minInsertLength = minLength; } 

	public void processBarcodeStatistics(String outputPath) throws IOException{
		IO_utils.printLineErr("Processing barcode statistics");
		_stats.processBarcodeStatistics(outputPath);
		//_stats.writeGlobalAdapterStats(outputPath+".global");
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
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("Minimum insert size (nt) following barcode removal [default: 16]").create("min"));
		options.addOption(OptionBuilder.withArgName("barcode3p").withDescription("Look for barcode at 3' end of the read").create("3p"));
		options.addOption(OptionBuilder.withArgName("barcode5p").withDescription("Look for barcode at 5' end of the read").create("5p"));
		options.addOption(OptionBuilder.withArgName("path").hasArg().withDescription("Calculate statistics on the random barcodes and write to the given file").create("stats"));
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
		//args = new String[]{"ProcessFastqWithRandomBarcode", "-n","4", "--3p", "--5p", "--stats","/Users/robk/Downloads/randomBarcodeTest.fq.stats", "/Users/robk/Downloads/randomBarcodeTest.fq"};

		//args = new String[]{"ProcessFastqWithRandomBarcode", "-n", "4", "--3p", "--5p", "--stats", "/Users/robk/Downloads/randomBarcodeTest.fq"};

		//"/Users/robk/Downloads/miR-Pool-5pmoles-TruSeqRecJ-4N-2_S8_L001_R1_001_TEST_RandomBarcode_miR-Pool-5pmoles-TruSeqRecJ-4N-2_S8_L001_R1_001.clipped.fastq"

		/*args = new String[]{"ProcessFastqWithRandomBarcode", 
				"-n", "4", 
				"--3p", 
				"--5p", 
				"--stats", "/Users/robk/Downloads/test.fq.stats.txt", 
				"-o","/Users/robk/Downloads/miR-Pool.clipped.norand.fastq", 
				"/Users/robk/Downloads/miR-Pool-5pmoles-TruSeqRecJ-4N-2_S8_L001_R1_001_TEST_RandomBarcode_miR-Pool-5pmoles-TruSeqRecJ-4N-2_S8_L001_R1_001.clipped.fastq"};
		 */

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
			
			// if the user has specified a different minimum insert size
			if(cmdArgs.hasOption("min"))
				engine.setMinInsertLength(Integer.valueOf(cmdArgs.getOptionValue("min")).intValue());
			

			if(cmdArgs.hasOption(Thunder.OPT_PATH_OUTPUT)){
				BufferedWriter bw = new BufferedWriter(new FileWriter(cmdArgs.getOptionValue(Thunder.OPT_PATH_OUTPUT)));
				engine.processFastqFile(bw);
				bw.flush();
				bw.close();
			}else{
				engine.processFastqFile();
			}

			if(cmdArgs.hasOption("stats")){
				engine.processBarcodeStatistics(cmdArgs.getOptionValue("stats"));
			}

			IO_utils.printLineErr("All done");
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" ProcessFastqWithRandomBarcode [options] <sequenceFile (or '-' for stdin)>", "", getCmdLineOptions(), "");
			System.out.println();
		}
	}
}	



