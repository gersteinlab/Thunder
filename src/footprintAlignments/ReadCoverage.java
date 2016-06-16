package footprintAlignments;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import objects.SAMRecordReduced;
import objects.Transcript;
import samTools.SAMReader;
import statsTools.Stats;
import utils.IO_utils;
import utils.Object_utils;
import annotation.ReadGTF;
import annotation.TranscriptAnnotation;

public class ReadCoverage {

	private HashMap<String, TranscriptCoverage> _referenceID2RefObj = new HashMap<String, TranscriptCoverage>();
	private TranscriptAnnotation _annotation;

	public ReadCoverage(String annotationPath, int totalBins) throws Exception{
		_annotation = ReadGTF.readGTF(annotationPath);
		//_annotation.getTranscript("ENSMUST00000202381.1").printDebugInfo();
		//_annotation.getTranscript("ENSMUST00000058793.10").printDebugInfo();
		getUTRandCDSbins(totalBins);
	}


	private double _averageLength_CDS = 0.0, _averageLength_UTR5 = 0.0, _averageLength_UTR3 = 0.0;
	private int _nBinsCDS, _nBins_UTR5, _nBins_UTR3;

	/**
	 * 
	 * @param totalBins
	 */
	private void getUTRandCDSbins(int totalBins){
		int count = 0;
		ArrayList<Integer> tmp_lengthsCDS = new ArrayList<Integer>();
		ArrayList<Integer> tmp_lengthsUTR5 = new ArrayList<Integer>();
		ArrayList<Integer> tmp_lengthsUTR3 = new ArrayList<Integer>();

		for(String transcriptID: _annotation.getTranscripts().keySet()){
			Transcript tmp = _annotation.getTranscript(transcriptID);
			if(tmp.getTotalCodingExonLength() > 0){
				count ++;
				_averageLength_CDS += tmp.getTotalCodingExonLength();
				_averageLength_UTR5 += tmp.getCDSStartFromTxStartInTxCoords()-1;
				_averageLength_UTR3 += tmp.getTotalExonLength() - (tmp.getTotalCodingExonLength() + tmp.getCDSStartFromTxStartInTxCoords() - 1);

				tmp_lengthsCDS.add(tmp.getTotalCodingExonLength());
				tmp_lengthsUTR5.add(tmp.getCDSStartFromTxStartInTxCoords()-1);
				tmp_lengthsUTR3.add(tmp.getTotalExonLength() - (tmp.getTotalCodingExonLength() + tmp.getCDSStartFromTxStartInTxCoords() - 1));
			}
		}
		
		// for MEAN lengths:
		//_averageLength_CDS = _averageLength_CDS / (count+0.0);
		//_averageLength_UTR5 = _averageLength_UTR5 / (count+0.0);
		//_averageLength_UTR3 = _averageLength_UTR3 / (count+0.0);

		// for MEDIAN lengths:
		_averageLength_CDS = Stats.median(tmp_lengthsCDS);
		_averageLength_UTR5 = Stats.median(tmp_lengthsUTR5);
		_averageLength_UTR3 = Stats.median(tmp_lengthsUTR3);

		double totalLength = _averageLength_CDS + _averageLength_UTR5 + _averageLength_UTR3;
		double basesPerBin = totalLength/(totalBins+0.0);
		_nBins_UTR5 = (int)Math.floor(_averageLength_UTR5/(basesPerBin+0.0));
		_nBins_UTR3 = (int)Math.floor(_averageLength_UTR3/(basesPerBin+0.0));
		_nBinsCDS = totalBins - (_nBins_UTR5 + _nBins_UTR3);

		IO_utils.printLineErr("Average annotated lengths:  5'UTR="+Math.floor(_averageLength_UTR5)+"  CDS="+Math.round(_averageLength_CDS)+"  3'UTR="+Math.floor(_averageLength_UTR3));
		IO_utils.printLineErr("nBins:  5'UTR="+_nBins_UTR5+"  CDS="+_nBinsCDS+"  3'UTR="+_nBins_UTR3);
	}


	private boolean _inputSourceExceRpt = false;
	private String reformatExceRptReferenceID(String id){
		String newID = null;
		if(id.startsWith("gencode"))
			newID = id.split(":")[1];
		return newID;
	}
	
	
	/**
	 * Reads *transcriptome* alignments in the given file
	 * @param alignmentFile
	 * @throws IOException 
	 */
	//private void readAlignments(File alignmentFile, int nBinsPerTranscript, boolean chooseOneTranscriptPerGene){
	//private void readAlignments(File alignmentFile, int nBinsPerTranscript, File output, String colDelim) throws IOException{
	private void readAlignments(File alignmentFile, File output, String colDelim) throws IOException{

		//_nBins_UTR5 = (int)Math.round(nBinsPerTranscript/20.0);
		//_nBins_UTR3 = (int)Math.round(nBinsPerTranscript/10.0);

		//IO_utils.printLineErr("Reading alignments...");
		SAMReader engine = new SAMReader(alignmentFile);

		//if(verbose)
		//	System.err.println("BAM sorted by readID: "+engine.isSortedByReadID());
		if( ! engine.isSortedByReadID()){
			System.err.println("BAM must be sorted by readID ('queryname')\nThis BAM is sorted by: "+engine.getSortOrder().toString());
			System.exit(0);
		}

		// read reference sequence lengths
		IO_utils.printLineErr("Reading reference sequences (transcripts)...");
		HashMap<String, Integer> referenceID2Length = engine.getReferenceLengths();
		Iterator<String> it = referenceID2Length.keySet().iterator();
		while(it.hasNext()){
			String refID = it.next();
			
			// if this BAM is from exceRpt, we need to modify the reference IDs to remove the extra crap
			if(_inputSourceExceRpt){
				String newID = reformatExceRptReferenceID(refID);
				if(newID != null)
					refID = newID;
			}
			
			try{
				_referenceID2RefObj.put(refID, new TranscriptCoverage(refID, _nBinsCDS, _nBins_UTR5, _nBins_UTR3, _annotation.getTranscript(refID)));
			}catch(NullPointerException e){
				//IO_utils.printLineErr("WARNING: unable to locate transcript \'"+refID+"\' in the annotation");				
			}
		}
		IO_utils.printLineErr("N reference sequences (transcripts): "+referenceID2Length.size());


		IO_utils.printLineErr("Adding read alignments...");
		ArrayList<SAMRecordReduced> thisRead;
		// Loop through all reads
		while((thisRead = engine.getAlignmentsForNextRead()) != null){
			for(SAMRecordReduced alignment: thisRead){
				
				String refID = alignment.getReferenceName();
				// if this BAM is from exceRpt, we need to modify the reference IDs to remove the extra crap
				if(_inputSourceExceRpt){
					String newID = reformatExceRptReferenceID(refID);
					if(newID != null)
						refID = newID;
				}
				
				if(_referenceID2RefObj.containsKey(refID)){
					int readMid = (int) Math.floor(alignment.getAlignmentStart()+((alignment.getAlignmentEnd()-alignment.getAlignmentStart())/2.0));
					_referenceID2RefObj.get(refID).addRead(readMid);
				}
			}
		}

		IO_utils.printLineErr("Writing read bins...");
		BufferedWriter out = new BufferedWriter(new FileWriter(output));
		out.write("GeneID"+colDelim+"TranscriptID"+colDelim+"TranscriptBiotype"+colDelim+"Length_Transcript"+colDelim+"Length_CDS"+colDelim+"ReadCount_Transcript"+colDelim+"ReadCount_UTR5"+colDelim+"ReadCount_CDS"+colDelim+"ReadCount_UTR3"+colDelim+"Entropy_Transcript"+colDelim+"Entropy_CDS");

		//out.write(outputColumnDelimiter+"utr5_0");
		for(int i=0;i<_nBins_UTR5;i++)
			out.write(colDelim+"utr5_"+i);

		//out.write(outputColumnDelimiter+"cds_0");
		for(int i=0;i<_nBinsCDS;i++)
			out.write(colDelim+"cds_"+i);

		//out.write(outputColumnDelimiter+"utr3_0");
		for(int i=0;i<_nBins_UTR3;i++)
			out.write(colDelim+"utr3_"+i);

		out.write("\n");

		it = _annotation.getMap_gene2transcript().keySet().iterator();
		//it = _referenceID2RefObj.keySet().iterator();
		while(it.hasNext()){
			String geneID = it.next();
			ArrayList<String> transcripts = _annotation.getMap_gene2transcript().get(geneID);
			for(String transcript: transcripts){
				if(_referenceID2RefObj.containsKey(transcript)){
					TranscriptCoverage thisTranscript = _referenceID2RefObj.get(transcript);
					if(thisTranscript.getNReads() > 0){
						String cdsEntropy = "NA";
						if(thisTranscript.hasCDS())
							cdsEntropy = thisTranscript.calcCDSEntropy()+"";

						out.write(geneID +colDelim+ thisTranscript.getID() +colDelim+ thisTranscript.getBiotype() +colDelim+ thisTranscript.getLength_Transcript() +colDelim+ thisTranscript.getLength_CDS() +colDelim+
								thisTranscript.getNReads() +colDelim+ thisTranscript.getNReads_5p() +colDelim+ thisTranscript.getNReads_CDS() +colDelim+ thisTranscript.getNReads_3p() +colDelim+
								thisTranscript.calcTranscriptEntropy() +colDelim+ cdsEntropy +colDelim+ 
								thisTranscript.toString(colDelim)+"\n");
					}
				}
			}
		}
		out.flush();
		out.close();

		IO_utils.printLineErr("Done");
	}




	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		//options.addOption(new Option("verbose", "write the rejected sequences to stderr"));
		options.addOption(OptionBuilder.withArgName("annotationPath").hasArg().withDescription("Path to the GTF file containing the transcript/gene relationship").create(Thunder.OPT_PATH_ANNOTATION));
		options.addOption(OptionBuilder.withArgName("read alignments").hasArg().withDescription("SAM/BAM file containing alignments against the transcriptome").create("f"));
		options.addOption(OptionBuilder.withLongOpt("outputPrefix").withArgName("outputPath").hasArg().withDescription("[optional] Output prefix for results [default: same as input file]").create("o"));
		options.addOption(OptionBuilder.withArgName("nBins").hasArg().withDescription("[optional] Transcript must have more than this minimum number of footprints to be included [default: 200]").create("N"));
		options.addOption(OptionBuilder.withDescription("[optional] The alignments are output from the exceRpt pipeline").create("exceRpt"));
		return options;
	}

	public static void main(String[] args) throws Exception {
		/*String annotation = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/MOUSE/gencode.vM6.primary_assembly.annotation.gtf";
		String transcriptomeMappedReads = "/Users/robk/Downloads/D2-S1_Genome_Aligned.toTranscriptome.sorted.bam";
		String outputPath = "/Users/robk/Downloads/D2-S1_Genome_Aligned.toTranscriptome.sorted.bam.bins.txt";
		int nBins = 200;*/
		//boolean chooseOneTranscriptPerGene = false;
		
		/*args = new String[]{
				"-f","/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Footprinting/FP_original_2_Genome_Aligned.toTranscriptome.sorted.bam",
				"-"+Thunder.OPT_PATH_ANNOTATION,"/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation_noSelenocysteine.gtf"};
		*/
		
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption(Thunder.OPT_PATH_ANNOTATION)  &&  cmdArgs.hasOption("f")){

			String outPath = cmdArgs.getOptionValue("f")+".coverage.txt";
			if(cmdArgs.hasOption("o"))
				outPath = cmdArgs.getOptionValue("o");
			
			int nBins = 200;	
			if(cmdArgs.hasOption("o"))
				nBins = Integer.valueOf(cmdArgs.getOptionValue("N")).intValue();
			
			ReadCoverage engine = new ReadCoverage(cmdArgs.getOptionValue(Thunder.OPT_PATH_ANNOTATION), nBins);
			
			// do we need to modify the reference IDs (from the exceRpt output)
			if(cmdArgs.hasOption("exceRpt"))
				engine._inputSourceExceRpt = true;
			
			//engine.readAlignments(new File(transcriptomeMappedReads), nBins, chooseOneTranscriptPerGene);
			engine.readAlignments(new File(cmdArgs.getOptionValue("f")), new File(outPath), ",");
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ReadCoverage", getCmdLineOptions());
			System.err.println();
		}
	}
}





class TranscriptCoverage{
	private String _id;
	Transcript _annotation;
	private int _lengthExon, _lengthCDS, _length5p, _length3p;
	private double _basesPerBin_Exon, _basesPerBin_CDS, _basesPerBin_5p, _basesPerBin_3p;
	private int _nReads_ALL=0, _nReads_CDS=0, _nReads_5p=0, _nReads_3p=0;
	private int[] _binsExon, _binsCDS, _bins5p, _bins3p;
	private HashMap<String, Integer[]> _readPositions = new HashMap<String, Integer[]>();

	TranscriptCoverage(String id, int nBins_CDS, int nBins_utr5, int nBins_utr3, Transcript transcriptAnnotation){

		_id = id;
		_annotation = transcriptAnnotation;
		_binsExon = new int[nBins_CDS+nBins_utr5+nBins_utr3];
		_binsCDS = new int[nBins_CDS];
		_bins5p = new int[nBins_utr5];
		_bins3p = new int[nBins_utr3];

		_lengthExon = transcriptAnnotation.getTotalExonLength();
		_lengthCDS = transcriptAnnotation.getTotalCodingExonLength();
		_length5p = transcriptAnnotation.getCDSStartFromTxStartInTxCoords();
		_length3p = transcriptAnnotation.getTotalExonLength()-(_length5p+_lengthCDS);

		_basesPerBin_Exon = (_lengthExon+0.0)/(nBins_CDS+nBins_utr5+nBins_utr3+0.0);
		_basesPerBin_CDS = (_lengthCDS+0.0)/(nBins_CDS+0.0);
		_basesPerBin_5p = (_length5p+0.0)/(nBins_utr5);
		_basesPerBin_3p = (_length3p+0.0)/(nBins_utr3+0.0);

		_readPositions.put("exon", Object_utils.initArray_Integer(_lengthExon));
		if(_lengthCDS > 0){
			_readPositions.put("UTR5", Object_utils.initArray_Integer(_length5p));
			_readPositions.put("CDS", Object_utils.initArray_Integer(_lengthCDS));
			_readPositions.put("UTR3", Object_utils.initArray_Integer(_length3p));
		}
	}


	/**
	 * 
	 * @param position
	 */
	void addRead(int position){
		_nReads_ALL ++;
		_binsExon[(int) Math.floor(position/_basesPerBin_Exon)] ++;
		_readPositions.get("exon")[position-1] += 1;

		if(_lengthCDS > 0){
			int cdsStart = _annotation.getCDSStartFromTxStartInTxCoords();

			if(position < cdsStart-1){
				_nReads_5p ++;
				_bins5p[(int) Math.floor(position/_basesPerBin_5p)] ++;
				_readPositions.get("UTR5")[position-1] += 1;
			}else if(position >= (cdsStart+_lengthCDS+1)){
				try{
					_nReads_3p ++;
					_bins3p[(int) Math.floor((position-(cdsStart+_lengthCDS+1))/_basesPerBin_3p)] ++;
					_readPositions.get("UTR3")[(position-(cdsStart+_lengthCDS))] += 1;
				}catch(ArrayIndexOutOfBoundsException e){ 
					System.out.println("Error in transcript:"+_id+" position:"+position+" cdsStart:"+cdsStart+" cdsLength:"+_lengthCDS+" exonLength:"+_lengthExon); 
				}
			}else if(position >= cdsStart  &&  position < (cdsStart+_lengthCDS)){
				_nReads_CDS ++;
				_binsCDS[(int) Math.floor((position-cdsStart)/_basesPerBin_CDS)] ++;
				_readPositions.get("CDS")[(position-cdsStart)] += 1;
			}
		}
	}

	String getID(){ return _id; }
	String getSymbol(){ return _annotation.getTranscriptSymbol(); }
	String getBiotype(){ return _annotation.getTranscriptBiotype(); }
	boolean hasCDS(){ return _lengthCDS > 0; }
	int getNReads(){ return _nReads_ALL; }
	int getNReads_CDS(){ return _nReads_CDS; }
	int getNReads_5p(){ return _nReads_5p; }
	int getNReads_3p(){ return _nReads_3p; }
	int[] getBins_CDS(){ return _binsCDS; }
	int[] getBins_5p(){ return _bins5p; }
	int[] getBins_3p(){ return _bins3p; }
	int getLength_Transcript(){ return _lengthExon; }
	int getLength_CDS(){ return _lengthCDS; }

	double calcTranscriptEntropy(){
		Integer[] readCounts = _readPositions.get("exon");
		double[] positionProbs = new double[readCounts.length];
		for(int i=0;i<readCounts.length;i++)
			positionProbs[i] = (readCounts[i]+0.0)/(_nReads_ALL+0.0);

		double expectedProbs = 1.0/(_lengthExon+0.0);
		if(_nReads_ALL < _lengthExon)
			expectedProbs = 1.0/(_nReads_ALL+0.0);

		return Stats.calculateKS(positionProbs);
		//return Stats.calculateKL(positionProbs, expectedProbs);
		//return Stats.calculateEntropy(positionProbs);
	}

	double calcCDSEntropy(){
		Integer[] readCounts = _readPositions.get("CDS");
		double[] positionProbs = new double[readCounts.length];
		for(int i=0;i<readCounts.length;i++)
			positionProbs[i] = (readCounts[i]+0.0)/(_nReads_CDS+0.0);

		double expectedProbs = 1.0/(_lengthCDS+0.0);
		if(_nReads_CDS < _lengthCDS)
			expectedProbs = 1.0/(_nReads_CDS+0.0);

		return Stats.calculateKS(positionProbs);
		//return Stats.calculateKL(positionProbs, expectedProbs);
		//return Stats.calculateEntropy(positionProbs);
	}

	/**
	 * 
	 * @param sep
	 * @return
	 */
	String toString(String sep){
		if(_lengthCDS > 0){
			String out = _bins5p[0]+"";
			for(int i=1;i<_bins5p.length;i++)
				out = out.concat(sep+_bins5p[i]);

			//String out = _binsCDS[0]+"";
			for(int i=0;i<_binsCDS.length;i++)
				out = out.concat(sep+_binsCDS[i]);

			for(int i=0;i<_bins3p.length;i++)
				out = out.concat(sep+_bins3p[i]);

			return out;
		}else{
			//if this is a non-coding transcript
			String out = _binsExon[0]+"";
			for(int i=1;i<_binsExon.length;i++)
				out = out.concat(sep+_binsExon[i]);

			return out;
		}
	}
}

