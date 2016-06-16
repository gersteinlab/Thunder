package samTools;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Iterator;
import java.util.List;

import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMSequenceRecord;
import objects.SAMRecordReduced;

public class SAMReader {

	SAMFileReader _inputSAM;
	Iterator<SAMRecord> _recordIterator;
	public SAMReader(File inputFile){
		//SAMFileReader.setDefaultValidationStringency(ValidationStringency.LENIENT);
		_inputSAM = new SAMFileReader(inputFile);
		_inputSAM.setValidationStringency(ValidationStringency.LENIENT);
		_recordIterator = _inputSAM.iterator();
	}
	
	
	/**
	 * Get the nucleotide lengths of all *reference* sequences in this BAM 
	 * @return
	 */
	public HashMap<String, Integer> getReferenceLengths(){
		HashMap<String, Integer> referenceLengths = new HashMap<String, Integer>();
		List<SAMSequenceRecord> references = _inputSAM.getFileHeader().getSequenceDictionary().getSequences();
		for(SAMSequenceRecord thisReference: references){
			referenceLengths.put(thisReference.getSequenceName(), thisReference.getSequenceLength());
		}
		return referenceLengths;
	}

	/**
	 * Is this file sorted by readID?
	 * @return
	 */
	public boolean isSortedByReadID(){
		return _inputSAM.getFileHeader().getSortOrder().compareTo(SortOrder.queryname) == 0;
	}

	/**
	 * Get the sort order of this SAM/BAM file
	 * @return
	 */
	public SortOrder getSortOrder(){
		return _inputSAM.getFileHeader().getSortOrder();
	}

	
	/**
	 * Variables for reading and keeping track of alignments
	 */
	private SAMRecord _currentRecord = null;
	boolean _EOF = false;
	
	/**
	 * Reads the next read alignment (if there is one)
	 * @return
	 */
	public ArrayList<SAMRecordReduced> getAlignmentsForNextRead(){

		if(_EOF){
			return null;

		}else{

			ArrayList<SAMRecordReduced> alignments = new ArrayList<SAMRecordReduced>();
			SAMRecord thisAlignment;

			// if there is an alignment for the previous read, add it.
			if(_currentRecord != null)
				alignments.add(new SAMRecordReduced(_currentRecord));

			boolean keepReading = true;
			while(keepReading  &&  _recordIterator.hasNext()){
				thisAlignment = _recordIterator.next();

				if(_currentRecord != null){
					// this is not the first read

					if(thisAlignment.getReadName().equals(_currentRecord.getReadName())){
						// readID matches the last one, keep adding
						alignments.add(new SAMRecordReduced(thisAlignment));
					}else{
						// this is a new read ID
						keepReading = false;
					}
				}else{
					// this is the first read
					alignments.add(new SAMRecordReduced(thisAlignment));
				}
				_currentRecord = thisAlignment;
			}

			// if this is the last end of the file, return the current read but set the flag so no more reads are read
			if( ! _recordIterator.hasNext())
				_EOF = true;

			return alignments;
		}
	}
	
	/**
	 * Close the reader
	 */
	public void close(){
		_inputSAM.close();
	}
	
	
	/**
	 * Compute the number of distinct genic loci mapped by this gene
	 * @param thisRead
	 * @param transcriptID_2_geneID
	 * @return
	 */
	public static int getNumberOfGenesForRead(ArrayList<SAMRecordReduced> thisRead, HashMap<String, String> transcriptID_2_geneID){

		HashSet<String> geneIDs = new HashSet<String>();
		String thisTranscriptID = "";

		// Loop through all alignments for this read to identify if this multimaps across more than one gene
		for(SAMRecordReduced thisAlignment: thisRead){

			// get the transcriptID for this alignment:
			thisTranscriptID = thisAlignment.getReferenceName();

			// if this transcript has a gene in the annotation:
			if(transcriptID_2_geneID.containsKey(thisTranscriptID)){
				String thisGeneID = transcriptID_2_geneID.get(thisTranscriptID);

				// add the geneID to the list, if more than one gene we can stop later as this is a multimapper 
				geneIDs.add(thisGeneID);
			}else{
				System.err.println("ERROR, something went really wrong: transcript \'"+thisTranscriptID+"\' does not appear in the annotation.");
			}
		}
		
		return geneIDs.size();
	}
	
	
	/**
	 * Does this gene gave a unique alignment to the genome?
	 * @param thisRead
	 * @param transcriptID_2_geneID
	 * @return
	 */
	public static boolean isReadUniquelyMappedToOneGene(ArrayList<SAMRecordReduced> thisRead, HashMap<String, String> transcriptID_2_geneID){
		if(getNumberOfGenesForRead(thisRead, transcriptID_2_geneID) == 1)
			return true;
		else
			return false;
	}
	


	public static void main(String[] args) throws IOException {
		// TODO Auto-generated method stub
		

		SAMReader engine = new SAMReader(new File("/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Footprinting/FP_original_1_Genome_Aligned.toTranscriptome.sorted.bam"));
		System.out.println("Sorted by readID: "+engine.isSortedByReadID());

		/*
		
		ArrayList<SAMRecord> alignments;
		System.out.println();
		alignments = engine.getAlignmentsForNextRead();
		for(SAMRecord thisAlignment: alignments)
			System.out.println(thisAlignment.getReadName()+"\t"+thisAlignment.getReferenceName());

		System.out.println();
		alignments = engine.getAlignmentsForNextRead();
		for(SAMRecord thisAlignment: alignments)
			System.out.println(thisAlignment.getReadName()+"\t"+thisAlignment.getReferenceName());

		System.out.println();
		alignments = engine.getAlignmentsForNextRead();
		for(SAMRecord thisAlignment: alignments)
			System.out.println(thisAlignment.getReadName()+"\t"+thisAlignment.getReferenceName());
		 */

		BufferedWriter out = new BufferedWriter(new FileWriter(new File("/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Footprinting/OUT_TEST.txt")));
		ArrayList<SAMRecordReduced> thisRead;
		int readCount = 0;
		// Loop through all reads
		while((thisRead = engine.getAlignmentsForNextRead()) != null){
			readCount ++;
			out.write("\n");
			for(SAMRecordReduced thisAlignment: thisRead)
				out.write(thisAlignment.getReadName()+"\t"+thisAlignment.getReferenceName()+"\n");
		}
		out.flush();
		out.close();
		System.out.println("N reads: "+readCount);

		engine.close();
	}

}
