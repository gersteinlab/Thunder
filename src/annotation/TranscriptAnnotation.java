package annotation;

import java.util.ArrayList;
import java.util.HashMap;

import objects.GenomicCoordinate;
import objects.Transcript;


public class TranscriptAnnotation {

	private HashMap<String, String> _transcript2gene = new HashMap<String, String>();
	private HashMap<String, ArrayList<String>> _gene2transcript = new HashMap<String, ArrayList<String>>();
	private HashMap<String, Transcript> _transcripts = new HashMap<String, Transcript>();
	
	TranscriptAnnotation(){}
	
	/**
	 * 
	 * @param chr
	 * @param source
	 * @param strand
	 * @param geneID
	 * @param transcriptID
	 */
	private void addTranscript(String chr, String source, String strand, String geneID, String transcriptID, String transcriptBiotype){
		if(!_transcripts.containsKey(transcriptID)){
			_transcripts.put(transcriptID, new Transcript(transcriptID, chr, strand, transcriptBiotype, geneID));
			_transcript2gene.put(transcriptID, geneID);
		}
		if(!_gene2transcript.containsKey(geneID))
			_gene2transcript.put(geneID, new ArrayList<String>());
		if(!_gene2transcript.get(geneID).contains(transcriptID))
			_gene2transcript.get(geneID).add(transcriptID);
	}
	
	
	/**
	 * 
	 * @param chr
	 * @param source
	 * @param start
	 * @param stop
	 * @param strand
	 * @param geneID
	 * @param transcriptID
	 */
	public void addExon(String chr, String source, GenomicCoordinate coords, String strand){//, String geneID, String transcriptID, String transcriptBiotype){
		addTranscript(chr, source, strand, coords.getAttribute("gene_id"), coords.getAttribute("transcript_id"), coords.getAttribute("transcript_type"));
		_transcripts.get(coords.getAttribute("transcript_id")).addExon(coords);
	}
	
	/**
	 * 
	 * @param chr
	 * @param source
	 * @param start
	 * @param stop
	 * @param strand
	 * @param geneID
	 * @param transcriptID
	 */
	public void addCDS(String chr, String source, GenomicCoordinate coords, String strand){//, String geneID, String transcriptID, String transcriptBiotype){
		//addTranscript(chr, source, strand, geneID, transcriptID, transcriptBiotype);
		//_transcripts.get(transcriptID).addCDS(start, stop);
		
		addTranscript(chr, source, strand, coords.getAttribute("gene_id"), coords.getAttribute("transcript_id"), coords.getAttribute("transcript_type"));
		_transcripts.get(coords.getAttribute("transcript_id")).addCDS(coords);
	}
	
	/**
	 * Access transcript info
	 * @return
	 */
	public HashMap<String, Transcript> getTranscripts(){ return _transcripts; }
	public Transcript getTranscript(String transcriptID){ return _transcripts.get(transcriptID); }
	
	
	/**
	 * Access transcript -> gene ID map
	 * @return
	 */
	public HashMap<String, String> getMap_transcript2gene(){ return _transcript2gene; }
	public String getGeneForTranscript(String transcriptID){ return _transcript2gene.get(transcriptID); }
	
	/**
	 * Access gene -> transcript ID map
	 * @return
	 */
	public HashMap<String, ArrayList<String>> getMap_gene2transcript(){ return _gene2transcript; }
	public ArrayList<String> getTranscriptsForGene(String geneID){ return _gene2transcript.get(geneID); }
	
	
}





