package expectationMaximisation;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

public class EM_Gene{
	private String _id;
	private HashMap<String, EM_Transcript> _transcripts = new HashMap<String, EM_Transcript>();

	public EM_Gene(String id){ _id = id; }
	public String getID(){ return _id; }

	public boolean containsTranscript(String id){ return(_transcripts.containsKey(id)); }
	public void addTranscript(String id){ 
		_transcripts.put(id, new EM_Transcript(id));
		if(!_transcriptReadCounts.containsKey(id))
			_transcriptReadCounts.put(id, 0);
	}
	
	private boolean _restrictToCDS = false;
	public void setRestrictToCDS(boolean val){ _restrictToCDS = val; }
	public boolean IsRestrictedToCDS(){ return _restrictToCDS; }

	private boolean _hasCDS = false;
	public void setProteinCoding(boolean val){ _hasCDS = val; }
	public boolean hasCDS(){ return _hasCDS; }

	public int transcriptCount(){ return(_transcripts.size()); }
	public EM_Transcript getTranscript(String id){ 
		return(_transcripts.get(id)); 
	}


	/**
	 * 
	 * @return
	 */
	public ArrayList<String> getTranscriptIDs(){
		ArrayList<String> ids = new ArrayList<String>();
		Iterator<String> it = _transcripts.keySet().iterator();
		while(it.hasNext()){
			ids.add(it.next());
		}
		Collections.sort(ids);
		return(ids);
	}


	/**
	 * 
	 */
	public void setPriors(){
		double totalTPM = 0.0;
		Iterator<String> it = _transcripts.keySet().iterator();
		while(it.hasNext()){
			totalTPM += _transcripts.get(it.next()).getTPM();
		}
		it = _transcripts.keySet().iterator();
		String id = "";
		while(it.hasNext()){
			id = it.next();
			if(totalTPM == 0.0){
				_transcripts.get(id).setPrior(1.0 / (_transcripts.size()+0.0));
				_hasFlatPrior = true;
			}else{
				_transcripts.get(id).setPrior(_transcripts.get(id).getTPM() / totalTPM);
				_hasFlatPrior = false;
			}
		}
	}

	private boolean _hasFlatPrior = true;
	public boolean hasFlatPrior(){ return _hasFlatPrior; }
	
	
	
	/**
	 * 
	 * @return
	 */
	public HashMap<String, BigDecimal> getPriors(){
		HashMap<String, BigDecimal> priors = new HashMap<String, BigDecimal>();
		Iterator<String> it = _transcripts.keySet().iterator();
		String id;
		while(it.hasNext()){
			id = it.next();
			priors.put(id, new BigDecimal(_transcripts.get(id).getPrior()));
		}
		return priors;
	}


	/**
	 * 
	 * @return
	 */
	public HashMap<String, Double> getAbundances(){
		HashMap<String, Double> expressions = new HashMap<String, Double>();
		Iterator<String> it = _transcripts.keySet().iterator();
		String id;
		while(it.hasNext()){
			id = it.next();
			expressions.put(id, _transcripts.get(id).getTPM());
		}
		return expressions;
	}

	/**
	 * 
	 * @return
	 */
	public HashMap<String, Integer> getLengths(){
		HashMap<String, Integer> lengths = new HashMap<String, Integer>();
		Iterator<String> it = _transcripts.keySet().iterator();
		String id;
		while(it.hasNext()){
			id = it.next();
			if(!_restrictToCDS)
				lengths.put(id, _transcripts.get(id).getExonLength());
			else
				lengths.put(id, _transcripts.get(id).getCodingExonLength());
		}
		return lengths;
	}



	private HashMap<String, ArrayList<String>> _reads = new HashMap<String, ArrayList<String>>();
	private HashMap<String, Integer> _readCounts = new HashMap<String, Integer>();
	private HashMap<String, Integer> _transcriptReadCounts = new HashMap<String, Integer>();
	
	public void addRead(String readID, ArrayList<String> hitsTranscripts){
		if(!_reads.containsKey(readID)){
			_reads.put(readID, hitsTranscripts); 
			_readCounts.put(readID, 1);
		
			// Count reads per transcript
			for(int i=0;i<hitsTranscripts.size();i++){
				String tmpID = hitsTranscripts.get(i);
				_transcriptReadCounts.put(tmpID, _transcriptReadCounts.get(tmpID)+1);
			}	
		}
	}

	
	/**
	 * Add a read and alignments (meant for situations where a 'read' is more like a peptide in that it has some intensity value, rather than being an atomic unit
	 * @param readID
	 * @param hitsTranscripts
	 * @param readCount
	 */
	public void addRead(String readID, ArrayList<String> hitsTranscripts, int readCount){ 
		if(!_reads.containsKey(readID)){
			_reads.put(readID, hitsTranscripts); 
			_readCounts.put(readID, readCount);
			
			// Count reads per transcript
			for(int i=0;i<hitsTranscripts.size();i++){
				String tmpID = hitsTranscripts.get(i);
				_transcriptReadCounts.put(tmpID, _transcriptReadCounts.get(tmpID)+1);
			}
		}else{
			_readCounts.put(readID, _readCounts.get(readID)+readCount);
		}
	}
	
	
	public HashMap<String, Integer> getTranscriptReadCounts(){ return _transcriptReadCounts; }
	public int getTranscriptReadCounts(String transcriptID){ return _transcriptReadCounts.get(transcriptID); }
	
	public HashMap<String, ArrayList<String>> getReads(){ return _reads; }
	public HashMap<String, Integer> getReadCounts(){ return _readCounts; }
	public int getReadCount(String readID){ return _readCounts.get(readID); }
}




