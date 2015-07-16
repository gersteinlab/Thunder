package expectationMaximisation;

import java.util.HashMap;

public class EM_Transcript{
	private String _id;
	
	
	public EM_Transcript(String id){ _id = id; }
	public String getID(){ return _id; }
	
	private String _biotype = "";
	public void setBiotype(String biotype){ _biotype = biotype; }
	public String getBiotype(){ return(_biotype); }
	
	private int _transcriptLength = 0;
	private int _cdsLength = 0;
	//public void addExonLength(int nt){ _transcriptLength += nt; }
	//public void addCodingExonLength(int nt){ _cdsLength += nt; }
	public void setExonLength(int nt){ _transcriptLength = nt; }
	public void setCodingExonLength(int nt){ _cdsLength = nt; }
	public int getExonLength(){ return _transcriptLength; }
	public int getCodingExonLength(){ return _cdsLength; }
	
	private int _startPos_CDS, _endPos_CDS, _startPos_transcript, _endPos_transcript = 0;
	public void addExon(int start, int stop){ 
		_transcriptLength += (stop-start);
		if(start < _startPos_transcript)
			_startPos_transcript = start;
		if(stop> _endPos_transcript)
			_endPos_transcript = stop;
	}
	public int getTranscriptStart(){ return _startPos_transcript; }
	public int getTranscriptStop(){ return _endPos_transcript; }
	
	private boolean _hasCDS = false;
	public void addCDS(int start, int stop){ 
		_hasCDS = true;
		_cdsLength += (stop-start);
		if(start < _startPos_CDS)
			_startPos_CDS = start;
		if(stop> _endPos_CDS)
			_endPos_CDS = stop;
	}
	public int getCDSStart(){ return _startPos_CDS; }
	public int getCDSStop(){ return _endPos_CDS; }
	public boolean hasCDS(){ return _hasCDS; }
	
	
	private boolean _hasRNAseq = false;
	public boolean hasRNAseq(){ return _hasRNAseq; }
	
	private double _TPM = 0.0;
	public void setTPM(double tpm){ _TPM = tpm; if(tpm > 0.0){ _hasRNAseq = true; } }
	public double getTPM(){ return(_TPM); }
	
	private double _prior = 0.0;
	public void setPrior(double prior){ _prior = prior; }
	public double getPrior(){ return(_prior); }
	
	
	// For the MS/MS EM:
	private boolean _isMSMS = false;
	public boolean isMassSpec(){ return _isMSMS; }
	private HashMap<Integer,Integer> _peptideCounts;
	public Integer getPeptideCount(int frame){ return(_peptideCounts.get(frame)); };
	public void setPeptideDigest(HashMap<Integer,Integer> counts){ _peptideCounts = counts;    _isMSMS = true; }
}
