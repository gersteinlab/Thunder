package expectationMaximisation;

import java.math.BigDecimal;
import java.util.HashMap;


public class EM_Result {
	
	private EM_Gene _gene;
	private int _nIterations;
	private String _flag;
	private HashMap<String, BigDecimal> _final_readCounts;
	private HashMap<String, BigDecimal> _final_transcriptLikelihoods;
	public EM_Result(EM_Gene gene, int nIterations, String flag, HashMap<String, BigDecimal> final_readCounts, HashMap<String, BigDecimal> final_transcriptLikelihoods){
		_gene = gene;
		_nIterations = nIterations;
		_flag = flag;
		_final_readCounts = final_readCounts;
		_final_transcriptLikelihoods = final_transcriptLikelihoods;
	}
	
	public EM_Gene getGene(){ return this._gene; }
	public int getNIterations(){ return this._nIterations; }
	public String getFlag(){ return this._flag; }
	public HashMap<String, BigDecimal> getFinalEffectiveReadCount(){ return this._final_readCounts; }
	public HashMap<String, BigDecimal> getFinalTranscriptLikelihoods(){ return this._final_transcriptLikelihoods; }
	
	//out_exprs.printf(_gene.getID()+"\t"+_transcriptIDs.get(0)+"\t"+_hasFlatPrior+"\t%e\t1\t0\t%f\t%e\t%e\n", _priors[0]+0.0, _reads.size()+0.0, _reads.size()/(_transcriptLengths[0]-_readLength+1.0), _priors[0]+0.0);
}
