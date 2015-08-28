package expectationMaximisation;

import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

public class EM_Core_original {

	private boolean _verbose = false;

	private EM_Gene _gene;
	private ArrayList<String> _transcriptIDs;
	private boolean _hasFlatPrior = false;

	double[] _priors;//, _finalTranscriptRatios, _finalAbundances;
	int[] _transcriptLengths;

	private HashMap<String, ArrayList<String>> _reads = new HashMap<String, ArrayList<String>>();
	private Object[] _readIDs;

	private static final double _TINY = 1.0/10000.0;



	public EM_Core_original(EM_Gene gene){
		_gene = gene;
		_transcriptIDs = _gene.getTranscriptIDs();
		Collections.sort(_transcriptIDs);

		EM_Transcript tmp;

		_priors = new double[_transcriptIDs.size()];
		_transcriptLengths = new int[_transcriptIDs.size()];

		for(int i=0;i<_transcriptIDs.size();i++){
			tmp = _gene.getTranscript(_transcriptIDs.get(i));
			//System.out.println("transcript: "+tmp.getID());
			_priors[i] = tmp.getPrior();
			_transcriptLengths[i] = tmp.getExonLength();

			//if(!_hasFlatPrior)
			_hasFlatPrior = !tmp.hasRNAseq();
		}
		_reads = _gene.getReads();
		_readIDs = _reads.keySet().toArray();
	}



	public void runEM(int maxIterations, double converganceDistance, PrintWriter out_exprs){
		//StringBuilder sb = new StringBuilder();

		//if(_verbose){
		/*if(_gene.getID().equals("ENSG00000091483.6")){
			System.err.println("N isoforms: "+_gene.transcriptCount());
			System.err.println("N reads: "+_reads.size());
			Iterator<String> it = _reads.keySet().iterator();
			Iterator<String> it2;
			String id = "";
			while(it.hasNext()){
				id = it.next();
				System.err.print(id+"\t");
				it2 = _reads.get(id).iterator();
				while(it2.hasNext())
					System.err.print(it2.next()+"\t");
				System.err.println();
			}
			System.err.println("\n");
		}*/

		//System.out.println("geneID\ttranscriptID\thasFlatPrior\tprior\tflag\tnIterations\teffectiveReadCount\treadDensity\ttranPercent");

		if(_gene.transcriptCount() == 1){									// single isoform gene!
			//sb.append(_gene.getID()+"\t"+_transcriptIDs.get(0)+"\t"+_hasFlatPrior+"\t"+_priors[0]+"\t"+1+"\t"+"NA"+"\t"+_reads.size()+"\t"+_reads.size()/(_transcriptLengths[0]-_readLength+1.0)+"\t"+_priors[0]+"\n");
			//sb.append(_gene.getID()+"\t"+_transcriptIDs.get(0)+"\t"+_hasFlatPrior+"\t"+formatSci.format(_priors[0])+"\t"+1+"\t"+"NA"+"\t"+formatDec.format(_reads.size())+"\t"+formatSci.format(_reads.size()/(_transcriptLengths[0]-_readLength+1.0))+"\t"+formatSci.format(_priors[0])+"\n");
			//System.out.printf(_gene.getID()+"\t"+_transcriptIDs.get(0)+"\t"+_hasFlatPrior+"\t%e\t1\t0\t%f\t%e\t%e\n", _priors[0]+0.0, _reads.size()+0.0, _reads.size()/(_transcriptLengths[0]-_readLength+1.0), _priors[0]+0.0);
			out_exprs.printf(_gene.getID()+"\t"+_transcriptIDs.get(0)+"\t"+_hasFlatPrior+"\t%e\t1\t0\t%f\t%e\t%e\n", _priors[0]+0.0, _reads.size()+0.0, _reads.size()/(_transcriptLengths[0]-_readLength+1.0), _priors[0]+0.0);
		}else if(_reads.size() == 0){										// no reads!
			for(int k=0;k<_transcriptIDs.size();k++){
				//sb.append(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t"+formatSci.format(_priors[k])+"\t"+1+"\t"+"NA"+"\t"+"0.000000"+"\t"+"0.000000"+"\t"+formatSci.format(_priors[k])+"\n");
				//System.out.printf(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t%e\t1\tNA\t%f\t%e\t%e\n", _priors[k]+0.0, 0.0, 0.0, _priors[k]+0.0);
				out_exprs.printf(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t%e\t1\tNA\t%f\t%e\t%e\n", _priors[k]+0.0, 0.0, 0.0, _priors[k]+0.0);
			}
		}else{																// multi-isoform gene with read support, do EM!

			_mu = _priors;
			double[] mu_new;
			double dist = 0.0;
			int flag = 0;

			int count = 0;
			for(int i=0;i<maxIterations;i++){

				// E-step
				expectationStep();

				// M-step
				mu_new = this.maximisationStep();

				// Distance calculation
				dist = this.calcConvergence(_mu, mu_new);

				if(_verbose){
					System.err.print(_gene.getID()+"\tk="+i+"\tmu=");
					for(int k=0;k<_mu.length;k++)
						System.err.print(_mu[k]+",");
					System.err.print("\tmu_new=");
					for(int k=0;k<_mu.length;k++)
						System.err.print(mu_new[k]+",");
					System.err.println("\tdis="+dist);
				}

				// Stop?
				if(dist < converganceDistance){
					flag = 1;
					break;
				}
				double sum = 0.0;
				for(int k=0;k<mu_new.length;k++)
					sum += mu_new[k];
				if(sum <= 0.0){
					flag = -1;
					break;
				}

				// Update
				_mu = mu_new;
				count ++;
			}

			// Finished with this gene
			for(int k=0;k<_transcriptIDs.size();k++){
				//System.out.println("geneID\ttranscriptID\thasFlatPrior\tprior\tflag\tnIterations\teffectiveReadCount\treadDensity\ttranPercent");
				//sb.append(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t"+_priors[k]+"\t"+flag+"\t"+count+"\t"+_readCounts[k]+"\t"+_readCounts[k]/(_transcriptLengths[k]-_readLength+1.0)+"\t"+_mu[k]+"\n");
				//System.out.printf(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t%e\t"+flag+"\t"+count+"\t%f\t%e\t%e\n", _priors[k]+0.0, _readCounts[k]+0.0, _readCounts[k]/(_transcriptLengths[k]-_readLength+1.0), _mu[k]+0.0);
				out_exprs.printf(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t%e\t"+flag+"\t"+count+"\t%f\t%e\t%e\n", _priors[k]+0.0, _readCounts[k]+0.0, _readCounts[k]/(_transcriptLengths[k]-_readLength+1.0), _mu[k]+0.0);

				//System.out.println("DONE: "+_gene.getID()+"\t"+_transcriptIDs.get(k)+"\tflag="+flag+"\tk="+count+"\tmu="+_mu[k]+"\tcount="+_readCounts[k]+"\texp="+_readCounts[k]/(_transcriptLengths[k]-_readLength+1.0));
			}
		}

		//return(sb.toString());
	}

	public HashMap<String, String> sampleReadAssignments(HashMap<String, String> readMap){
		Iterator<String> reads = _reads.keySet().iterator();
		String thisRead;
		//ArrayList<String> theseTranscripts;
		//HashMap<String, String> readMap = new HashMap<String, String>();

		while(reads.hasNext()){
			thisRead = reads.next();
			//theseTranscripts = _reads.get(thisRead);
			//double[] transcriptFracs = new double[_reads.get(thisRead).size()];
			double[] transcriptFracs = new double[_transcriptIDs.size()];
			double transcriptSum = 0.0;

			if(_transcriptIDs.size() == 1){
				transcriptFracs[0] = 1.0;
				transcriptSum = 1.0;
			}else{
				for(int k=0;k<_transcriptIDs.size();k++){
					String thisTranscriptID = _transcriptIDs.get(k);
					//if(theseTranscripts.contains(thisTranscriptID))
					if(_reads.get(thisRead).contains(thisTranscriptID)){
						transcriptFracs[k] = _readCounts[k];
						transcriptSum += _readCounts[k]+0.0;
					}
					else{
						transcriptFracs[k] = 0.0;
					}
				}
			}

			for(int k=0;k<_transcriptIDs.size();k++){
				transcriptFracs[k] = transcriptFracs[k] / (transcriptSum+0.0);
			}

			double rand = Math.random();
			double cumsum = 0.0;
			for(int k=0;k<_transcriptIDs.size();k++){
				cumsum += transcriptFracs[k];
				if(cumsum >= rand){
					readMap.put(thisRead, _transcriptIDs.get(k));
					break;
				}
			}
		}
		return(readMap);
	}


	private int _readLength = 29;
	public void setReadLength(int readLength){ this._readLength = readLength; }

	private double[] mu2tau(double[] priors){
		double[] tau = new double[priors.length];

		if(priors.length != _transcriptLengths.length){
			System.err.println("mu2tau: input vectors are not of the same length\n");
			System.exit(1);
		}

		double tauSum = 0.0;
		for(int k=0;k<priors.length;k++){
			tau[k] = priors[k] * Math.max(0, _transcriptLengths[k]-_readLength);
			tauSum += tau[k]; 
		}

		if(Math.abs(tauSum - _TINY) < 0){
			System.err.println("mu_2_tau: sum of prior is too small, please check parameters\n");
			System.exit(1);
		}

		for(int k=0;k<priors.length;k++){
			tau[k] = tau[k]/tauSum;
		}

		return(tau);
	}


	private double[] allocateRead(String readID, double[] tau, double[] readCounts){
		double[] addVal = new double[tau.length];
		double sumVal = 0.0;
		for(int k=0;k<_transcriptIDs.size();k++){
			String thisTranscriptID = _transcriptIDs.get(k);
			if(_reads.get(readID).contains(thisTranscriptID)){
				addVal[k] = tau[k];
			}else{
				addVal[k] = 0.0;
			}
			sumVal += addVal[k];
		}

		if(sumVal > 0.0){
			//*^^^ split the read proportionally to the transcripts^^^*/
			for(int k=0;k<_transcriptIDs.size();k++){
				//readCounts[k] = readCounts[k] + addVal[k]/sumVal;
				readCounts[k] = readCounts[k] + ((addVal[k]/sumVal) * _gene.getReadCount(readID));
			}
		}

		return(readCounts);
	}


	// mu are the priors
	private double[] _mu, _tau, _readCounts;

	/**
	 * 
	 */
	public void expectationStep(){

		_tau = mu2tau(_mu);
		_readCounts = new double[_tau.length];

		for(int k=0; k<_reads.size(); k++){
			_readCounts = allocateRead((String)_readIDs[k], _tau, _readCounts);
		}
	}

	/**
	 * 
	 * @return
	 */
	public double[] maximisationStep(){
		double[] mu_new = new double[_mu.length];
		double s = 0.0;

		for(int k=0;k<_readCounts.length;k++){
			mu_new[k] = _readCounts[k]/Math.max(1, _transcriptLengths[k]-_readLength);
			s += mu_new[k];
		}
		if(s > 0.0){
			for(int k=0;k<_readCounts.length;k++){
				mu_new[k] = mu_new[k]/s;
			}
		}
		return(mu_new);
	}

	/**
	 * 
	 * @param mu_old
	 * @param mu_new
	 * @return
	 */
	public double calcConvergence(double[] mu_old, double[] mu_new){
		if(mu_old.length != mu_new.length){
			System.err.println("cover_distance: two transcripts expression percent not of the same length\n");
			System.exit(1);
		}
		double dis = 0.0;
		for(int k=0;k<mu_old.length;k++){
			dis += Math.abs(mu_old[k]-mu_new[k]);
		}
		return(dis);
	}
}
