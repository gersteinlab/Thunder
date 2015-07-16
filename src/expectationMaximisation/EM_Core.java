package expectationMaximisation;

import java.math.BigDecimal;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

public class EM_Core {

	
	public static void main(String[] args){
		EM_Gene gene = new EM_Gene("robbles");
		gene.addTranscript("A");
		gene.addTranscript("B");
		gene.addTranscript("C");
		
		gene.getTranscript("A").addCDS(0, 1000); //gene.getTranscript("A").addExon(0, 1000);
		gene.getTranscript("B").addCDS(0, 900); //gene.getTranscript("B").addExon(0, 900);
		gene.getTranscript("C").addCDS(0, 600); //gene.getTranscript("C").addExon(0, 600);
		
		gene.setRestrictToCDS(true);
		
		ArrayList<String> tmp;
		tmp = new ArrayList<String>(); tmp.add("A");tmp.add("B");tmp.add("C");
		gene.addRead("r1", tmp);
		tmp = new ArrayList<String>(); tmp.add("A");tmp.add("B");
		gene.addRead("r2", tmp);
		tmp = new ArrayList<String>(); tmp.add("A");tmp.add("B");tmp.add("C");
		gene.addRead("r3", tmp);
		tmp = new ArrayList<String>(); tmp.add("A");tmp.add("B");tmp.add("C");
		gene.addRead("r4", tmp);
		
		//gene.setPriors();
		gene.getTranscript("A").setPrior(0.95);
		gene.getTranscript("B").setPrior(0.01);
		gene.getTranscript("C").setPrior(0.04);
		
		EM_Core engine = new EM_Core(gene);
		engine.runEM(10, 0.00001);
	}
	
	
	private boolean _verbose = false;

	private EM_Gene _gene;
	//private ArrayList<String> _transcriptIDs;
	//private boolean _hasFlatPrior = false;

	//double[] _priors;//, _finalTranscriptRatios, _finalAbundances;
	//int[] _transcriptLengths;

	//private HashMap<String, ArrayList<String>> _reads = new HashMap<String, ArrayList<String>>();
	//private Object[] _readIDs;

	private static final double _TINY = 1.0/10000.0;


	// mu are the priors
	//private double[] _mu, _tau, _readCounts;
	private HashMap<String, BigDecimal> _mu, _tau, _readCounts;
	//private HashMap<String, Integer> _transcriptLengths;
	

	public EM_Core(EM_Gene gene){
		_gene = gene;
	}



	public EM_Result runEM(int maxIterations, double converganceDistance){//, PrintWriter out_exprs){
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
			
			//String thisTranscriptID = _gene.getTranscriptIDs().get(0);
			//out_exprs.printf(_gene.getID()+"\t"+thisTranscriptID+"\t"+_hasFlatPrior+"\t%e\t1\t0\t%f\t%e\t%e\n", _gene.getPriors().get(thisTranscriptID)+0.0, _gene.getReads().size()+0.0, _gene.getReads().size()/(_gene.getLengths().get(thisTranscriptID)-_readLength+1.0), _gene.getPriors().get(thisTranscriptID)+0.0);
			HashMap<String,BigDecimal> tmpReadCount = new HashMap<String,BigDecimal>();
			tmpReadCount.put(_gene.getTranscriptIDs().get(0), new BigDecimal(_gene.getReads().size()+0.0));
			return(new EM_Result(_gene, 0, "1", tmpReadCount, _gene.getPriors()));
		}else if(_gene.getReads().size() == 0){										// no reads!
			HashMap<String,BigDecimal> tmpReadCount = new HashMap<String,BigDecimal>();
			for(int k=0;k<_gene.getTranscriptIDs().size();k++){
				tmpReadCount.put(_gene.getTranscriptIDs().get(k), new BigDecimal(0.0));
			}
			return(new EM_Result(_gene, 0, "NA", tmpReadCount, _gene.getPriors()));
			
			/*for(int k=0;k<_gene.getTranscriptIDs().size();k++){
				String thisTranscriptID = _gene.getTranscriptIDs().get(k);
				//sb.append(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t"+formatSci.format(_priors[k])+"\t"+1+"\t"+"NA"+"\t"+"0.000000"+"\t"+"0.000000"+"\t"+formatSci.format(_priors[k])+"\n");
				//System.out.printf(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t%e\t1\tNA\t%f\t%e\t%e\n", _priors[k]+0.0, 0.0, 0.0, _priors[k]+0.0);
				out_exprs.printf(_gene.getID()+"\t"+thisTranscriptID+"\t"+_hasFlatPrior+"\t%e\t1\tNA\t%f\t%e\t%e\n", _gene.getPriors().get(thisTranscriptID)+0.0, 0.0, 0.0, _gene.getPriors().get(thisTranscriptID)+0.0);
			}*/
		}else{																// multi-isoform gene with read support, do EM!

			_mu = _gene.getPriors();
			HashMap<String,BigDecimal> mu_new;
			double dist = 0.0;
			String flag = "0";

			int count = 0;
			for(int i=0;i<maxIterations;i++){

				// E-step
				expectationStep();

				// M-step
				mu_new = maximisationStep();

				// Distance calculation
				dist = calcConvergence(_mu, mu_new);

				if(_verbose){
					System.err.println(_gene.getID()+"\tk="+i+"\n\tmu:");
					for(int k=0;k<_gene.getTranscriptIDs().size();k++){
						String thisTranscriptID = _gene.getTranscriptIDs().get(k);
						System.err.println("\t"+thisTranscriptID+": "+_mu.get(thisTranscriptID));
					}
					
					System.err.println("\treads:");
					for(int k=0;k<_gene.getTranscriptIDs().size();k++){
						String thisTranscriptID = _gene.getTranscriptIDs().get(k);
						System.err.println("\t"+thisTranscriptID+": "+_readCounts.get(thisTranscriptID));
					}
					System.err.println("\tdis="+dist+"\n");
				}

				count ++;
				// Update
				this._mu = mu_new;
				
				
				// Stop?
				if(dist < converganceDistance){
					flag = "1";
					break;
				}
				double sum = 0.0;
				Iterator<String> it = mu_new.keySet().iterator();
				while(it.hasNext())
					sum += mu_new.get(it.next()).doubleValue();
				if(sum <= 0.0){
					flag = "-1";
					break;
				}


				
			}

			// Finished with this gene
			return(new EM_Result(_gene, count, flag, _readCounts, _mu));
			
			/*for(int k=0;k<_gene.getTranscriptIDs().size();k++){
				String thisTranscriptID = _gene.getTranscriptIDs().get(k);
				//System.out.println("geneID\ttranscriptID\thasFlatPrior\tprior\tflag\tnIterations\teffectiveReadCount\treadDensity\ttranPercent");
				//sb.append(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t"+_priors[k]+"\t"+flag+"\t"+count+"\t"+_readCounts[k]+"\t"+_readCounts[k]/(_transcriptLengths[k]-_readLength+1.0)+"\t"+_mu[k]+"\n");
				//System.out.printf(_gene.getID()+"\t"+_transcriptIDs.get(k)+"\t"+_hasFlatPrior+"\t%e\t"+flag+"\t"+count+"\t%f\t%e\t%e\n", _priors[k]+0.0, _readCounts[k]+0.0, _readCounts[k]/(_transcriptLengths[k]-_readLength+1.0), _mu[k]+0.0);
				out_exprs.printf(_gene.getID()+"\t"+thisTranscriptID+"\t"+_hasFlatPrior+"\t%e\t"+flag+"\t"+count+"\t%f\t%e\t%e\n", _gene.getPriors().get(thisTranscriptID)+0.0, _readCounts.get(thisTranscriptID)+0.0, _readCounts.get(thisTranscriptID)/(_gene.getLengths().get(thisTranscriptID)-_readLength+1.0), _mu.get(thisTranscriptID)+0.0);
				//out_exprs.printf(_gene.getID()+"\t"+thisTranscriptID+"\t"+_hasFlatPrior+"\t%e\t"+flag+"\t"+count+"\t%f\t%e\t%e\n", _priors[k]+0.0, _readCounts[k]+0.0, _readCounts[k]/(_transcriptLengths[k]-_readLength+1.0), _mu[k]+0.0);

				//System.out.println("DONE: "+_gene.getID()+"\t"+_transcriptIDs.get(k)+"\tflag="+flag+"\tk="+count+"\tmu="+_mu[k]+"\tcount="+_readCounts[k]+"\texp="+_readCounts[k]/(_transcriptLengths[k]-_readLength+1.0));
			}*/
		}

		//return(sb.toString());
	}

	public HashMap<String, String> sampleReadAssignments(HashMap<String, String> readMap){
		Iterator<String> reads = _gene.getReads().keySet().iterator();
		String thisRead;
		//ArrayList<String> theseTranscripts;
		//HashMap<String, String> readMap = new HashMap<String, String>();

		while(reads.hasNext()){
			thisRead = reads.next();
			//theseTranscripts = _reads.get(thisRead);
			//double[] transcriptFracs = new double[_reads.get(thisRead).size()];
			HashMap<String, BigDecimal> transcriptFracs = new HashMap<String, BigDecimal>();
			double transcriptSum = 0.0;

			if(_gene.getTranscriptIDs().size() == 1){
				transcriptFracs.put(_gene.getTranscriptIDs().get(0), new BigDecimal(1.0));
				transcriptSum = 1.0;
			}else{
				for(int k=0;k<_gene.getTranscriptIDs().size();k++){
					String thisTranscriptID = _gene.getTranscriptIDs().get(k);
					
					//if(theseTranscripts.contains(thisTranscriptID))
					if(_gene.getReads().get(thisRead).contains(thisTranscriptID)){
						transcriptFracs.put(thisTranscriptID, _readCounts.get(thisTranscriptID));
						transcriptSum += _readCounts.get(thisTranscriptID).doubleValue()+0.0;
					}
					else{
						transcriptFracs.put(thisTranscriptID, new BigDecimal(0.0));
					}
				}
			}

			for(int k=0;k<_gene.getTranscriptIDs().size();k++){
				String thisTranscriptID = _gene.getTranscriptIDs().get(k);
				transcriptFracs.put(thisTranscriptID, new BigDecimal(transcriptFracs.get(thisTranscriptID).doubleValue() / (transcriptSum+0.0)));
			}

			double rand = Math.random();
			double cumsum = 0.0;
			for(int k=0;k<_gene.getTranscriptIDs().size();k++){
				String thisTranscriptID = _gene.getTranscriptIDs().get(k);
				cumsum += transcriptFracs.get(thisTranscriptID).doubleValue();
				if(cumsum >= rand){
					readMap.put(thisRead, thisTranscriptID);
					break;
				}
			}
		}
		return(readMap);
	}


	private int _readLength = 29;
	public void setReadLength(int readLength){ this._readLength = readLength; }

	private HashMap<String,BigDecimal> mu2tau(HashMap<String,BigDecimal> priors){
		HashMap<String,BigDecimal> tau = new HashMap<String,BigDecimal>();

		if(priors.size() != _gene.getLengths().size()){
			System.err.println("mu2tau: input vectors are not of the same length\n");
			System.exit(1);
		}

		double tauSum = 0.0;
		for(int k=0;k<_gene.getTranscriptIDs().size();k++){
			String thisTranscriptID = _gene.getTranscriptIDs().get(k);
			tau.put(thisTranscriptID, new BigDecimal(priors.get(thisTranscriptID).doubleValue() * Math.max(0, _gene.getLengths().get(thisTranscriptID)-_readLength)));
			tauSum += tau.get(thisTranscriptID).doubleValue();
		}
		/*for(int k=0;k<priors..length;k++){
			tau[k] = priors[k] * Math.max(0, _transcriptLengths[k]-_readLength);
			tauSum += tau[k]; 
		}*/

		if(Math.abs(tauSum - _TINY) < 0){
			System.err.println("mu_2_tau: sum of prior is too small, please check parameters\n");
			System.exit(1);
		}

		for(int k=0;k<_gene.getTranscriptIDs().size();k++){
			String thisTranscriptID = _gene.getTranscriptIDs().get(k);
			tau.put(thisTranscriptID, new BigDecimal(tau.get(thisTranscriptID).doubleValue()/tauSum));
		}
		/*for(int k=0;k<priors.length;k++){
			tau[k] = tau[k]/tauSum;
		}*/

		return(tau);
	}


	private void allocateRead(String readID){//, double[] tau){//, double[] readCounts){
		//HashMap<String, Double> addVal = _tau;
		HashMap<String, BigDecimal> addVal = new HashMap<String, BigDecimal>();
		double sumVal = 0.0;
		for(int k=0;k<_gene.getTranscriptIDs().size();k++){
			String thisTranscriptID = _gene.getTranscriptIDs().get(k);
			if(_gene.getReads().get(readID).contains(thisTranscriptID)){
				addVal.put(thisTranscriptID, _tau.get(thisTranscriptID));
			}else{
				addVal.put(thisTranscriptID, new BigDecimal(0.0));
			}
			sumVal += addVal.get(thisTranscriptID).doubleValue();
		}

		if(sumVal > 0.0){
			//*^^^ split the read proportionally to the transcripts^^^*/
			for(int k=0;k<_gene.getTranscriptIDs().size();k++){
				String thisTranscriptID = _gene.getTranscriptIDs().get(k);
				//readCounts[k] = readCounts[k] + addVal[k]/sumVal;
				//readCounts[k] = readCounts[k] + ((addVal[k]/sumVal) * _gene.getReadCount(readID));
				_readCounts.put(thisTranscriptID, new BigDecimal(_readCounts.get(thisTranscriptID).doubleValue() + ((addVal.get(thisTranscriptID).doubleValue()/sumVal) * _gene.getReadCount(readID))));
			}
		}

		//return(readCounts);
	}


	

	/**
	 * 
	 */
	public void expectationStep(){

		_tau = mu2tau(_mu);
		
		//_readCounts = _tau;
		_readCounts = new HashMap<String, BigDecimal>();
		Iterator<String> it = _tau.keySet().iterator();
		while(it.hasNext()){
			String id = it.next();
			_readCounts.put(id, new BigDecimal(0.0));
		}

		it = _gene.getReads().keySet().iterator();
		while(it.hasNext())
			allocateRead(it.next());//, _tau, _readCounts);
		
		//for(int k=0; k<_reads.size(); k++){
			//_readCounts = allocateRead((String)_reads.get.keySet()_readIDs[k], _tau, _readCounts);
		//}
	}

	/**
	 * 
	 * @return
	 */
	public HashMap<String,BigDecimal> maximisationStep(){
		//HashMap<String,Double> mu_new = _mu;
		HashMap<String,BigDecimal> mu_new = new HashMap<String,BigDecimal>();
		double s = 0.0;

		for(int k=0;k<_gene.getTranscriptIDs().size();k++){
			String thisTranscriptID = _gene.getTranscriptIDs().get(k);
			mu_new.put(thisTranscriptID, new BigDecimal(_readCounts.get(thisTranscriptID).doubleValue() / Math.max(1, _gene.getLengths().get(thisTranscriptID).doubleValue()-_readLength)));
			s += mu_new.get(thisTranscriptID).doubleValue();
		}
		
		/*for(int k=0;k<_readCounts.length;k++){
			mu_new[k] = _readCounts[k]/Math.max(1, _transcriptLengths[k]-_readLength);
			s += mu_new[k];
		}*/
		
		
		if(s > 0.0){
			for(int k=0;k<_gene.getTranscriptIDs().size();k++){
				String thisTranscriptID = _gene.getTranscriptIDs().get(k);
				mu_new.put(thisTranscriptID, new BigDecimal(mu_new.get(thisTranscriptID).doubleValue()/s));
			}
			/*for(int k=0;k<_readCounts.length;k++){
				mu_new[k] = mu_new[k]/s;
			}*/
		}
		return(mu_new);
	}

	
	/**
	 * 
	 * @param mu_old
	 * @param mu_new
	 * @return
	 */
	public double calcConvergence(HashMap<String,BigDecimal> mu_old, HashMap<String,BigDecimal> mu_new){
		if(mu_old.size() != mu_new.size()){
			System.err.println("cover_distance: two transcripts expression percent not of the same length\n");
			System.exit(1);
		}
		double dis = 0.0;
		for(int k=0;k<_gene.getTranscriptIDs().size();k++){
			String thisTranscriptID = _gene.getTranscriptIDs().get(k);
			dis += Math.abs(mu_old.get(thisTranscriptID).doubleValue() - mu_new.get(thisTranscriptID).doubleValue()); 
		}
		/*for(int k=0;k<mu_old.length;k++){
			dis += Math.abs(mu_old[k]-mu_new[k]);
		}*/
		return(dis);
	}
}
