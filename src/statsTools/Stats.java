package statsTools;

import java.util.ArrayList;
import java.util.Collections;

import org.h2.util.New;

import utils.Object_utils;


public class Stats {

	/**
	 * Calculates the Kulbach-Liebler divergence for the given barcode counts
	 * @param numberOfObservationsOfALLBarcodeSequences
	 * @param numberOfObservationsOfThisBarcodeSequence
	 * @param nBarcodeBases
	 * @param nObservedBarcodes
	 * @return
	 */
	public static double calculateKL(double probability, int numberOfObservationsOfThisBarcodeSequence, int nBarcodeBases, int nObservedBarcodes){

		double frequency = (double) numberOfObservationsOfThisBarcodeSequence * probability;
		return frequency * (Math.log(frequency*nObservedBarcodes) / Math.log(2));
	}


	/**
	 * Calculates the Kullback-Leibler divergence of a set of observed probabilities vs a non-uniform null distribution
	 * @param observed
	 * @param expected
	 * @return
	 */
	public static double calculateKL(double[] observed, double[] expected){
		if(observed.length > 0  &&  observed.length == expected.length){
			double kl = 0.0;
			for(int i=0;i<observed.length;i++){
				if(observed[i] > 0.0)
					kl += observed[i]*log2(observed[i]/expected[i]);
			}
			return kl;
		}else{
			return 99.0;
		}
	}


	/**
	 * Calculates the Kullback-Leibler divergence of a set of observed probabilities vs a uniform null distribution where the probability of the latter is always 'expected'
	 * @param observed
	 * @param expected
	 * @return
	 */
	public static double calculateKL(double[] observed, double expected){
		if(observed.length > 0){
			double kl = 0.0;
			for(int i=0;i<observed.length;i++){
				if(observed[i] > 0.0)
					kl += (observed[i]+0.0)*log2((observed[i]+0.0)/(expected+0.0));
			}
			return kl;
		}else{
			return 99.0;
		}
	}

	public static double calculateEntropy(double[] observed){
		if(observed.length > 0){
			double h = 0.0;
			for(int i=0;i<observed.length;i++){
				if(observed[i] > 0.0)
					h += -1*((observed[i]+0.0)*log2(observed[i]+0.0));
			}
			return h;
		}else{
			return 99.0;
		}
	}


	/**
	 * Calculates Kolmogorov-Smirnov (KS) distance against an expected null distribution of the same length as the observed distribution
	 * @param observed
	 * @return
	 */
	public static double calculateKS(double[] observed, double[] expected){
		double result = 0.0;
		if(observed.length > 0  &&  observed.length == expected.length){
			//double[] observedCDF = Object_utils.initdoubleArray(observed.length);
			//double[] expectedCDF = Object_utils.initdoubleArray(expected.length);
			double sum_observed = 0.0;
			double sum_expected = 0.0;
			for(int i=0;i<observed.length;i++){
				sum_observed += observed[i];
				sum_expected += expected[i];
				double diff = Math.sqrt(Math.pow(sum_observed-sum_expected, 2));
				if(diff > result)
					result = diff;
			}
		}
		return result;
	}

	/**
	 * Calculates Kolmogorov-Smirnov (KS) distance against a uniform distribution of the same length as the observed distribution
	 * @param observed
	 * @return
	 */
	public static double calculateKS(double[] observed){
		int n = observed.length;
		double[] uniformDist = Object_utils.initArray_double(n);
		for(int i=0;i<n;i++)
			uniformDist[i] = 1.0/(n+0.0);
		return calculateKS(observed, uniformDist);
	}
	
	
	
	/*public static double median(ArrayList<Double> values){
		Collections.sort(values);
		return values.get((int)Math.round(values.size()/2.0));
	}*/
	public static int median(ArrayList<Integer> values){
		Collections.sort(values);
		return values.get((int)Math.round(values.size()/2.0));
	}

	public static double log2(double val){
		//Logb x = Loga x/Loga b
		return Math.log10(val) / Math.log10(2.0);
	}


	public static void main(String[] args) {
		// TODO Auto-generated method stub

	}

}
