package footprintAlignments;

import net.sf.samtools.*;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
import static java.lang.Math.*;
import java.math.*;
import java.util.ArrayList;
import java.util.Map;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.Set;
import java.util.logging.Level;
import java.util.logging.Logger;


/**
 *
 * @author J. David McPeek
 */
public class FootprintAlignmentSummary {
    
    File inputSortedBamFile = null;
    File outputSummaryFile = null;
    HashMap<String, String> transcriptTypeMap = null;
    // TODO: perhaps have some other interfaces if clients already have the gencode data
    
    public FootprintAlignmentSummary(File inputSortedBamFile, File outputSummaryFile, File gencodeDataFile) throws IOException {
        this.inputSortedBamFile = inputSortedBamFile;
        this.outputSummaryFile = outputSummaryFile;
        this.transcriptTypeMap = new gtfParser(gencodeDataFile).getTypeMap(); 
    }
    
    protected static class TranscriptSummary {
        String transcript;
        String type;
        int numberOfReads;
        int maxInFrame;  // k
        int midInFrame;  // m
        int minInFrame;  // m
        
        double alignmentPValue;
        
        private static HashMap<Integer, BigInteger> factorialCache = new HashMap<>();  // don't redo computations
        
        TranscriptSummary(String transcriptID, String transcriptType, int totalReads, int maxAligned, int midAligned, int minAligned) {
            transcript = transcriptID;
            type = transcriptType;
            numberOfReads = totalReads;
            maxInFrame = maxAligned;
            midInFrame = midAligned;
            minInFrame = minAligned;
            alignmentPValue = this.getAlignmentPValue();
        }
        
        double getAlignmentPValue() {
            // find the absolute number of aligned reads by multiplying the total numberOfReads with the percentage that are in frame
            BigDecimal sum = new BigDecimal("0");  // in the end, we subtract `sum` from 1 to get a probability
            BigDecimal factorialN = new BigDecimal(factorial(numberOfReads));  // factorial(n) converted to BigDecimal
            BigInteger threePowerN = BigInteger.valueOf(3).pow(numberOfReads);
            BigInteger denominator;  
            // TODO : HASH PREVIOUSLY COMPUTED VALUES OF FACTORIAL 

//            HashMap<Integer, BigInteger> factorialCache = new HashMap<>();

            BigInteger iprime;
            BigInteger jprime;
            BigInteger mprime;


    //        int m;
            for (int i = 0; i < maxInFrame; i++) {
                for (int j = 0; j < maxInFrame; j++) {
                    if (j + i > numberOfReads) break; 
                    for (int m = 0; m < maxInFrame; m++) {

                        // the reason why this optimization doesn't work is because it allows m to be 0 multiple times within what would be the same loop. 
                        // m should only equal 0 once.
                        // the way this problem could arise is if you have k = 7. Both (i = 5, j = 2) and (j = 5, i = 2) make = 0; same goes for (i = 4, j = 3) and vice-versa. There are tons 
                        // of combinations that allow m to be zero multiple times. 
                        // TODO: CHECK TO SEE IF THIS PROBLEM MAKES SENSE. I THINK THE PROBLEM IS THAT IT'S AN OVERESTIMATION. IF THAT'S THE CASE, THEN FIGURE OUT WAYS TO SWITCH-ON SWITCH-OFF m.
    //                    if (i + j <= n) {
    //                        m = n - (i + j);
    //                    } else {
    //                        break;
    //                    }

    //                   
                        if (i + j + m == numberOfReads) { 
                            // multiply denominator a number of times
                            if (factorialCache.containsKey(i)) {
                                iprime = factorialCache.get(i);
                            } else {
                                iprime = factorial(i);
                                factorialCache.put(i, iprime);
                            }

                            if (factorialCache.containsKey(j)) {
                                jprime = factorialCache.get(j);
                            } else {
                                jprime = factorial(j);
                                factorialCache.put(j, jprime);
                            }

                            if (factorialCache.containsKey(m)) {
                                mprime = factorialCache.get(m);
                            } else {
                                mprime = factorial(m);
                                factorialCache.put(m, mprime);
                            }

                            denominator = iprime.multiply(jprime);
                            denominator = denominator.multiply(mprime);
                            denominator = denominator.multiply(threePowerN);
                            BigDecimal denominatorDecimal = new BigDecimal(denominator);
                            denominatorDecimal = factorialN.divide(denominatorDecimal, 200, RoundingMode.HALF_UP);
                            sum = sum.add(denominatorDecimal);
                        } else if (i + j + m > numberOfReads) break;
                    }
                }
            }

    //        System.out.println(new BigDecimal("1").subtract(sum));

            return new BigDecimal("1").subtract(sum).doubleValue();
            
        }
        
        double getAlignmentPValueHelper(int maxAligned, int midAligned, int minAligned) {
            int total = maxAligned + midAligned + minAligned;
            BigDecimal likelihood = new BigDecimal("1");
            likelihood = likelihood.multiply(new BigDecimal(pow(1/3d, total) + "")); 
            BigInteger totalPrime  = factorial(total);
            BigInteger maxPrime = factorial(maxAligned);
            BigInteger midPrime = factorial(midAligned);
            BigInteger minPrime = factorial(minAligned);
            BigInteger denominator = maxPrime.multiply(midPrime.multiply(minPrime)); 
            BigInteger multinomialFactor =  totalPrime.divide(denominator);
            likelihood = likelihood.multiply(new BigDecimal(multinomialFactor));
            
            double d;  // d = 3! / (3 - numberOfDifferentlySizedBuckets)!
            if (maxAligned == midAligned && midAligned == minAligned) d = 1;
            else if (maxAligned != midAligned && midAligned == minAligned) d = 3;
            else if (maxAligned == midAligned && midAligned != minAligned) d = 3;
            else d = 6;
            
            likelihood = likelihood.multiply(new BigDecimal(d + "")); 

            return likelihood.doubleValue();
        }
        
       public void writeSummaryToFile(BufferedWriter outputFile) throws IOException {
            outputFile.write(transcript + "\t");
            outputFile.write(type + "\t");
            outputFile.write(numberOfReads + "\t");
            outputFile.write(String.valueOf(maxInFrame + "\t"));
            outputFile.write(String.valueOf(midInFrame + "\t"));
            outputFile.write(String.valueOf(minInFrame + "\t"));
            outputFile.write(String.valueOf(alignmentPValue));
            outputFile.newLine();
       }
        
        
    }
    
    public void createSummary() throws IOException {
        
        try (SAMFileReader inputBam = new SAMFileReader(inputSortedBamFile);
                BufferedWriter outputFile = new BufferedWriter(new FileWriter(this.outputSummaryFile))) {
            
            // large array of individual transcript summaries that we're eventually going to write to a file
            ArrayList<TranscriptSummary> transcriptSummaries = new ArrayList<>();
            
            Map<String, ArrayList<SAMRecord>> transcriptMap = new HashMap<>();
            
            for (final SAMRecord samRecord : inputBam) {
                // we're only looking at reads that match 29 base pairs
                if (samRecord.getCigar().getReferenceLength() != 29) {
                    continue;
                }
                String transcriptID = samRecord.getReferenceName();
                
                if (transcriptMap.containsKey(transcriptID)) {
                    transcriptMap.get(transcriptID).add(samRecord);
                } else {
                    ArrayList<SAMRecord> sameTranscriptReads = new ArrayList<>();
                    sameTranscriptReads.add(samRecord);
                    transcriptMap.put(transcriptID, sameTranscriptReads);
                }
                
            }
            
            Iterator<HashMap.Entry<String, ArrayList<SAMRecord>>> it = transcriptMap.entrySet().iterator();
            
            String fileHeader = "\"transcript.id\"\t\"transcript.type\"\t\"total.reads\"\t\"aligned.max\"\t\"alignment.mid\"\t\"alignment.min\"\t\"alignment.score\"";
            outputFile.write(fileHeader);
            outputFile.newLine();
            
            while (it.hasNext()) {
                
                Entry<String, ArrayList<SAMRecord>> thisTranscript = it.next();
                
                String transcriptID = thisTranscript.getKey();
                ArrayList<SAMRecord> reads = thisTranscript.getValue();
                int totalReads = reads.size();
                
                HashMap<Integer, Integer> readFrames = new HashMap<>();
                readFrames.put(0, 0);
                readFrames.put(1, 0);
                readFrames.put(2, 0);
                int readPosition;
                for (int i = 0; i < totalReads; i++) {
                    readPosition = reads.get(i).getAlignmentStart() % 3;
                    readFrames.put(readPosition, readFrames.get(readPosition) + 1);
                }
                int bin0 = readFrames.get(0);
                int bin1 = readFrames.get(1);
                int bin2 = readFrames.get(2);
                int maxAligned = Math.max(Math.max(bin0, bin1), bin2);
                int minAligned = Math.min(Math.min(bin0, bin1), bin2);
                int midAligned = (bin0 + bin1 + bin2) - maxAligned - minAligned;
                
                String transcriptType = transcriptTypeMap.get(transcriptID);
                TranscriptSummary thisSummary = new TranscriptSummary(transcriptID, transcriptType, totalReads, maxAligned, midAligned, minAligned);
                transcriptSummaries.add(thisSummary);
                
                thisSummary.writeSummaryToFile(outputFile);
                
            }
            
        } catch (IOException ex) {
            Logger.getLogger(FootprintAlignmentSummary.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            System.out.println("You finished, yay!!!!");
        }
        
        
    }
    
    /**
     *
     * @param stringLength the length of the original string we're taking combinations of
     * @param k the size of each combinations group
     * @return the BigInteger representation of the number of combinations described by stringLength-choose-k
     */
    public static BigInteger numberOfCombinations(int stringLength, int k) {
        
        try {
            if (k < 1 || k > stringLength) {
                throw new IllegalArgumentException();
            }
        } catch (IllegalArgumentException e) {
            System.out.println("Illegal argument exception: " + e.getMessage());
            System.out.println("You cannot find the combinations of k symbols when k is above the given string length or below 1");
            System.exit(1);
        }
        
        // numerator factorial
        BigInteger numPermutations = new BigInteger("1");
        
        // implement classic permutations with `i` number of slots left over
        for (int i = stringLength; i > stringLength - k; i--) {
            numPermutations = numPermutations.multiply(new BigInteger(i + ""));
        }
        BigInteger kPrime = new BigInteger("1");
        for (int i = k; i > 0; i--) {
            kPrime = kPrime.multiply(new BigInteger(i + "")); 
        }
        
        return numPermutations.divide(kPrime);
        
    }
    
    public static BigInteger factorial(int n) {
        BigInteger result = BigInteger.ONE;
        while (n != 0) {
            result = result.multiply(BigInteger.valueOf(n));
            n--;
        }
        return result;
    } 
    

        static double helper(int maxAligned, int midAligned, int minAligned) {
            int total = maxAligned + midAligned + minAligned;
            BigDecimal likelihood = new BigDecimal("1");
            likelihood = likelihood.multiply(new BigDecimal(pow(1/3d, total) + "")); 
            BigInteger totalPrime  = factorial(total);
            BigInteger maxPrime = factorial(maxAligned);
            BigInteger midPrime = factorial(midAligned);
            BigInteger minPrime = factorial(minAligned);
            BigInteger denominator = maxPrime.multiply(midPrime.multiply(minPrime)); 
            BigInteger multinomialFactor =  totalPrime.divide(denominator);
            likelihood = likelihood.multiply(new BigDecimal(multinomialFactor));
            
            double d;  // d = 3! / (3 - numberOfDifferentlySizedBuckets)!
            if (maxAligned == midAligned && midAligned == minAligned) d = 1;
            else if (maxAligned != midAligned && midAligned == minAligned) d = 3;
            else if (maxAligned == midAligned && midAligned != minAligned) d = 3;
            else d = 6;
            
            likelihood = likelihood.multiply(new BigDecimal(d + "")); 

            return likelihood.doubleValue();
        }
    
    
    


    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {

//        File inputBamFile = new File("/Users/jdmcpeek/FootprintAlignments/data/sorted.bam");
//        File outputFile = new File("./data/outTest.txt");
//        File gencodeData = new File("./data/gencode.v21.transcripts.gtf.txt"); 
//        
//        FootprintAlignmentSummary summary = new FootprintAlignmentSummary(inputBamFile, outputFile, gencodeData);
//        summary.createSummary();
//        
//        System.out.print("finished.");
        
         
        
        int num = 400;
//          
          double f = helper(num, 0, 0);
          System.out.println(f);
//          

        System.out.println(belowK2(5, 2));
        
//        double sum = helper(2, 1, 1) + helper(2, 0, 2) + helper(2, 2, 0) + helper(3,0,1) + helper(3, 1,0) + helper(4, 0, 0);
//        double sum2 = helper(2, 1, 1) + helper(2, 2, 0) + helper(3,0,1) + helper(4, 0, 0);
//        
//        System.out.println(sum2);
//        
        double sum = helper(3, 1, 1) + helper(3, 0, 2) + helper(4,0,1) + helper(5, 0, 0) + helper(2, 2, 1);
        System.out.println(sum);
        
        
        System.out.println(belowK2(8, 4));
        
        sum = helper(4, 2, 2) + helper(4, 1, 3) + helper(4, 0, 4)/2 + helper(5, 1, 2) + helper(5, 0, 3) + helper(6, 1, 1) + helper(6, 0, 2) + helper(7, 0, 1) + helper(8,0, 0); 
        System.out.println(sum);
        
        
        System.out.println(belowK2(3, 1)); 
        sum = helper(1, 1, 1) + helper(2, 1, 0) + helper(3, 0 ,0); 
        System.out.println(sum);
        
        
        
        System.out.println("========================");
        
        // this method works for cases when k > 0.4N
        
        int n = 10;
        int k = 5;
        
        System.out.println(belowK2(n, k)); 
        // find the probability of getting k similar alignments OR BETTER.
        int bigBucket = n;
        sum = 0;
        
        while (bigBucket >= k) {
            // what's left for the other buckets
            int left = n - bigBucket;
            int middleBucket = left;
            int smallBucket = 0;
            
//            if (left > k) {
//                break;
//            }
            
            if (left < 2) {
                sum += helper(bigBucket, left, 0);
            } else if (left % 2 == 0) {  // even cases
                sum += helper(bigBucket, left/2, left/2);  // the last case. We're doing it here to avoid a weird while loop
                while (middleBucket != smallBucket) {
                    double term = helper(bigBucket, middleBucket, smallBucket);
                    if (bigBucket == middleBucket) {
                        sum += term/2;  // avoid duplicates
                    } else {
                        sum += term;
                    }
                    middleBucket--;
                    smallBucket++; 
                }
            } else if (left % 2 == 1) {  // odd cases
                sum += helper(bigBucket, 1 + left/2, left/2);
                while (middleBucket - 1 != smallBucket) {
                    double term = helper(bigBucket, middleBucket, smallBucket);
                    if (bigBucket == middleBucket) {
                        sum += term/2;  // avoid duplicates
                    } else {
                        sum += term;
                    }
                    middleBucket--;
                    smallBucket++; 
                }
            }
            
            bigBucket--;
        }
        
        System.out.println(sum);
        
        
        
        
//        double test = helper(6, 0, 0) + helper(5, 1, 0) + helper(4,2,0) + helper(4,1,1) + helper(3, 3, 0) +helper(3,2,1);
//        System.out.println(test);
//          
        
        
        
    
    }
    
    static double belowK2(int n, int k) {
            // find the absolute number of aligned reads by multiplying the total numberOfReads with the percentage that are in frame
            BigDecimal sum = new BigDecimal("0");  // in the end, we subtract `sum` from 1 to get a probability
            BigDecimal factorialN = new BigDecimal(factorial(n));  // factorial(n) converted to BigDecimal
            BigInteger threePowerN = BigInteger.valueOf(3).pow(n);
            BigInteger denominator;  
            // TODO : HASH PREVIOUSLY COMPUTED VALUES OF FACTORIAL 

            HashMap<Integer, BigInteger> factorialCache = new HashMap<>();

            BigInteger iprime;
            BigInteger jprime;
            BigInteger mprime;


    //        int m;
            for (int i = 0; i < k; i++) {
                for (int j = 0; j < k; j++) {
                    if (j + i > n) break; 
                    for (int m = 0; m < k; m++) {

                        // the reason why this optimization doesn't work is because it allows m to be 0 multiple times within what would be the same loop. 
                        // m should only equal 0 once.
                        // the way this problem could arise is if you have k = 7. Both (i = 5, j = 2) and (j = 5, i = 2) make = 0; same goes for (i = 4, j = 3) and vice-versa. There are tons 
                        // of combinations that allow m to be zero multiple times. 
                        // TODO: CHECK TO SEE IF THIS PROBLEM MAKES SENSE. I THINK THE PROBLEM IS THAT IT'S AN OVERESTIMATION. IF THAT'S THE CASE, THEN FIGURE OUT WAYS TO SWITCH-ON SWITCH-OFF m.
    //                    if (i + j <= n) {
    //                        m = n - (i + j);
    //                    } else {
    //                        break;
    //                    }

    //                   
                        if (i + j + m == n) { 
                            // multiply denominator a number of times
                            if (factorialCache.containsKey(i)) {
                                iprime = factorialCache.get(i);
                            } else {
                                iprime = factorial(i);
                                factorialCache.put(i, iprime);
                            }

                            if (factorialCache.containsKey(j)) {
                                jprime = factorialCache.get(j);
                            } else {
                                jprime = factorial(j);
                                factorialCache.put(j, jprime);
                            }

                            if (factorialCache.containsKey(m)) {
                                mprime = factorialCache.get(m);
                            } else {
                                mprime = factorial(m);
                                factorialCache.put(m, mprime);
                            }

                            denominator = iprime.multiply(jprime);
                            denominator = denominator.multiply(mprime);
                            denominator = denominator.multiply(threePowerN);
                            BigDecimal denominatorDecimal = new BigDecimal(denominator);
                            denominatorDecimal = factorialN.divide(denominatorDecimal, 200, RoundingMode.HALF_UP);
                            sum = sum.add(denominatorDecimal);
                        } else if (i + j + m > n) break;
                    }
                }
            }

    //        System.out.println(new BigDecimal("1").subtract(sum));

            return new BigDecimal("1").subtract(sum).doubleValue();
            
        }
}
