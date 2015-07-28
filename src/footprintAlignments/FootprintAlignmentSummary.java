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
    private static final HashMap<Integer, BigInteger> factorialCache = new HashMap<Integer, BigInteger>();  // don't redo computations
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
        int maxInFrame;  
        int midInFrame;  
        int minInFrame;  
        
        double alignmentPValue;
        double adjustedScore;
        
        
        
        TranscriptSummary(String transcriptID, String transcriptType, int totalReads, int maxAligned, int midAligned, int minAligned, int numberOfTests) {
            transcript = transcriptID;
            type = transcriptType;
            numberOfReads = totalReads;
            maxInFrame = maxAligned;
            midInFrame = midAligned;
            minInFrame = minAligned;
            alignmentPValue = this.getAlignmentPValue();
            adjustedScore = alignmentPValue * numberOfTests;
        }
        
        double getAlignmentPValue() {
            // find the absolute number of aligned reads by multiplying the total numberOfReads with the percentage that are in frame
            BigInteger sum = new BigInteger("0");  // in the end, we subtract `sum` from 1 to get a probability

            for (int i = 0; i < maxInFrame; i++) {
                for (int j = 0; j < maxInFrame; j++) {
                    for (int m = 0; m < maxInFrame; m++) {
                        if (i + j + m == numberOfReads) { 
//                            denominatorDecimal = factorialN.divide(denominatorDecimal, 200, RoundingMode.HALF_UP);
                            sum = sum.add(binPatterns(numberOfReads, i, j, m));
                        } 
                    }
                }
            }

            BigInteger totalPossibilities = BigInteger.valueOf(3).pow(numberOfReads);
            BigInteger matches = totalPossibilities.subtract(sum);
            
            return new BigDecimal(matches).divide(new BigDecimal(totalPossibilities), 200, RoundingMode.HALF_UP).doubleValue(); 
            
        }
        
        private BigInteger binPatterns(int n, int j, int m, int t) {
            return getFactorial(n).divide(getFactorial(j).multiply(getFactorial(m)).multiply(getFactorial(t)));
        }
        
        private BigInteger getFactorial(int n) {
            if (factorialCache.containsKey(n)) {
                return factorialCache.get(n);
            } else {
                BigInteger factorialN = factorial(n);
                factorialCache.put(n, factorialN);
                return factorialN;
            }
        }
        
        
        
//        double getAlignmentPValueHelper(int maxAligned, int midAligned, int minAligned) {
//            int total = maxAligned + midAligned + minAligned;
//            BigDecimal likelihood = new BigDecimal("1");
//            likelihood = likelihood.multiply(new BigDecimal(pow(1/3d, total) + "")); 
//            BigInteger totalPrime  = factorial(total);
//            BigInteger maxPrime = factorial(maxAligned);
//            BigInteger midPrime = factorial(midAligned);
//            BigInteger minPrime = factorial(minAligned);
//            BigInteger denominator = maxPrime.multiply(midPrime.multiply(minPrime)); 
//            BigInteger multinomialFactor =  totalPrime.divide(denominator);
//            likelihood = likelihood.multiply(new BigDecimal(multinomialFactor));
//            
//            double d;  // d = 3! / (3 - numberOfDifferentlySizedBuckets)!
//            if (maxAligned == midAligned && midAligned == minAligned) d = 1;
//            else if (maxAligned != midAligned && midAligned == minAligned) d = 3;
//            else if (maxAligned == midAligned && midAligned != minAligned) d = 3;
//            else d = 6;
//            
//            likelihood = likelihood.multiply(new BigDecimal(d + "")); 
//
//            return likelihood.doubleValue();
//        }
        
       
        

       public void writeSummaryToFile(BufferedWriter outputFile) throws IOException {
            outputFile.write(transcript + "\t");
            outputFile.write(type + "\t");
            outputFile.write(numberOfReads + "\t");
            outputFile.write(String.valueOf(maxInFrame + "\t"));
            outputFile.write(String.valueOf(midInFrame + "\t"));
            outputFile.write(String.valueOf(minInFrame + "\t"));
            outputFile.write(String.valueOf(alignmentPValue + "\t"));
            outputFile.write(String.valueOf(adjustedScore)); 
            outputFile.newLine();
       }
        
        
    }
    
    public void createSummary() throws IOException {
        
        SAMFileReader inputBam = null;
        BufferedWriter outputFile = null;
        
        try {
            
            inputBam = new SAMFileReader(inputSortedBamFile);
            outputFile = new BufferedWriter(new FileWriter(this.outputSummaryFile));
            
            // large array of individual transcript summaries that we're eventually going to write to a file
            ArrayList<TranscriptSummary> transcriptSummaries = new ArrayList<TranscriptSummary>();
            
            Map<String, ArrayList<SAMRecord>> transcriptMap = new HashMap<String, ArrayList<SAMRecord>>();
            
            for (final SAMRecord samRecord : inputBam) {
                // we're only looking at reads that match 29 base pairs
                if (samRecord.getCigar().getReferenceLength() != 29) {
                    continue;
                }
                String transcriptID = samRecord.getReferenceName();
                
                if (transcriptMap.containsKey(transcriptID)) {
                    transcriptMap.get(transcriptID).add(samRecord);
                } else {
                    ArrayList<SAMRecord> sameTranscriptReads = new ArrayList<SAMRecord>();
                    sameTranscriptReads.add(samRecord);
                    transcriptMap.put(transcriptID, sameTranscriptReads);
                }
                
            }
            
            Iterator<HashMap.Entry<String, ArrayList<SAMRecord>>> it = transcriptMap.entrySet().iterator();
            
            int numberOfTranscripts = transcriptMap.size();
            
            String fileHeader = "\"transcript.id\"\t\"transcript.type\"\t\"total.reads\"\t\"aligned.max\"\t\"alignment.mid\"\t\"alignment.min\"\t\"alignment.score\"\t\"adjusted.score\"";
            outputFile.write(fileHeader);
            outputFile.newLine();
            
            while (it.hasNext()) {
                
                Entry<String, ArrayList<SAMRecord>> thisTranscript = it.next();
                
                String transcriptID = thisTranscript.getKey();
                ArrayList<SAMRecord> reads = thisTranscript.getValue();
                int totalReads = reads.size();
                
                HashMap<Integer, Integer> readFrames = new HashMap<Integer, Integer>();
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
                
                
                
                TranscriptSummary thisSummary = new TranscriptSummary(transcriptID, transcriptType, totalReads, maxAligned, midAligned, minAligned, numberOfTranscripts);
                transcriptSummaries.add(thisSummary);
                
                thisSummary.writeSummaryToFile(outputFile);
                
            }
            
        } catch (IOException ex) {
            Logger.getLogger(FootprintAlignmentSummary.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            if (inputBam != null) inputBam.close();
            if (outputFile != null) outputFile.close();
               
            System.out.println("You finished, yay!!!!");
        }
        
        
    }
    

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
    
    public static void writeMapToFile(Map map, File file) throws IOException {
        BufferedWriter outputFile = null;
        try {
            outputFile = new BufferedWriter(new FileWriter(file));
            Iterator it = map.entrySet().iterator();
            outputFile.write("\"transcript.id\"\t\"transcript.type\"");
            outputFile.newLine();
            while (it.hasNext()) {
                Map.Entry pair = (Map.Entry)it.next();
                outputFile.write(pair.getKey() + "\t" + pair.getValue());
                outputFile.newLine();
                
            }
        } catch (IOException ex) {
            Logger.getLogger(FootprintAlignmentSummary.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            if (outputFile != null) outputFile.close();
        }
                

    }

 
    
    // create a function that computes the significance of a set of alignments.
    // returns true/false if it falls below the critical value, the necessary level of significance.
    // public boolean chiSquaredTest();

    

    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {

        File inputBamFile = new File("/Users/jdmcpeek/FootprintAlignments/data/sorted.bam");
        File outputFile = new File("./data/outTest3.txt");
        File gencodeData = new File("./data/gencode.v21.transcripts.gtf.txt"); 
//        
        FootprintAlignmentSummary summary = new FootprintAlignmentSummary(inputBamFile, outputFile, gencodeData);
        summary.createSummary();
//        
        
        
        
//        writeMapToFile(summary.transcriptTypeMap, new File("./data/transcriptTypeMap.txt"));
     

    
    }
    

}
