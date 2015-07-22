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
            double likelihood = getAlignmentPValueHelper(maxInFrame, midInFrame, minInFrame);
            
//            while (aligned < numberOfReads) {
//                aligned++;
//                likelihood += getAlignmentPValueHelper(aligned, (numberOfReads - aligned));
//            }
            return likelihood;
            
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
            result = result.multiply(new BigInteger(n + ""));
            n--;
        }
        return result;
    } 
    

    
    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {

        File inputBamFile = new File("/Users/jdmcpeek/FootprintAlignments/data/sorted.bam");
        File outputFile = new File("./data/outTest.txt");
        File gencodeData = new File("./data/gencode.v21.transcripts.gtf.txt"); 
        
        FootprintAlignmentSummary summary = new FootprintAlignmentSummary(inputBamFile, outputFile, gencodeData);
        summary.createSummary();
//        
        TranscriptSummary sample = new TranscriptSummary("test", "type", 3, 3, 2, 2);
//        System.out.println(sample.getAlignmentPValueHelper(1, 1, 0));
        System.out.println(sample.alignmentPValue);
        
        for (int i = 7; i > -1; i--) {
            System.out.println("factorial " + i + " = " + factorial(i));
        }
        
        
    
    }
}
