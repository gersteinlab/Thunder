package footprintAlignments;

import net.sf.samtools.*;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
import java.io.IOException;
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
        protected String transcript;
        protected String type;
        protected int numberOfReads;
        protected double maxInFrame;        
        
        TranscriptSummary(String transcriptID, String transcriptType, int totalReads, double maxAligned) {
            transcript = transcriptID;
            type = transcriptType;
            numberOfReads = totalReads;
            maxInFrame = maxAligned;
        }
        
       public void writeSummaryToFile(BufferedWriter outputFile) throws IOException {
            outputFile.write(transcript + "\t");
            outputFile.write(type + "\t");
            outputFile.write(numberOfReads + "\t");
            outputFile.write(String.valueOf(maxInFrame + "\t"));
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
            
            String fileHeader = "\"transcript.id\"\t\"transcript.type\"\t\"total.reads\"\t\"aligned.max\"";
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

                String transcriptType = transcriptTypeMap.get(transcriptID); 
                TranscriptSummary thisSummary = new TranscriptSummary(transcriptID, transcriptType, totalReads, (double) maxAligned/totalReads);
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
     * @param args the command line arguments
     */
    public static void main(String[] args) throws IOException {

        File inputBamFile = new File("/Users/jdmcpeek/FootprintAlignments/data/sorted.bam");
        File outputFile = new File("./data/outTest.txt");
        File gencodeData = new File("./data/gencode.v21.transcripts.gtf.txt"); 
        
        FootprintAlignmentSummary summary = new FootprintAlignmentSummary(inputBamFile, outputFile, gencodeData);
        summary.createSummary();
    
    }
}
