package footprintAlignments;

import java.io.BufferedOutputStream;
import java.io.IOException;
import net.sf.samtools.*;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileOutputStream;
import java.io.FileWriter;
import java.io.PrintStream;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Map;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.logging.Level;
import java.util.logging.Logger;
import java.util.Collections;
import java.util.HashSet;
import java.util.Set;

/**
 *
 * @author jdmcpeek
 */
public class bestAlignmentPatternIdentifier {
    private final File inputSortedBamFile;
    private final Map<String, String> transcriptToGeneMap;
    
    public bestAlignmentPatternIdentifier(File inputSortedBamFile, File gencodeDataFile) throws IOException {
        this.inputSortedBamFile = inputSortedBamFile;
        this.transcriptToGeneMap = new gtfParser(gencodeDataFile).getTranscriptGenes();
    }
    
    
    public Map<String, Map<String, Double>> getAlignmentMatrix() throws IOException {
        Map<String, ArrayList<String>> readToTranscriptMap = getReadToTranscriptMap();
        Map<String, ArrayList<SAMRecord>> transcriptToReadMap = getTranscriptToReadMap();
        Map<String, Map<String, Double>> alignmentMatrix = new HashMap();
        for (Map.Entry pair : readToTranscriptMap.entrySet()) {
            String readInQuestion = (String) pair.getKey();
            ArrayList<String> transcripts = (ArrayList<String>) pair.getValue();
            
            alignmentMatrix.put(readInQuestion, new HashMap<String, Double>());
            // i somehow need to get the offset value for this read on this specific transcript. how do i do that?
            // iterate through `transcripts`
            // get the offset of every read in each transcript
            //      for each read, tally up the number of other reads that align. Divide by the total. If there are no other reads on that transcript, this value should be... null?
            for (String transcript : transcripts) {
                // for every transcript that this read might belong to
                
                int[] bins = {0, 0, 0};  // tally the mod-3 values (buckets 0, 1, and 2)
                // other reads on that transcript
                ArrayList<SAMRecord> otherReads = transcriptToReadMap.get(transcript);
                int desiredBin = -1;  // initialize to something
                int readFrameBin;
                for (SAMRecord otherRead : otherReads) {
                    readFrameBin = otherRead.getAlignmentStart() % 3;
                    if (otherRead.getReadName() == null ? readInQuestion == null : otherRead.getReadName().equals(readInQuestion)) {
                        desiredBin = readFrameBin;
                    }
                    bins[readFrameBin] = bins[readFrameBin] + 1;
                }
                double fractionAligned = -1; 
                // if fractionAligned hasn't changed from 1, we're going to have an array exception, and we're going to put -1 on that cell of the matrix
                try {
                    // we're counting fractionAligned by how many reads were in the desired bucket, and dividing that value by the total number of reads. So if there's only 1 read, that's always 100% aligned
                    fractionAligned = bins[desiredBin] / (double) (bins[0] + bins[1] + bins[2]);
                } catch (ArrayIndexOutOfBoundsException ex) {
                    // we want to make sure out catch block doesn't actually quit the program. Null values (-1) are acceptable in our matrix because that signifies that the transcript doesn't contain that read!
                    Logger.getLogger(FootprintAlignmentSummary.class.getName()).log(Level.SEVERE, null, ex);
                    System.out.println("Something is wrong with your logic. There should exist a read that corresponds to your read in question.");
                } finally {
                    alignmentMatrix.get(readInQuestion).put(transcript, fractionAligned);
                }
            }
        }
        
        
        return alignmentMatrix; 
    
    }
   
    
    // the values of this hash table are arraylists that hold SAMRecords. 
    public Map<String, ArrayList<SAMRecord>> getTranscriptToReadMap() throws IOException {
        Map<String, ArrayList<SAMRecord>> transcriptToReadMap = new HashMap<String, ArrayList<SAMRecord>>();
        SAMFileReader inputBam = new SAMFileReader(inputSortedBamFile);
        
        for (final SAMRecord samRecord : inputBam) {
            String transcriptID = samRecord.getReferenceName();
            
            if (transcriptToReadMap.containsKey(transcriptID)) {
                transcriptToReadMap.get(transcriptID).add(samRecord);
            } else {
                transcriptToReadMap.put(transcriptID, new ArrayList<SAMRecord>(Arrays.asList(samRecord))); 
            }
        }
        return transcriptToReadMap;
    }
    
    public Map<String, ArrayList<String>> getReadToTranscriptMap() throws IOException {
        Map<String, ArrayList<String>> readToTranscriptMap = new HashMap<String, ArrayList<String>>();
        SAMFileReader inputBam = new SAMFileReader(inputSortedBamFile);
        for (final SAMRecord samRecord : inputBam) {
            String transcriptID = samRecord.getReferenceName();
            String readName = samRecord.getReadName();
            if (readToTranscriptMap.containsKey(readName)) {
                readToTranscriptMap.get(readName).add(transcriptID);
            } else {
                readToTranscriptMap.put(readName, new ArrayList<String>(Arrays.asList(transcriptID)));
            }            
        }
        // if a read matches transcripts that live on different genes, just remove it from the map, we don't need it.        
        removeCrossGeneReads(readToTranscriptMap); 
        return readToTranscriptMap;
    }
    
    // when are we supposed to remove a read? If there is any non-uniformity in genes that the transcripts belong to? I feel like that totally winnows down the set.
    public void removeCrossGeneReads(Map<String, ArrayList<String>> readToTranscriptMap) {
        Iterator it = readToTranscriptMap.entrySet().iterator();
        while (it.hasNext()) {
            Map.Entry pair = (Map.Entry)it.next();
            ArrayList<String> transcripts = (ArrayList<String>) pair.getValue();
            
            Set<String> uniqueGenes = new HashSet<String>();
            
            for (String transcript : transcripts) {
                String nextGene = transcriptToGeneMap.get(transcript);
                uniqueGenes.add(nextGene);
            }
            
            // If the read matches transcripts that together belong to multifarious genes, toss that read out. It is of no use to us. 
            if (uniqueGenes.size() > 1) {
                it.remove();
            }
            
        }
        
    }
    
    
    public static void main(String[] args) throws IOException { 
        System.out.println("bestAlignmentPatternIdentifier.java");
        
        bestAlignmentPatternIdentifier test;
        test = new bestAlignmentPatternIdentifier(new File("./data/fullFootPrintAlignment.sorted.bam"), new File("./data/gencode.v21.transcripts.gtf.txt"));
        
        Map<String, Map<String, Double>> alignmentMatrix = test.getAlignmentMatrix();
        
        System.setOut(new PrintStream(new BufferedOutputStream(new FileOutputStream("./data/test.txt"))));
        
        
        for (Map.Entry<String, Map<String, Double>> entry : alignmentMatrix.entrySet()) {
            System.out.println(entry.getKey());
            for (Map.Entry<String, Double> pair : entry.getValue().entrySet()) {
               System.out.println("\t" + pair.getKey() + " -> " + pair.getValue());
            }
            System.out.println();
            
        }
    }
}

