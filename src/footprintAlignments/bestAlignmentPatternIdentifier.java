package footprintAlignments;

import java.io.IOException;
import net.sf.samtools.*;
import java.io.File;
import java.io.BufferedWriter;
import java.io.FileWriter;
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
    
    public Map<String, ArrayList<String>> getReadTranscriptMap() throws IOException {
        Map<String, ArrayList<String>> readToTranscriptMap; 
        readToTranscriptMap = new HashMap<String, ArrayList<String>>();
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
            
            for (int i = 0; i < transcripts.size(); i++) {
                String nextGene = transcriptToGeneMap.get(transcripts.get(i));
                uniqueGenes.add(nextGene);
            }
            
            if (uniqueGenes.size() > 1) {
                it.remove();
            }
            
        }
        
    }
    
    
    public static void main(String[] args) throws IOException { 
        System.out.println("test");
        
        bestAlignmentPatternIdentifier test;
        test = new bestAlignmentPatternIdentifier(new File("./data/fullFootPrintAlignment.sorted.bam"), new File("./data/gencode.v21.transcripts.gtf.txt"));
        
        Map<String, ArrayList<String>> readTranscriptMap = test.getReadTranscriptMap(); 
        
        
        Iterator it = readTranscriptMap.entrySet().iterator();
        
        BufferedWriter testAgain;
        testAgain = new BufferedWriter(new FileWriter("./data/testAgain.txt"));
        
        while(it.hasNext()){
            Map.Entry pair = (Map.Entry)it.next();
            testAgain.write(pair.getKey() + "  -  " + pair.getValue()); 
            testAgain.newLine();
        }
        
        
    }
    
}

