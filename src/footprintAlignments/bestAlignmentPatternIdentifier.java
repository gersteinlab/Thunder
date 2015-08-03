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
            if (readToTranscriptMap.containsKey(transcriptID)) {
                readToTranscriptMap.get(transcriptID).add(readName);
            } else {
                readToTranscriptMap.put(transcriptID, new ArrayList<String>(Arrays.asList(readName)));
            }            
        }
        return readToTranscriptMap;
    }
    
    
    
    public static void main(String[] args) throws IOException { 
        System.out.println("test");
        
        bestAlignmentPatternIdentifier test;
        test = new bestAlignmentPatternIdentifier(new File("./data/fullFootPrintAlignment.sorted.bam"), new File("./data/gencode.v21.transcripts.gtf.txt"));
        
        test.getReadTranscriptMap();
    }
    
}

