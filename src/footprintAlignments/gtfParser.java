package footprintAlignments;

import java.io.File;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Arrays;
import java.util.Iterator;
import java.util.Map;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author J. David McPeek
 */
public class gtfParser {
    private final File gtfFile;
<<<<<<< Updated upstream
   
=======
    private final HashMap<String, String> transcriptTypeMap = new HashMap<String, String>();
>>>>>>> Stashed changes
    
    public gtfParser(File gtf) throws IOException {
        this.gtfFile = gtf;
    }
    
    public HashMap<String, String> getTranscriptGenes() throws IOException {
        final HashMap<String, String> transcriptToGeneMap = new HashMap<String, String>();
        String thisEntry;
        String transcriptID;
        String geneID;
        BufferedReader gtf = null;
        
                    
        try {
            gtf = new BufferedReader(new FileReader(gtfFile));
            
            while((thisEntry = gtf.readLine()) != null) {

                String[] tokenizedEntry = thisEntry.split("\"");
                
                transcriptID = tokenizedEntry[3];  // index of the id
                geneID = tokenizedEntry[1];  // index of the gene id

                transcriptToGeneMap.put(transcriptID, geneID);
                
            } 
            
        } catch (IOException ex) {
            Logger.getLogger(gtfParser.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            if (gtf != null) gtf.close(); 
            
            return transcriptToGeneMap;
        } 
        
        
        
        
    }
    
    public HashMap<String, String> getTranscriptTypes() throws IOException { 
        
        final HashMap<String, String> transcriptTypeMap = new HashMap<String, String>();

        String thisEntry;
        String transcriptID;
        String transcriptType;
        BufferedReader gtf = null;
            
        try {
            gtf = new BufferedReader(new FileReader(gtfFile));
            
            while((thisEntry = gtf.readLine()) != null) {

                String[] tokenizedEntry = thisEntry.split("\"");

                transcriptID = tokenizedEntry[3];  // index of the id
                transcriptType = tokenizedEntry[11];  // index of the type

                transcriptTypeMap.put(transcriptID, transcriptType);
                
            } 
            
        } catch (IOException ex) {
            Logger.getLogger(gtfParser.class.getName()).log(Level.SEVERE, null, ex);
        } finally {
            if (gtf != null) gtf.close(); 
            
            return transcriptTypeMap;
        } 
        
        
    }
    
    public static void writeMapToFile(Map map, File file, String value) throws IOException {
        BufferedWriter outputFile = null;
        try {
            outputFile = new BufferedWriter(new FileWriter(file));
            Iterator it = map.entrySet().iterator();
            outputFile.write("\"transcript.id\"\t\"" + value + "\"");
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
    
    
    public static void main(String[] args) throws IOException {
            gtfParser gtf = new gtfParser(new File("./data/gencode.v21.transcripts.gtf.txt")); 
            HashMap test;
            test = gtf.getTranscriptGenes();
            writeMapToFile(test, new File("./data/transcriptToGene"), "gene.id"); 
            
            
    }
    

}
