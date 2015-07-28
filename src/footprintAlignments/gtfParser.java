package footprintAlignments;

import java.io.File;
import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;
import java.util.logging.Level;
import java.util.logging.Logger;

/**
 *
 * @author J. David McPeek
 */
public class gtfParser {
    private final File gtfFile;
    private final HashMap<String, String> transcriptTypeMap = new HashMap<String, String>();
    
    
    // Wait, do I actually NEED unique transcripts here? No. I'll just store all transcripts in this hashmap.
    public gtfParser(File gtf) throws IOException {
        this.gtfFile = gtf;
        this.getTranscriptTypes();  // fills the transcriptTypeMap
    }
    
    
    public HashMap<String, String> getTypeMap() {
        return transcriptTypeMap; 
    }

    
    private void getTranscriptTypes() throws IOException { 

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
        } 
        
        
    }
    
    
    public static void main(String[] args) throws IOException {
            gtfParser gtf = new gtfParser(new File("./gencode.v21.transcripts.gtf.txt")); 
            gtf.getTranscriptTypes();
    }
    

}
