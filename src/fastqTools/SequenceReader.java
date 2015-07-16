package fastqTools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class SequenceReader {

	public static void main(String[] args) throws IOException {
		SequenceReader reader = new SequenceReader("/Users/robk/WORK/YALE_offline/ANNOTATIONS/test.fa");
		//SequenceReader reader = new SequenceReader("/Users/robk/Desktop/Box Documents/Astrid_HEEBOarrays/HEEBO_probes.fq");
		//SequenceReader reader = new SequenceReader(System.in);
		
		SequenceRecord tmp;
		while((tmp=reader.readNextRecord()) != null){
			System.out.println(tmp.toString());
		}

	}
	
	
	//
	private String thisSequenceID = null;
	private int nullCount = 0;
	private BufferedReader buffer;
	
	

	/**
	 * 
	 * @param inputFilePath
	 * @throws IOException
	 */
	public SequenceReader(String inputFilePath) throws IOException{
		buffer = new BufferedReader(new FileReader(inputFilePath));
	}
	
	public SequenceReader(InputStream in) throws IOException{
		buffer = new BufferedReader(new InputStreamReader(in));
	}
	
	
	

	/**
	 * 
	 * @return
	 * @throws IOException
	 */
	public SequenceRecord readNextRecord() throws IOException{
		String line;
		SequenceRecord toReturn = null;
		
		// Reads the first header in the file
		if(thisSequenceID == null){
			line=buffer.readLine();
			if(line.startsWith("@") || line.startsWith(">")){
				thisSequenceID = line.substring(1).trim();
			}
		}

		// Create a new object to hold the sequence about to be read
		toReturn = new SequenceRecord(thisSequenceID);
		
		boolean readingQual = false;
		
		// Read lines until the next fasta header
		while((line=buffer.readLine()) != null){
			if(line.startsWith("@") || line.startsWith(">")){
				thisSequenceID = line.substring(1).trim();
				break;
			}else if(line.startsWith("+")){
				readingQual = true;	
			}else{
				if(!readingQual)
					toReturn.addSequenceString(line.trim());
				else
					toReturn.addQualityString(line.trim());
			}
		}
		
		// hack to output the final sequence, then output null for the N+1th sequence
		if(line == null){
			if(nullCount == 0){
				nullCount ++;
			}else{
				toReturn = null;
			}
		}

		return toReturn;
	}

}
