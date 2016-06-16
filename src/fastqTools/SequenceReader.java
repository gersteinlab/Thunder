package fastqTools;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStream;
import java.io.InputStreamReader;

public class SequenceReader {

	public static void main(String[] args) throws IOException {
		//SequenceReader reader = new SequenceReader("/Users/robk/WORK/YALE_offline/ANNOTATIONS/test.fa");
		SequenceReader reader = new SequenceReader("/Users/robk/Downloads/test.fq");
		//SequenceReader reader = new SequenceReader(System.in);

		SequenceRecord tmp;
		while((tmp=reader.readNextRecord()) != null){
			System.out.println(tmp.toString());
		}

	}


	//
	private String thisSequenceID = null;
	private int nullCount = 0;
	private BufferedReader _buffer;



	/**
	 * 
	 * @param inputFilePath
	 * @throws IOException
	 */
	public SequenceReader(String inputFilePath) throws IOException{
		if(inputFilePath.trim().equals("-"))
			_buffer = new BufferedReader(new InputStreamReader(System.in));
		else
			_buffer = new BufferedReader(new FileReader(inputFilePath));
	}

	public SequenceReader(InputStream in) throws IOException{
		_buffer = new BufferedReader(new InputStreamReader(in));
	}


	public void close() throws IOException{
		_buffer.close();
	}

	private boolean isFastq = true;
	private int lineCount = 0;
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
			line=_buffer.readLine();
			//if(line.startsWith("@") || line.startsWith(">")){
			if(line.startsWith("@")){
				isFastq = true;
				lineCount ++;
			}else if(line.startsWith(">"))
				isFastq = false;

			if(line.startsWith("@") || line.startsWith(">")){
				thisSequenceID = line.substring(1).trim();
			}
		}

		// Create a new object to hold the sequence about to be read
		toReturn = new SequenceRecord(thisSequenceID);

		boolean readingQual = false;

		// Read lines until the next header
		while((line=_buffer.readLine()) != null){
			lineCount ++;
			//if(line.startsWith("@") || line.startsWith(">")){
			if(isFastq){
				if(lineCount == 5){
					thisSequenceID = line.substring(1).trim();
					lineCount = 1;
					break;
				}else if(line.startsWith("+")  &&  lineCount == 3){
					readingQual = true;	
				}else{
					if(!readingQual)
						toReturn.addSequenceString(line.trim());
					else
						toReturn.addQualityString(line.trim());
				}
			}else{
				if(line.startsWith(">")){
					thisSequenceID = line.substring(1).trim();
					break;
				}else{
					toReturn.addSequenceString(line.trim());
				}
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
