package objects;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileInputStream;
import java.io.FileReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.util.zip.GZIPInputStream;

public class FastX {

	public static final int FILETYPE_FASTA = 1;
	public static final int FILETYPE_FASTQ = 2;

	//private File file = null;
	private int fileType = 0;



	public FastX(String filePath) throws IOException{
		this.open(new File(filePath));
	}
	public FastX(File file) throws IOException{
		this.open(file);
	}



	private BufferedReader in;
	private void open(File file) throws IOException{
		//this.file = file;
		
		// detect if this is a compressed file
		String[] tmp = file.getName().split("\\.");
		if(tmp[tmp.length-1].equalsIgnoreCase("gz")){
			//System.out.println("Reading compressed file");
			this.fileType = getFileType(new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file)))));
			in = new BufferedReader(new InputStreamReader(new GZIPInputStream(new FileInputStream(file))));
		}else{
			this.fileType = getFileType(new BufferedReader(new FileReader(file)));
			in = new BufferedReader(new FileReader(file));
		}		
	}
	public void close() throws IOException{
		this.in.close();
	}


	/**
	 * Read the next fasta/fastq record (if none exists, return null)
	 * 
	 * @return
	 * @throws IOException
	 */
	private String headerLine = null;
	boolean fileFinished = false;
	public FastX_Record readNext() throws IOException{
		String line = "";
		String newHeaderLine = null;
		String seq = "";
		String qual = "";
		FastX_Record result = null;
		//boolean readingQual = false;

		if(!fileFinished){

			if(this.fileType == FILETYPE_FASTA){
				// for the very first line in the file
				if(headerLine == null){
					headerLine = in.readLine().substring(1);
				}

				while((line=in.readLine()) != null){
					if(line.startsWith("@") || line.startsWith(">")){  // break
						newHeaderLine = line.substring(1);
						break;
					}else{
						seq = seq.concat(line);
					}
				}
				result = new FastX_Record(headerLine, seq);

				// indicate EOF to prevent infinite looping!
				if(line == null){
					fileFinished = true;
				}

				// remember this header line for the next sequence and return the current record
				headerLine = newHeaderLine;

			}
			else if(this.fileType == FILETYPE_FASTQ){
				while((line=in.readLine()) != null){
					if(line.startsWith("@") || line.startsWith(">")){  // break
						newHeaderLine = line.substring(1);
						//break;
						seq = in.readLine();
						in.readLine();
						qual = in.readLine();
						result = new FastX_Record(newHeaderLine, seq, qual);
						break;
					}/*else if(line.startsWith("+")){
					readingQual = true;
				}else{
					if(readingQual)
						qual = qual.concat(line);
					else
						seq = seq.concat(line);
				}*/
				}
				
				if(line == null){
					fileFinished = true;
				}
				
			}
			
			else{
				// don't know what kind of file this is
				return null;
			}


			return(result);


		}else{
			// EOF
			return null;
		}

	}








	/**
	 * Read the next fasta/fastq record (if none exists, return null)
	 * 
	 * @return
	 * @throws IOException
	 *//*
	public FastX_Record readNext() throws IOException{
		String line, id, seq, qual;
		if((line=in.readLine()) != null){
			if(this.fileType == FILETYPE_FASTA){
				if(line.startsWith("@") || line.startsWith(">")){
					id = line.substring(1);
					seq = in.readLine();
					return(new FastX_Record(id, seq));
				}else
					return null;

			}else if(this.fileType == FILETYPE_FASTQ){
				if(line.startsWith("@") || line.startsWith(">")){
					id = line.substring(1);
					seq = in.readLine();
					in.readLine();
					qual = in.readLine();
					return(new FastX_Record(id, seq, qual));
				}else
					return null;
			}else{
				return null;
			}
		}else{
			return null;
		}
	}*/



	/**
	 * 
	 * @param in
	 * @return
	 * @throws IOException
	 */
	public int getFileType(BufferedReader in) throws IOException{
		int fileType = 0;
		String line = in.readLine();
		String line2;
		if(line!=null){
			if(line.startsWith("@") || line.startsWith(">")){
				if((in.readLine()) != null){
					if((line2 = in.readLine()) != null){
						if(line2.startsWith("@") || line.startsWith(">")){
							/** is a FastA file... **/
							fileType = FILETYPE_FASTA;
						}else if (line2.startsWith("+")){
							/** is a FastQ file... **/
							fileType = FILETYPE_FASTQ;
						}else{
							System.err.println("Failed to ascertain the format of this file.");
							System.exit(0);
						}
					}
				}
			}
		}
		in.close();
		return fileType;
	}


}
