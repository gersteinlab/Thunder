package annotation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;

import objects.GTF;
import objects.GenomicCoordinate;

public class ReadGTF {

	/**
	 * Ascertain the source of the information in this GTF file -- relates to how to parse the columns
	 * @param gtfFile
	 * @return
	 */
	public static int getSourceOfGTF(String gtfFile){

		int thisType = 0;
		try{
			BufferedReader in = new BufferedReader(new FileReader(gtfFile));
			String line = "";
			while((line=in.readLine()) != null){
				if(line.startsWith("##")){
					if(line.startsWith("##provider: GENCODE")){
						thisType = GTF.GTF_TYPE_GENCODE;
					}
				}else{
					String[] bits = parseLine(line);
					if(bits[1].trim().equalsIgnoreCase("Cufflinks")){
						thisType = GTF.GTF_TYPE_CUFFMERGE;
					}

					break;
				}
			}
			in.close();
		}catch(IOException e){
			e.printStackTrace();
		}

		return thisType;
	}


	/**
	 * Read the records contained in this GTF file
	 * @param gtfType
	 * @param gtfFile
	 * @param keepAttributes
	 * @throws Exception
	 */
	public static GTF readGTF(String gtfFile, boolean keepAttributes, ArrayList<String> featureType, boolean collapseByTranscriptID, boolean suppressNs, ArrayList<String> additionalAttributes) throws Exception{
		System.out.print("Reading GTF file...");

		GTF thisGTF = new GTF();
		thisGTF.setType(getSourceOfGTF(gtfFile));
		
		BufferedReader in = new BufferedReader(new FileReader(gtfFile));

		String line = "";
		int count = 0;
		while((line=in.readLine()) != null){
			if(!line.startsWith("##")){
				if(parseLine(line, keepAttributes, thisGTF, featureType, suppressNs, additionalAttributes))
					count ++;
			}
		}
		in.close();
		System.out.println("Done- read "+count+" GTF entries.");

		if(collapseByTranscriptID){
			System.out.print("Sorting and collapsing GTF entries by transcriptID...");
			GTF tmp = GTF.collapseEntriesByTranscriptID(thisGTF);
			System.out.println("Done");
			return tmp;
		}
		
		return thisGTF;
	}


	/**
	 * Split each GTF record based on the file delimiter
	 * @param line
	 * @return
	 */
	private static String[] parseLine(String line){
		return line.split(" |\t");
	}

	
	/**
	 * Parses the useful information in each GTF record
	 * @param line
	 * @param keepAttributes
	 */
	private static boolean parseLine(String line, boolean keepAttributes, GTF thisGTF, ArrayList<String> featureType, boolean suppressNs, ArrayList<String> additionalAttributes){
		String tmpID = "";
		String[] bits = parseLine(line);

		// allow filtering by particular feature types (exon, transcript, CDS, etc.)
		if(featureType.size() == 0  ||  featureType.contains(bits[2])){
			tmpID = trimAttribute(bits[11].trim());

//			GenomicCoordinate tmp = new GenomicCoordinate(tmpID, bits[0], bits[6], Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
			
			//GenomicCoordinate(String id, String chrom, String source, String featureType, int start, int end, double score, String strand, int frame)
			GenomicCoordinate tmp = new GenomicCoordinate(tmpID, bits[0], bits[1], bits[2], Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue(), 0.0, bits[6], 0);
			// System.out.println(tmpID+" "+bits[0]+" "+bits[6]+" "+Integer.valueOf(bits[3]).intValue()+" "+Integer.valueOf(bits[4]).intValue());
			
			if(suppressNs == false)
				tmp.suppressNs = false;
			
			tmp.setAttributesToAddToFasta(additionalAttributes);
			
			if(keepAttributes)
				thisGTF.addCoordinate(addAttributes(tmp, bits));
			else
				thisGTF.addCoordinate(tmp);
			
			
			return true;
		}else{
			return false;
		}
	}



	/**
	 * Add information on each of the optional attributes for this GTF record
	 * @param tmp
	 * @param bits
	 * @return
	 */
	private static GenomicCoordinate addAttributes(GenomicCoordinate tmp, String[] bits){
		for(int i=8;i<bits.length;i+=2){
			tmp.addAttribute(bits[i].trim(), trimAttribute(bits[i+1].trim()));
		}
		return tmp;
	}



	/**
	 * Remove the leading '"' and trailing '";' from each attribute of this GTF line 
	 * @param in
	 * @return
	 */
	public static String trimAttribute(String in){
		//		return in.substring(1, in.length()-2);
		return in.replaceAll("^\"|\";$", "");
	}


}
