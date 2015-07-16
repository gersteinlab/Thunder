package annotation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.InputStreamReader;
import java.util.ArrayList;

import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import objects.GTF;

public class ConvertCoordinates {

	private static boolean verbose = true;
	
	public static void transcriptome2Genome(String coordsPath, String gtfPath) throws Exception{
		ArrayList<String> features = new ArrayList<String>();
		features.add("exon");
		GTF annotation = ReadGTF.readGTF(gtfPath, true, features, true, false, null);

		/*
		 *  Read coordinates
		 *  assume columns are:
		 *  TranscriptID	chromosome	start	stop
		 */
		//HashMap<String, ArrayList<GenomicCoordinate>> coords = new HashMap<String, ArrayList<GenomicCoordinate>>();
		BufferedReader in = new BufferedReader(new FileReader(coordsPath));
		String line = "";
		String[] bits;
		//int count = 0;
		while((line=in.readLine()) != null){
			bits = line.split("\t");

			//if( ! coords.containsKey(bits[1]) ){
			//	coords.put(bits[1], new ArrayList<GenomicCoordinate>());
			//}
			//coords.get(bits[1]).add(new GenomicCoordinate(bits[0], bits[1], Integer.valueOf(bits[2]).intValue(), Integer.valueOf(bits[23]).intValue()));

			//annotation.getCoordinates().get(bits[1]).get
		}



	}


	/**
	 * 
	 * @param gtfPath
	 * @throws Exception
	 */
	public static void genome2Transcriptome(String gtfPath) throws Exception{
		ArrayList<String> features = new ArrayList<String>();
		features.add("exon");
		GTF annotation = ReadGTF.readGTF(gtfPath, true, features, true, false, null);

		BufferedReader in = new BufferedReader(new InputStreamReader(System.in));
		String line;
		String[] bits;
		while ((line = in.readLine()) != null && line.length() != 0){
			System.out.println(line);
			bits = line.split("\t");
			
			System.out.println(annotation.getCoordinates().containsKey(bits[3]));
		}
		// An empty line or Ctrl-Z terminates the program
	}


	
	/**
	 * 
	 * @param alignments
	 * @param gtfPath
	 * @throws Exception
	 */
	public static void genome2Transcriptome(File alignments, String gtfPath) throws Exception{
		ArrayList<String> features = new ArrayList<String>();
		features.add("exon");
		GTF annotation = ReadGTF.readGTF(gtfPath, true, features, true, false, null);

		SAMFileReader inputSam = new SAMFileReader(alignments);
		SAMRecord thisRecord;
		SAMRecordIterator it = inputSam.iterator();
		int count = 0;
		while(it.hasNext()){
			thisRecord = it.next();
			if(thisRecord.getIntegerAttribute("NH").equals(1)){
				if(verbose) System.out.print(thisRecord.getReferenceName()+":"+thisRecord.getAlignmentStart()+"-"+thisRecord.getAlignmentEnd()+"   "+thisRecord.getReadString());
				if(verbose) System.out.println("   NM:"+thisRecord.getIntegerAttribute("NM")+"   NH:"+thisRecord.getIntegerAttribute("NH")+"   "+thisRecord.getCigarString());
			}
			count++;
			if(count > 10)
				break;
		}



	}


	public static void main(String[] args) throws Exception {

		String gtfPath = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v18.annotation.gtf";
		//String samPath = "";

		genome2Transcriptome(gtfPath);


	}
	

}
