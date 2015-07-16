package transcriptome;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import objects.GTF;
import objects.GenomicCoordinate;
import objects.GenomicCoordinateComparator;

public class Sextractor {

	public Sextractor(){}

	
	/**
	 * 
	 * @param coords
	 * @return
	 */
	public static GenomicCoordinate parseCoordString(String coords){
		return parseCoordString("no_id", coords);
	}
	
	
	
	/**
	 * 
	 * @param id
	 * @param coords
	 * @return
	 */
	public static GenomicCoordinate parseCoordString(String id, String coords){
		String chrom = "";
		int start = -1;
		int end = -1;

		String[] bits = coords.split(":");

		chrom = bits[0]; 

		bits = bits[1].split("-");
		start = Integer.valueOf(bits[0]).intValue();
		end = Integer.valueOf(bits[1]).intValue();

		return new GenomicCoordinate(id, chrom, start, end);
	}



	/**
	 * 
	 */
	private HashMap<String, ArrayList<GenomicCoordinate>> map_coordsToExtract = new HashMap<String, ArrayList<GenomicCoordinate>>();
	public void setCoords(HashMap<String, ArrayList<GenomicCoordinate>> newCoords){
		this.map_coordsToExtract = newCoords;
	}

	
	
	/**
	 * 
	 * @param id
	 * @param chromosome
	 * @param strand
	 * @param start
	 * @param end
	 */
	public void addCoord(String id, String chromosome, String strand, int start, int end){
		addCoord(new GenomicCoordinate(id, chromosome, strand, start, end));
	}


	/**
	 * 
	 * @param coord
	 */
	public void addCoord(GenomicCoordinate coord){
		if( ! map_coordsToExtract.containsKey(coord.getChrom())){
			map_coordsToExtract.put(coord.getChrom(), new ArrayList<GenomicCoordinate>());
		}
		map_coordsToExtract.get(coord.getChrom()).add(coord);
	}
	


	/**
	 * 
	 */
//	public void collapseEntriesByTranscriptID(){
//		Iterator<String> mapIterator = map_coordsToExtract.keySet().iterator();
//		while(mapIterator.hasNext()){
//			String thisChrom = mapIterator.next();
//			this.map_coordsToExtract.put(thisChrom, collapseEntriesByTranscriptID(map_coordsToExtract.get(thisChrom), thisChrom));
//		}
//
//	}

	/**
	 * 
	 * @param coordsOnThisChromosome
	 * @param thisChrom
	 * @return
	 */
//	public ArrayList<GenomicCoordinate> collapseEntriesByTranscriptID(ArrayList<GenomicCoordinate> coordsOnThisChromosome, String thisChrom){
//
//		System.out.println("Collapsing "+coordsOnThisChromosome.size()+" sequences on chromosome: "+thisChrom);
//
//		Collections.sort(coordsOnThisChromosome, new CoordinateIDComparator());
//
//		ArrayList<GenomicCoordinate> coordsOnThisChromosome_collapsed = new ArrayList<GenomicCoordinate>();
//
//		GenomicCoordinate thisCoord;
//		coordsOnThisChromosome_collapsed.add(coordsOnThisChromosome.get(0));
//
//		String lastCoordinateID = coordsOnThisChromosome.get(0).getCoordinateID();
//
//		int usingCoordinateAtIndex = 0;
//		for(int i=1;i<coordsOnThisChromosome.size();i++){
//			thisCoord = coordsOnThisChromosome.get(i);
//
//			if(thisCoord.getCoordinateID().equals(lastCoordinateID)){
//				// update existing entry
//				coordsOnThisChromosome_collapsed.get(usingCoordinateAtIndex).mergeCoordinates(thisCoord);
//			}else{
//				// retain entry -- update trackers and add entry to new list
//				coordsOnThisChromosome_collapsed.add(thisCoord);
//				lastCoordinateID = thisCoord.getCoordinateID();
//				usingCoordinateAtIndex ++;
//			}
//
//		}
//
//		return coordsOnThisChromosome_collapsed;
//	}


	/**
	 * 
	 * @param outPath
	 * @throws IOException
	 */
	public void writeExtractedSequences(String outPath) throws IOException{
		Writer out = null;
		if(outPath != null)
			out = new BufferedWriter(new FileWriter(outPath));

		Iterator<String> mapIterator = map_coordsToExtract.keySet().iterator();
		while(mapIterator.hasNext()){
			String writeChrom = mapIterator.next();
			Iterator<GenomicCoordinate> tmpIt = map_coordsToExtract.get(writeChrom).iterator();
			System.out.println("Writing "+map_coordsToExtract.get(writeChrom).size()+" sequences on chromosome: "+writeChrom);
			while(tmpIt.hasNext()){
				GenomicCoordinate thisCoord = tmpIt.next();
				if(outPath != null){
					out.write(thisCoord.toString()+"\n");
				}else{
					System.out.println(thisCoord.toString());
				}
			}
		}

		if(outPath != null){
			out.flush();
			out.close();
		}
	}



	/**
	 * 
	 * @param fastaPath
	 * @param collapseByID
	 * @param outputToPath
	 * @throws IOException
	 */
	public void extractFromFasta(String fastaPath, boolean collapseByID, String outputToPath) throws IOException{
//	public void extractFromFasta(String fastaPath, String outputToPath) throws IOException{
		Writer out = null;
		if(outputToPath != null){
			out = new BufferedWriter(new FileWriter(outputToPath));
			out.flush();
			out.close();
		}
		
		BufferedReader in = new BufferedReader(new FileReader(fastaPath));
		String line = "";
		String thisChromosome = "";
		boolean coordInThisChromosome = false;
		int currentPosition = 0;
		int startIndex = 0;

		ArrayList<GenomicCoordinate> coordsOnThisChromosome = new ArrayList<GenomicCoordinate>();

		while((line=in.readLine()) != null){
			
			if(line.startsWith(">")){
				
				if(coordsOnThisChromosome.size() > 0){
					

					// collapse transcripts on this chromosome?
					if(collapseByID){
//						coordsOnThisChromosome = collapseEntriesByTranscriptID(coordsOnThisChromosome, thisChromosome);
						System.out.println("Collapsing "+coordsOnThisChromosome.size()+" sequences on chromosome: "+thisChromosome);
						coordsOnThisChromosome = GTF.collapseEntriesByTranscriptID(coordsOnThisChromosome);
					}

					// write transcripts to fasta?
					if(outputToPath != null){
						map_coordsToExtract.remove(thisChromosome);
						
						out = new BufferedWriter(new FileWriter(outputToPath, true));
						
						Iterator<GenomicCoordinate> tmpIt = coordsOnThisChromosome.iterator();
						System.out.println("Writing "+coordsOnThisChromosome.size()+" sequences on chromosome: "+thisChromosome);
						while(tmpIt.hasNext()){
//							System.out.println(tmpIt.next().toString());
							out.write(tmpIt.next().toString()+"\n");
						}
						out.flush();
						out.close();
						
						coordsOnThisChromosome = null;
						System.gc();
					}else{
						map_coordsToExtract.put(thisChromosome, coordsOnThisChromosome);
					}
				}

				// update new chromosome name
				thisChromosome = line.substring(1);
				System.out.println("Reading chromosome: "+thisChromosome);

				if(map_coordsToExtract.containsKey(thisChromosome)){
					coordInThisChromosome = true;

					currentPosition = 1;
					coordsOnThisChromosome = this.map_coordsToExtract.get(thisChromosome);
					startIndex = 0;
					System.out.println("Extracting "+coordsOnThisChromosome.size()+" sequences from chromosome: "+thisChromosome);

					Collections.sort(coordsOnThisChromosome, new GenomicCoordinateComparator());

				}else{
					coordInThisChromosome = false;
					coordsOnThisChromosome = new ArrayList<GenomicCoordinate>();
				}
			}
			else if(coordInThisChromosome){

 
				for(int i=startIndex;i<coordsOnThisChromosome.size();i++){
					//				while(tmpIt.hasNext()){
					//					GenomicCoordinate thisCoord = tmpIt.next();
					GenomicCoordinate thisCoord = coordsOnThisChromosome.get(i);
					if(thisCoord.overlapsWith(currentPosition, currentPosition+line.length())){
						//						System.out.println("startIndex = "+startIndex+"\t"+coordsOnThisChromosome.get(startIndex));
						//						System.out.println(thisCoord.toString()+" -- "+line);
						int[] tmp = thisCoord.intersectWith(currentPosition, line.length());
						//						System.out.println(tmp[0]+" -- "+tmp[1]+":\t"+line.substring(tmp[0]-1, tmp[1]));
						try{
							thisCoord.appendSequence(line.substring(tmp[0]-1, tmp[1]));
						}catch(Exception e){
							System.err.println("thisCoord.getCoordinateID()="+thisCoord.getCoordinateID());
							System.err.println("thisCoord.getStart()="+thisCoord.getStart());
							System.err.println("thisCoord.getStop()="+thisCoord.getStop());
							
							//System.err.println("line="+line);
							
							System.err.println("currentPosition="+currentPosition);
							System.err.println("tmp[0]="+tmp[0]);
							System.err.println("tmp[1]="+tmp[1]);
							
							e.printStackTrace();
							System.exit(0);
						}

					}else if (thisCoord.getStart() > currentPosition+line.length()){
						break;
					}

				}

				if(startIndex < coordsOnThisChromosome.size()){
					while(coordsOnThisChromosome.get(startIndex).hasCompleteSequence()){
						startIndex ++;
						if(startIndex >= coordsOnThisChromosome.size())
							break;
					}
				}

				currentPosition += line.length();
			}	


		}

		in.close();
		
		
		
		
		/*
		 * Do last chromosome extraction
		 */
		if(coordsOnThisChromosome.size() > 0){
			
			// collapse transcripts on this chromosome:?
			if(collapseByID){
//				coordsOnThisChromosome = collapseEntriesByTranscriptID(coordsOnThisChromosome, thisChromosome);
				System.out.println("Collapsing "+coordsOnThisChromosome.size()+" sequences on chromosome: "+thisChromosome);
				coordsOnThisChromosome = GTF.collapseEntriesByTranscriptID(coordsOnThisChromosome);
			}

			// write transcripts to fasta?
			if(outputToPath != null){
				map_coordsToExtract.remove(thisChromosome);
				
				out = new BufferedWriter(new FileWriter(outputToPath, true));
				
				Iterator<GenomicCoordinate> tmpIt = coordsOnThisChromosome.iterator();
				System.out.println("Writing "+coordsOnThisChromosome.size()+" sequences on chromosome: "+thisChromosome);
				while(tmpIt.hasNext()){
					out.write(tmpIt.next().toString()+"\n");
				}
				out.flush();
				out.close();
				
				coordsOnThisChromosome = null;
				System.gc();
			}else{
				map_coordsToExtract.put(thisChromosome, coordsOnThisChromosome);
			}
		}
		
		
		
		if(outputToPath != null){
			out.close();
		}
		
	}






	/**
	 * 
	 * @param args
	 * @throws IOException
	 */
	public static void main(String[] args) throws IOException {

		String coordinate = "chr14:90178855-90179167";
		String fastaFile = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/hg19.fa";

		//		String coordinate = "chr1:10-30";
		//		String fastaFile = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/test.fa";

		GenomicCoordinate coords = parseCoordString("TCONS_00087140|ENST00000363442.1|7SK", coordinate);
		//		System.out.println(coords.toString());

		Sextractor thisInstance = new Sextractor();
		thisInstance.addCoord(coords);

		thisInstance.addCoord(parseCoordString("TCONS_00289398|ENST00000581458.1|Metazoa_SRP", "chr9:9442060-9442347"));
		thisInstance.addCoord(parseCoordString("TCONS_00045155|ENST00000501122.2|NEAT1", "chr11:65190269-65213011"));


		thisInstance.extractFromFasta(fastaFile, false, null);


		thisInstance.writeExtractedSequences(null);

	}


}



