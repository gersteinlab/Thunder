package objects;

import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

public class GTF {

	public static final int GTF_TYPE_GENCODE = 1;
	public static final int GTF_TYPE_CUFFMERGE = 2;
	private int gtfType = 0;
	public void setType(int type){ this.gtfType = type; }
	public int getType(){ return this.gtfType; }
	
	//private HashMap<Integer, >
	
	// transcripts by chromosome:
	private HashMap<String, ArrayList<GenomicCoordinate>> list_GTFcoords = new HashMap<String, ArrayList<GenomicCoordinate>>();
	
	
	public void setCoordinates(HashMap<String, ArrayList<GenomicCoordinate>> gtfCoords){ this.list_GTFcoords = gtfCoords; }
	public HashMap<String, ArrayList<GenomicCoordinate>> getCoordinates(){ return this.list_GTFcoords; }

	
	private HashMap<String, GenomicCoordinate> list_GTFcoords_byTranscript = new HashMap<String, GenomicCoordinate>();
	public HashMap<String, GenomicCoordinate> getCoordinates_byTranscript(){ return this.list_GTFcoords_byTranscript; }

	/**
	 * Add a GenomicCoordinate object to the hashmap
	 * @param coord
	 */
	public void addCoordinate(GenomicCoordinate coord){

		if( ! this.list_GTFcoords.containsKey(coord.getChrom())){
			this.list_GTFcoords.put(coord.getChrom(), new ArrayList<GenomicCoordinate>());
		}
		this.list_GTFcoords.get(coord.getChrom()).add(coord);

		this.list_GTFcoords_byTranscript.put(coord.getCoordinateID(), coord);

	}






	/**
	 * Sort and collapse transcripts on each chromosome by TRANSCRIPT ID
	 * 
	 * @param allCoords
	 * @return
	 */
	public static GTF collapseEntriesByTranscriptID(GTF gtfContents){
		Iterator<String> it = gtfContents.getCoordinates().keySet().iterator();
		HashMap<String, ArrayList<GenomicCoordinate>> newCoords = new HashMap<String, ArrayList<GenomicCoordinate>>();
		String tmpChromosome;
		while(it.hasNext()){
			tmpChromosome = it.next();
			//System.out.println("Collapsing "+allCoords.get(tmpChromosome).size()+" sequences on chromosome: "+tmpChromosome);
			newCoords.put(tmpChromosome, collapseEntriesByTranscriptID(gtfContents.getCoordinates().get(tmpChromosome)));
		}

		GTF newGTF = new GTF();
		newGTF.setType(gtfContents.getType());
		newGTF.setCoordinates(newCoords);
		return newGTF;
	}



	/**
	 * Sort and collapse transcripts on ONE chromosome by TRANSCRIPT ID
	 * 
	 * @param coordsOnThisChromosome
	 * @return
	 */
	public static ArrayList<GenomicCoordinate> collapseEntriesByTranscriptID(ArrayList<GenomicCoordinate> coordsOnThisChromosome){

		Collections.sort(coordsOnThisChromosome, new CoordinateIDComparator());

		ArrayList<GenomicCoordinate> coordsOnThisChromosome_collapsed = new ArrayList<GenomicCoordinate>();

		GenomicCoordinate thisCoord;
		coordsOnThisChromosome_collapsed.add(coordsOnThisChromosome.get(0));

		String lastCoordinateID = coordsOnThisChromosome.get(0).getAttribute("transcript_id");

		int usingCoordinateAtIndex = 0;
		for(int i=1;i<coordsOnThisChromosome.size();i++){
			thisCoord = coordsOnThisChromosome.get(i);

			if(thisCoord.getAttribute("transcript_id").equals(lastCoordinateID)){
				// update existing entry
				coordsOnThisChromosome_collapsed.get(usingCoordinateAtIndex).mergeCoordinates(thisCoord);
			}else{
				// retain entry -- update trackers and add entry to new list
				coordsOnThisChromosome_collapsed.add(thisCoord);
				lastCoordinateID = thisCoord.getAttribute("transcript_id");
				usingCoordinateAtIndex ++;
			}
		}
		return coordsOnThisChromosome_collapsed;
	}






	/**
	 * Sort and collapse transcripts on each chromosome by GENE ID
	 * 
	 * @param allCoords
	 * @return
	 */
	public static GTF collapseEntriesByGeneID(GTF gtfContents){
		Iterator<String> it = gtfContents.getCoordinates().keySet().iterator();
		HashMap<String, ArrayList<GenomicCoordinate>> newCoords = new HashMap<String, ArrayList<GenomicCoordinate>>();
		String tmpChromosome;
		while(it.hasNext()){
			tmpChromosome = it.next();
			//			System.out.println("Collapsing "+allCoords.get(tmpChromosome).size()+" sequences on chromosome: "+tmpChromosome);
			newCoords.put(tmpChromosome, collapseEntriesByGeneID(gtfContents.getCoordinates().get(tmpChromosome)));
		}

		GTF newGTF = new GTF();
		newGTF.setType(gtfContents.getType());
		newGTF.setCoordinates(newCoords);
		return newGTF;

	}



	/**
	 * Sort and collapse transcripts on ONE chromosome by GENE ID
	 * 
	 * @param coordsOnThisChrom
	 * @return
	 */
	public static ArrayList<GenomicCoordinate> collapseEntriesByGeneID(ArrayList<GenomicCoordinate> coordsOnThisChrom){

		Collections.sort(coordsOnThisChrom, new GeneIDComparator());

		ArrayList<GenomicCoordinate> coordsOnThisChromosome_collapsed = new ArrayList<GenomicCoordinate>();

		GenomicCoordinate thisCoord;
		coordsOnThisChromosome_collapsed.add(coordsOnThisChrom.get(0));

		String lastCoordinateID = coordsOnThisChrom.get(0).getAttribute("gene_id");

		int usingCoordinateAtIndex = 0;
		for(int i=1;i<coordsOnThisChrom.size();i++){
			thisCoord = coordsOnThisChrom.get(i);

			if(thisCoord.getAttribute("gene_id").equals(lastCoordinateID)){
				// update existing entry
				coordsOnThisChromosome_collapsed.get(usingCoordinateAtIndex).mergeCoordinates(thisCoord);
			}else{
				// retain entry -- update trackers and add entry to new list
				coordsOnThisChromosome_collapsed.add(thisCoord);
				lastCoordinateID = thisCoord.getAttribute("gene_id");
				usingCoordinateAtIndex ++;
			}
		}
		return coordsOnThisChromosome_collapsed;
	}


}
