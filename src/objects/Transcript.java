package objects;

import java.util.ArrayList;
import java.util.Iterator;

public class Transcript {

	private String _id, _chromosome, _strand, _transcriptBiotype, _geneID;
	private int _exonCount = 0;
	private int _cdsCount = 0;
	private ArrayList<GenomicCoordinate> _exons = new ArrayList<GenomicCoordinate>();
	private ArrayList<GenomicCoordinate> _cds = new ArrayList<GenomicCoordinate>();

	public Transcript(String id, String chromosome, String strand, String transcriptBiotype, String geneID){
		_id = id;
		_chromosome = chromosome;
		_strand = strand;
		_transcriptBiotype = transcriptBiotype;
		_geneID = geneID;
	}

	public void addExon(int start, int stop){
		_exonCount ++;
		_exons.add(new GenomicCoordinate("exon_"+_exonCount, _chromosome, start, stop));
	}
	public void addExon(GenomicCoordinate coords){
		_exonCount ++;
		_exons.add(new GenomicCoordinate("exon_"+_exonCount, _chromosome, coords.getStart(), coords.getStop()));
	}

	public void addCDS(int start, int stop){
		_cdsCount ++;
		_cds.add(new GenomicCoordinate("CDS_"+_cdsCount, _chromosome, start, stop));
	}
	public void addCDS(GenomicCoordinate coords){
		_cdsCount ++;
		_cds.add(new GenomicCoordinate("CDS_"+_exonCount, _chromosome, coords.getStart(), coords.getStop()));
	}


	public String getID(){ return _id; }
	public String getGeneID(){ return _geneID; }
	public String getStrand(){ return _strand; }
	public String getChromosome(){ return _chromosome; }
	public String getTranscriptBiotype(){ return _transcriptBiotype; }
	public int getNumberOfExons(){ return _exonCount; }
	public int getNumberOfCodingExons(){ return _cdsCount; }


	private int getTotalLength(ArrayList<GenomicCoordinate> coords){
		int result = 0;
		Iterator<GenomicCoordinate> it = coords.iterator();
		while(it.hasNext()){
			GenomicCoordinate thisExon = it.next();
			result += thisExon.getStop() - thisExon.getStart();
		}
		return result;
	}
	/**
	 * Compute the total length (nt) of this transcript
	 * @return
	 */
	public int getTotalExonLength(){ return getTotalLength(_exons); }
	/**
	 * Compute the total length (nt) of the CDS of this transcript
	 * @return
	 */
	public int getTotalCodingExonLength(){ return getTotalLength(_cds); }


	/**
	 * Get the START position of the given exon type in genomic coordinates
	 * @return
	 */
	private int getStart(ArrayList<GenomicCoordinate> coords){
		int start = -1;
		if(coords.size() > 0){
			GenomicCoordinate first = null;
			GenomicCoordinate tmp;
			Iterator<GenomicCoordinate> it = coords.iterator();
			while(it.hasNext()){
				tmp = it.next();
				if(first == null){
					first = tmp;
				}else if(_strand.equals("+")  &&  tmp.getStart() < first.getStart()){
					first = tmp;
				}else if(_strand.equals("-")  &&  tmp.getStop() > first.getStop()){
					first = tmp;
				}
			}
			if(_strand.equals("+"))
				start = first.getStart();
			else if(_strand.equals("-"))
				start = first.getStop();
		}
		return start;
	}

	/**
	 * Get the STOP position of the given exon type in genomic coordinates
	 * @return
	 */
	private int getStop(ArrayList<GenomicCoordinate> coords){
		int stop = -1;
		if(coords.size() > 0){
			GenomicCoordinate last = null;
			GenomicCoordinate tmp;
			Iterator<GenomicCoordinate> it = coords.iterator();
			while(it.hasNext()){
				tmp = it.next();
				if(last == null){
					last = tmp;
				}else if(_strand.equals("+")  &&  tmp.getStop() > last.getStop()){
					last = tmp;
				}else if(_strand.equals("-")  &&  tmp.getStart() < last.getStart()){
					last = tmp;
				}
			}
			if(_strand.equals("+"))
				stop = last.getStop();
			else if(_strand.equals("-"))
				stop = last.getStart();
		}
		return stop;
	}

	/**
	 * Get the genomic position of the START of this transcript
	 * @return
	 */
	public int getTranscriptStart(){ return getStart(_exons); }
	/**
	 * Get the genomic position of the END of this transcript
	 * @return
	 */
	public int getTranscriptStop(){ return getStop(_exons); }
	/**
	 * Get the genomic position of the START of the CDS of this transcript
	 * @return
	 */
	public int getCDSStart(){ return getStart(_cds); }
	/**
	 * Get the genomic position of the END of the CDS of this transcript
	 * @return
	 */
	public int getCDSStop(){ return getStop(_cds); }


	/**
	 * Compute the distance (in transcript coordinates) from the START of the FIRST exon to the START of the CDS
	 * @return
	 */
	public int getCDSStartFromTxStartInTxCoords(){
		int distance = -1;
		if(_cdsCount > 0){
			distance = 0;
			int cdsStart = getCDSStart();
			Iterator<GenomicCoordinate> it = _exons.iterator();

			if(_strand.equals("+")){
				while(it.hasNext()){
					GenomicCoordinate tmp = it.next();
					if(tmp.getStop() < cdsStart){ // CDS starts after this exon, include the whole thing
						distance += tmp.getStop()-tmp.getStart();
					}else if(tmp.getStart() < cdsStart){ // CDS starts inside this exon, include only the upstream part
						distance += cdsStart-tmp.getStart();
					}
				}
			}else if(_strand.equals("-")){
				while(it.hasNext()){
					GenomicCoordinate tmp = it.next();
					if(tmp.getStart() > cdsStart){ // CDS starts after this exon, include the whole thing
						distance += tmp.getStop()-tmp.getStart();
					}else if(tmp.getStop() > cdsStart){ // CDS starts inside this exon, include only the upstream part
						distance += tmp.getStop()-cdsStart;
					}
				}
			}
		}
		return distance;
	}

	/**
	 * Compute the distance (in transcript coordinates) from the END of the LAST exon to the END of the CDS
	 * @return
	 */
	public int getCDSStopFromTxStopInTxCoords(){
		int distance = -1;
		if(_cdsCount > 0){
			distance = 0;
			int cdsStop = getCDSStop();
			Iterator<GenomicCoordinate> it = _exons.iterator();

			if(_strand.equals("+")){
				while(it.hasNext()){
					GenomicCoordinate tmp = it.next();
					if(tmp.getStart() > cdsStop){ // CDS starts after this exon, include the whole thing
						distance += tmp.getStop()-tmp.getStart();
					}else if(tmp.getStop() > cdsStop){ // CDS starts inside this exon, include only the upstream part
						distance += tmp.getStop()-cdsStop;
					}
				}
			}else if(_strand.equals("-")){
				while(it.hasNext()){
					GenomicCoordinate tmp = it.next();
					if(tmp.getStop() < cdsStop){ // CDS starts after this exon, include the whole thing
						distance += tmp.getStop()-tmp.getStart();
					}else if(tmp.getStart() < cdsStop){ // CDS starts inside this exon, include only the upstream part
						distance += cdsStop-tmp.getStart();
					}
				}
			} 
		}
		return distance;
	}



	/**
	 * Get tab-delimited summary information about this transcript
	 */
	public String toString(){ return toString("\t"); }
	public String toString(String sep){
		String out = _id +sep+ _geneID +sep+ getTranscriptBiotype() +sep+ getChromosome() +sep+ getStrand() +sep+ getNumberOfExons() +sep+ getNumberOfCodingExons() +sep+ getTotalExonLength() +sep+ getTotalCodingExonLength();
		return out;
	}




}
