package objects;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

public class GenomicCoordinate{
	int start, end;
	String id = "NULL";
	String chrom;
	String strand = "+";
	String nucleotideSequence = "";

	String source = "";
	String featureType = "";
	double score;
	int frame;

	private ArrayList<String> attributesToAddToFasta = null;
	public void setAttributesToAddToFasta(ArrayList<String> attributeIDs){ 
		attributesToAddToFasta = attributeIDs;
	}


	public GenomicCoordinate(String chrom, int start, int end){
		this.chrom = chrom;
		this.start = start;
		this.end = end;
	}

	public GenomicCoordinate(String id, String chrom, int start, int end){
		this.id = id;
		this.chrom = chrom;
		this.start = start;
		this.end = end;
	}

	public GenomicCoordinate(String id, String chrom, String strand, int start, int end){
		this.id = id;
		this.chrom = chrom;
		this.strand = strand;
		this.start = start;
		this.end = end;
	}

	public GenomicCoordinate(String id, String chrom, String source, String featureType, int start, int end, double score, String strand, int frame){
		this.id = id;
		this.chrom = chrom;
		this.source = source;
		this.featureType = featureType;
		this.start = start;
		this.end = end;
		this.score = score;
		this.strand = strand;
		this.frame = frame;
	}

	public String getCoordinateID(){ return this.id; }
	public String getChrom(){ return this.chrom; }
	public String getSource(){ return this.source; }
	public String getFeatureType(){ return this.featureType; }
	public Integer getStart(){ return this.start; }
	public Integer getStop(){ return this.end; }
	public Double getScore(){ return this.score; }
	public String getStrand(){ return this.strand; }
	public Integer getFrame(){ return this.frame; }




	private HashMap<String, String> attributes = new HashMap<String, String>();
	public void addAttribute(String key, String value){
		this.attributes.put(key, value);
	}
	public boolean hasAttribute(String key){
		if(this.attributes.containsKey(key))
			return true;
		else
			return false;
	}
	public String getAttribute(String key){
		if(this.attributes.containsKey(key)){
			return this.attributes.get(key);
		}else{
			return null;
		}
	}
	public void setAttribute(String key, String newValue){
		if(hasAttribute(key)){
			this.attributes.put(key, newValue);
		}
	}



	/**
	 * Appends new sequence to the one currently stored for this entry
	 * @param seq
	 */
	public void appendSequence(String seq){ this.nucleotideSequence += seq; }


	/**
	 * Returns the sequence for this entry (can be the reverse complement if on the negative strand)
	 * 
	 * @param requireStrandedness
	 * @return
	 */
	public String getSequence(boolean requireStrandedness){
		if(requireStrandedness && this.strand.equals("-")){
			return reverseComplement(this.nucleotideSequence);
		}else{
			return this.nucleotideSequence;
		}
	}


	/**
	 * Calculates the reverse complement of the given sequence
	 * 
	 * @param seq
	 * @return
	 */
	static String reverseComplement(String seq){
		final StringBuilder out = new StringBuilder(seq.length());
		for (int i = seq.length(); i > 0; i--)
			try{
				out.append(replace(seq.charAt(i-1)));
			}catch(IllegalArgumentException e){
				out.append("N");
			}
		return out.toString();
	}

	/**
	 * Specifies the reverse complement of each nucleotide 
	 * 
	 * @param in
	 * @return
	 */
	public static char replace(char in) {
		switch (in) {
		case 'A': return 'T';
		case 'T': return 'A';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'N': return 'N';
		case 'a': return 't';
		case 't': return 'a';
		case 'c': return 'g';
		case 'g': return 'c';
		case 'n': return 'n';
		default: throw new IllegalArgumentException();
		}
	}

	public boolean suppressNs = true;

	/**
	 * Returns a fasta record for this entry, in a format that is guaranteed to be compatible with X! Tandem
	 */
	public String toString(){
		String thisSeq = this.getSequence(true);
		if(suppressNs  &&  thisSeq.toUpperCase().contains("N")){
			System.err.println(this.id+" contains Ns, ignoring.");
			return "";
		}else if(thisSeq.length() >= (512*1024)){
			System.err.println(this.id+" sequence too long for X! Tandem, ignoring (length: "+thisSeq.length()+"nt)");
			return "";
			//			return ">"+id+" GN="+this.getAttribute("gene_id")+"\n"+this.getSequence(true);
		}else{
			String geneName = this.getAttribute("gene_id");

			String header = ">"+id;
			if(attributesToAddToFasta != null){
				Iterator<String> it = attributesToAddToFasta.iterator();
				while(it.hasNext())
					header += ":"+this.getAttribute(it.next());
			}
			header += " GN="+geneName;

			//return header+"\n"+this.getSequence(true);
			return header+System.getProperty("line.separator")+this.getSequence(true);
		}
	}


	/**
	 * Check to see if the length of the nucleotide sequence for this entry matches the span of the annotated region from the GTF
	 * 
	 * @return
	 */
	public boolean hasCompleteSequence(){
		return (this.end - this.start + 1) == this.nucleotideSequence.length();
	}



	public boolean overlapsWith(int x1, int x2){
		int y1 = this.start;
		int y2 = this.end;
		return (x1 >= y1 && x1 <= y2) ||
				(x2 >= y1 && x2 <= y2) ||
				(y1 >= x1 && y1 <= x2) ||
				(y2 >= x1 && y2 <= x2);
	}

	public int[] intersectWith(int start, int lineLength){
		return(intersectWith(start, start+lineLength, lineLength));
	}


	public int[] intersectWith(int start, int end, int lineLength){
		int subseqStart, subseqEnd;

		if((subseqStart = this.start - start + 1) < 1){ subseqStart = 1; }
		if((subseqEnd = this.end - start + 1) > lineLength){ subseqEnd = lineLength; }

		return new int[]{subseqStart, subseqEnd};
	}




	public void mergeCoordinates(GenomicCoordinate newCoords){

		if(newCoords.hasAttribute("oId")){	
			if((!this.getAttribute("oId").startsWith("ENST")) && newCoords.getAttribute("oId").startsWith("ENST")){
				this.addAttribute("oId", newCoords.getAttribute("oId"));
			}
		}

		//		if(this.start < newCoords.getStart() && this.end < newCoords.getStart()){
		if(this.end < newCoords.getStop()){
			this.end = newCoords.getStop();
			this.appendSequence(newCoords.getSequence(false));
			//		}else if(newCoords.getStart() < this.start && newCoords.getStop() < this.start){
		}else if(newCoords.getStart() < this.start){
			this.start = newCoords.getStart();
			newCoords.appendSequence(this.nucleotideSequence);
			this.nucleotideSequence = newCoords.getSequence(false);
		}//else{
		// coordinates overlap - can't merge sequences - not a canonical transcript
		//	System.err.println("ERROR: cannot merge ranges:\n\t"+this.toString()+"\n\t"+newCoords.toString()+"\n");
		//}
	}

}
