package fastqTools;

public class SequenceRecord {
	private String id;
	private StringBuffer seq = new StringBuffer();
	private StringBuffer qual = new StringBuffer();
	
	public SequenceRecord(String sequenceID){
		this.id=sequenceID;	
	}
	
	public void addSequenceString(String sequenceString){ this.seq.append(sequenceString); }
	public void addQualityString(String qualString){ this.qual.append(qualString); }
	
	public String getSequence(){ return this.seq.toString().trim(); }
	public String getQuality(){ return this.qual.toString().trim(); }
	
	public int getSequenceLength(){ return this.seq.toString().trim().length(); }
	public String getSequenceID(){ return this.id; }
	public void setSequenceID(String newID){ this.id = newID; }
	
	public void setQuality(String qual){ this.qual = new StringBuffer(qual); }
	
	public String toString(){
		if(this.qual.length() == 0)
			return ">"+this.id+"\n"+this.seq.toString();
		else{
			return "@"+this.id+"\n"+this.seq.toString()+"\n+\n"+this.qual.toString();
		}
	}
}
