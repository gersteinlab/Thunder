package objects;

public class FastX_Record {

	private String id, seq, qual;
	public FastX_Record(String id, String seq, String qual){
		this.id = id;
		this.seq = seq;
		this.qual = qual;
	}
	public FastX_Record(String id, String seq){
		this.id = id;
		this.seq = seq;
		this.qual = null;
	}
	
	public String getID(){ return this.id; }
	public String getSequence(){ return this.seq; }
	public String getQuality(){ return this.qual; }
	
	public String toString(){
		if(this.qual == null){
			return ">"+this.id+"\n"+this.seq;
		}else{
			return "@"+this.id+"\n"+this.seq+"\n+\n"+this.qual;
		}
	}
	
}
