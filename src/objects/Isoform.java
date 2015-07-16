package objects;

import java.util.HashMap;
import java.util.Iterator;

import org.xml.sax.Attributes;

//public class Isoform extends Feature_DOM{
//public class Isoform extends Feature_SAX{
public class Isoform{
	
	private HashMap<String, String> private_attributes = new HashMap<String, String>();
	private HashMap<String, Peptide> peptides = new HashMap<String, Peptide>(); 
//	private int isoformID;
	private String proteinRefFile;
	private String isoformSequence = "";
	
	
	public static final String ATTRIBUTE_ID = "uid";
	//public static final String ATTRIBUTE_ID = "id";
	public static final String ATTRIBUTE_EXPECT = "expect";
	public static final String ATTRIBUTE_LABEL = "label";
	public static final String ATTRIBUTE_SUM_I = "sumI";

	
//	public Isoform(NamedNodeMap attributes){
//		super(attributes);
//		this.isoformID = Integer.valueOf(attributes.getNamedItem(Isoform.ATTRIBUTE_ID).getNodeValue());  
//	}
	public Isoform(Attributes isoform_attributes){
		this.setAttributes(isoform_attributes);
//		this.isoformID = Integer.valueOf(isoform_attributes.getValue(Isoform.ATTRIBUTE_ID));
	}
	
//	public String getID(){ return this.getAttribute(Isoform.ATTRIBUTE_ID); }
	
	
	
	public void addPeptide(Peptide p){
		if( !this.peptides.containsKey(p.getSequence()) )
				this.peptides.put(p.getSequence(), p); 
	}
	
	
	public int getNumberOfPeptides(){ return this.peptides.size(); }
	public HashMap<String, Peptide> getAllPeptides(){ return this.peptides; }
	public Peptide getPeptide(Peptide peptide){ return this.peptides.get(peptide.getSequence()); }
	public Peptide getPeptide(String peptideID){ return this.peptides.get(this.peptides.get(peptideID)); }
	
	public void setProteinReferenceURL(String proteinRefFile){ this.proteinRefFile = proteinRefFile; }
	public void setIsoformSequence(String proteinSequence){ this.isoformSequence = proteinSequence; }
	public void appendIsoformSequence(String proteinSequence){ this.isoformSequence += proteinSequence; }
	public String getProteinReferenceURL(){ return this.proteinRefFile; }
	public String getIsoformSequence(){ return this.isoformSequence; }
	public String getIsoformID(){ return this.getAttribute(ATTRIBUTE_LABEL); }
	
	
	
	public void setAttributes(Attributes attributes){
		for(int i=0;i<attributes.getLength();i++){
			this.private_attributes.put(attributes.getLocalName(i), attributes.getValue(i));
		}
	}
	public String getAttribute(String attributeID){
		try{
			return this.private_attributes.get(attributeID);
		}catch(Exception e){
			return null;
		}
	}
	public void printAttributes(){
		printAttributes("",false);
	}
	public void printAttributes(String prepend, boolean newLines){
		Iterator<String> it = this.private_attributes.keySet().iterator();
		String tmp;
		while(it.hasNext()){
			tmp = it.next();
			if(newLines){
				System.out.println(prepend+" \'"+tmp+"\':\t"+this.private_attributes.get(tmp));
			}else{
				System.out.print("\t\'"+tmp+"\'="+this.private_attributes.get(tmp));
			}
		}
	}
	
}
