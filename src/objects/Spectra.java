package objects;

import java.util.HashMap;
import java.util.Iterator;

import org.xml.sax.Attributes;

//public class Spectra extends Feature_SAX {
public class Spectra {
	
	private HashMap<String, String> private_attributes = new HashMap<String, String>();
//	private Attributes private_attributes;
//	private String spectraID;
	
	public static final String ATTRIBUTE_ID = "id";
	public static final String ATTRIBUTE_ACT = "act";
	public static final String ATTRIBUTE_EXPECT = "expect";
	public static final String ATTRIBUTE_F_I = "fI";
	public static final String ATTRIBUTE_LABEL = "label";
	public static final String ATTRIBUTE_MAX_I = "maxI";
	public static final String ATTRIBUTE_MH = "mh";
	public static final String ATTRIBUTE_RT = "rt";
	public static final String ATTRIBUTE_SUM_I = "sumI";
	public static final String ATTRIBUTE_TYPE = "type";
	public static final String ATTRIBUTE_Z = "z";
	
	
	
	
	
	
	public void setAttributes(Attributes attributes){
		for(int i=0;i<attributes.getLength();i++){
			this.private_attributes.put(attributes.getLocalName(i), attributes.getValue(i));
		}
	}
	
	public Spectra(Attributes spectra_attributes){
//		super(spectra_attributes);
		this.setAttributes(spectra_attributes);
//		this.spectraID = spectra_attributes.getValue(Spectra.ATTRIBUTE_ID);
	}

//	public String getID(){ return this.spectraID; }

	
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
