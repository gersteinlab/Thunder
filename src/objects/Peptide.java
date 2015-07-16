package objects;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.List;

import org.xml.sax.Attributes;

//public class Peptide extends Feature_DOM{
//public class Peptide extends Feature_SAX{
public class Peptide{

	private HashMap<String, String> private_attributes = new HashMap<String, String>();
	private String peptideSequence;
	private ArrayList<Spectra> spectra = new ArrayList<Spectra>();


	//	public static final String ATTRIBUTE_ID = "seq";
	public static final String ATTRIBUTE_ID = "id";
	public static final String ATTRIBUTE_B_IONS = "b_ions";
	public static final String ATTRIBUTE_B_SCORE = "b_score";
	public static final String ATTRIBUTE_DELTA = "delta";
	public static final String ATTRIBUTE_END = "end";
	public static final String ATTRIBUTE_EXPECT = "expect";
	public static final String ATTRIBUTE_HYPERSCORE = "hyperscore";
	public static final String ATTRIBUTE_MH = "mh";
	public static final String ATTRIBUTE_MISSED_CLEAVAGES = "missed_cleavages";
	public static final String ATTRIBUTE_NEXTSCORE = "nextscore";
	public static final String ATTRIBUTE_POST = "post";
	public static final String ATTRIBUTE_PRE = "pre";
	public static final String ATTRIBUTE_SEQ = "seq";
	public static final String ATTRIBUTE_START = "start";
	public static final String ATTRIBUTE_Y_IONS = "y_ions";
	public static final String ATTRIBUTE_Y_SCORE = "y_score";



	//	public Peptide(NamedNodeMap attributes){
	//		super(attributes);
	//		this.peptideSequence = attributes.getNamedItem(Peptide.ATTRIBUTE_ID).getNodeValue();  
	//	}
	public Peptide(Attributes peptide_attributes){
		//		super(peptide_attributes);
		//		this.setAttributes(new Attributes(peptide_attributes));

		this.setAttributes(peptide_attributes);

		//		System.out.println("ID: "+getAttribute(Peptide.ATTRIBUTE_ID));

		//		this.peptideSequence = peptide_attributes.getValue(Peptide.ATTRIBUTE_ID);
		this.peptideSequence = peptide_attributes.getValue(Peptide.ATTRIBUTE_SEQ);
	}

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

	public void addSpectra(Spectra spectra){ this.spectra.add(spectra); }
	public String getSequence(){ return this.peptideSequence; }

	private List<HashMap<String, String>> modifications = new ArrayList<HashMap<String, String>>();
	public void addModification(Attributes modification_attributes){
		HashMap<String, String> tmp = new HashMap<String, String>();
		for(int i=0;i<modification_attributes.getLength();i++){	
			tmp.put(modification_attributes.getLocalName(i), modification_attributes.getValue(i));
		}
		this.modifications.add(tmp);
	}
	public List<HashMap<String, String>> getModifications(){
		return this.modifications;
	}
	public String getModificationsString(){
		Iterator<HashMap<String, String>> it = getModifications().iterator();
		String modificationString = "";
		int count = 0;
		while(it.hasNext()){
			if(count > 0)
				modificationString += " | ";
			modificationString += formatModificationString(it.next());
		}
		return modificationString;
	}
	private String formatModificationString(HashMap<String, String> mods){
		String toReturn = "";
		toReturn += mods.get(Peptide.MODIFICATIONS_ATTRIBUTE_TYPE);
		toReturn += ":"+mods.get(Peptide.MODIFICATIONS_ATTRIBUTE_AT);
		toReturn += ":"+mods.get(Peptide.MODIFICATIONS_ATTRIBUTE_MODIFIED);
		return(toReturn);
	}

	public static final String MODIFICATIONS_ATTRIBUTE_TYPE = "type";
	public static final String MODIFICATIONS_ATTRIBUTE_AT = "at";
	public static final String MODIFICATIONS_ATTRIBUTE_MODIFIED = "modified";

}
