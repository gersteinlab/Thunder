package proteome;

import objects.Isoform;
import objects.Peptide;
import objects.Spectra;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import utils.IO_utils;

/**
 * Defines how to deal with each line of XML from the X!Tandem output file
 * 
 */
public class SAXHandler_XTandem extends DefaultHandler{

	private Spectra thisSpectra;
	private Peptide thisPeptide;
	private Isoform thisIsoform;
	//private boolean parsingPeptide = false;
	//private boolean parsingNote = false;

	private int spectraCount = 0;
	private void updateSpectraCount(){
		this.spectraCount ++;
		/*if(this.spectraCount % 100000 == 0){ System.err.println(this.spectraCount); }
		else if(this.spectraCount % 10000 == 0){ System.err.print(this.spectraCount); }
		else if(this.spectraCount % 1000 == 0){ System.err.print("."); }*/
		if(this.spectraCount % 1000 == 0){ IO_utils.printErr("Read "+this.spectraCount+" spectra\r"); }
	}
	
	SpectraAlignmentEngine spectraEngine;
	public void attachSpectraAlignmentEngine(SpectraAlignmentEngine spectraEngine){
		this.spectraEngine = spectraEngine;		
	}


	/**
	 * What to do when encountering the START of a new XML tag named 'qName'
	 */
	public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
		//			System.out.println("Start Element :" + qName);
		if(qName.equals("group")){
			thisSpectra = new Spectra(attributes);
			if(attributes.getValue("type").equals("model")){
				updateSpectraCount();
			}
		}else if(qName.equals("protein")){
			thisIsoform = new Isoform(attributes);
		}else if(qName.equals("note")){
		//	parsingNote = true;
		}else if(qName.equals("peptide")){
			//parsingPeptide = true;
		}else if(qName.equals("domain")){
			//parsingPeptide = false;	
			thisPeptide = new Peptide(attributes);

		}else if(qName.equals("aa")){
			//System.out.print("FOUND AA: "+attributes.getLength()+"\t");
			thisPeptide.addModification(attributes);
			//System.out.println(thisPeptide.getModifications().get(0).getLength());
		}
	}

	/**
	 * What to do when encountering the END of a new XML tag named 'qName'
	 */
	public void endElement(String uri, String localName, String qName) throws SAXException {
		if(qName.equals("protein")){
		}else if(qName.equals("peptide")){
			//parsingPeptide = false;
			try {
				//dataset.addSpectra(thisSpectra, thisPeptide, thisIsoform);
				this.spectraEngine.addSpectraAlignment(thisSpectra, thisPeptide, thisIsoform);
			} catch (NumberFormatException e) {
				e.printStackTrace();
			} catch(Exception e){
				System.err.println(thisPeptide.getAttribute(Peptide.ATTRIBUTE_ID));
				System.err.println(thisIsoform.getAttribute(Isoform.ATTRIBUTE_LABEL));
				e.printStackTrace();
			}
		}else if(qName.equals("note")){
			//parsingNote = false;
		}
	}

	/**
	 *  Gets the node's value (only used here to obtain the peptide sequence):
	 */
	public void characters(char ch[], int start, int length) throws SAXException {
		//if(parsingNote){
		//	thisIsoform.setIsoformID(new String(ch, start, length));
		//}
		//else if(parsingPeptide){
		//	thisIsoform.appendIsoformSequence((new String(ch, start, length)).replaceAll(" |\t|\n", "").trim());
		//}
	}

}


