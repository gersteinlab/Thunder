package proteome;

import java.util.HashMap;

import objects.MS1_scan;

import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

public class SAXHandler_mzXML extends DefaultHandler{

	HashMap<String, MS1_scan> ms1Intensities = new HashMap<String, MS1_scan>();
	public HashMap<String, MS1_scan> getMS1Intensities(){
		return(ms1Intensities);
	}
	
	private MS1_scan thisScan;
	private boolean parsingPrecursorInfo = false;
	private boolean parsingPeaks = false;

	//private int spectraCount = 0;
	private void updateSpectraCount(){
		/*this.spectraCount ++;
		if(this.spectraCount % 100000 == 0){ System.err.println(this.spectraCount); }
		else if(this.spectraCount % 10000 == 0){ System.err.print(this.spectraCount); }
		else if(this.spectraCount % 1000 == 0){ System.err.print("."); }*/
	}


	/**
	 * What to do when encountering the START of a new XML tag named 'qName'
	 */
	public void startElement(String uri, String localName, String qName, Attributes attributes) throws SAXException {
		if(qName.equals("scan")){
			thisScan = new MS1_scan(attributes);
			updateSpectraCount();
		}else if(qName.equals("precursorMz")){
			thisScan.addPrecursorAttributes(attributes);
			parsingPrecursorInfo = true;
		}else if(qName.equals("peaks")){
			this.thisScan.addPeaksAttributes(attributes);
			parsingPeaks = true;
		}
	}

	/**
	 * What to do when encountering the END of a new XML tag named 'qName'
	 */
	public void endElement(String uri, String localName, String qName) throws SAXException {
		if(qName.equals("scan")){

		}else if(qName.equals("precursorMz")){
			parsingPrecursorInfo = false;
		}else if(qName.equals("peaks")){
			parsingPeaks = false;
			if(thisScan.getAttribute(MS1_scan.ATTRIBUTE_MS_LEVEL).equals("2")){
				try {
					this.ms1Intensities.put(thisScan.getAttribute(MS1_scan.ATTRIBUTE_SCAN_NUMBER), thisScan);
					//this.dataset.addSpectra_mzXML(thisScan);
				} catch (NumberFormatException e) {
					e.printStackTrace();
				}
			}
		}
	}

	/**
	 *  Gets the node's value (only used here to obtain the precursor peak M/Z):
	 */
	public void characters(char ch[], int start, int length) throws SAXException {
		if(parsingPrecursorInfo){
			thisScan.addPrecursorMZ(new String(ch, start, length));
		}else if(parsingPeaks){
			// Not currently needed - removed for speed (saves having to decode the peaks string)
			//thisScan.addPrecursorMZ(new String(ch, start, length));
		}
	}

	
}

