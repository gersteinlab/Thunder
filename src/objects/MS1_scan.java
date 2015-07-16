package objects;

import java.nio.ByteBuffer;
import java.nio.ByteOrder;
import java.nio.FloatBuffer;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.codec.binary.Base64;
import org.xml.sax.Attributes;

public class MS1_scan {

	/**
	 * 
	 * @param atts
	 */
	public MS1_scan(Attributes atts){
		setAttributes(atts);
	}

	private HashMap<String, String> private_attributes = new HashMap<String, String>();	
	public static final String ATTRIBUTE_SCAN_NUMBER = "num";
	public static final String ATTRIBUTE_MS_LEVEL = "msLevel";
	public static final String ATTRIBUTE_POLARITY = "polarity";
	public static final String ATTRIBUTE_SCAN_TYPE = "scanType";
	public static final String ATTRIBUTE_FILTER_LINE = "filterLine";
	public static final String ATTRIBUTE_RETENTION_TIME = "retentionTime";
	public static final String ATTRIBUTE_LOW_MZ = "lowMz";
	public static final String ATTRIBUTE_HIGH_MZ = "highMz";
	public static final String ATTRIBUTE_BASE_PEAK_MZ = "basePeakMz";
	public static final String ATTRIBUTE_BASE_PEAK_INTENSITY = "basePeakIntensity";
	public static final String ATTRIBUTE_TOTAL_ION_CURRENT = "totIonCurrent";
	public static final String ATTRIBUTE_COLLISION_ENERGY = "collisionEnergy";
	
	/**
	 * 
	 * @param attributes
	 */
	public void setAttributes(Attributes attributes){
		for(int i=0;i<attributes.getLength();i++){
			this.private_attributes.put(attributes.getLocalName(i), attributes.getValue(i));
		}
	}
	
	/**
	 * 
	 * @param attributeID
	 * @return
	 */
	public boolean hasAttribute(String attributeID){
		return this.private_attributes.containsKey(attributeID);
	}
	
	
	/**
	 * 
	 * @param attributeID
	 * @return
	 */
	public String getAttribute(String attributeID){
		try{
			return this.private_attributes.get(attributeID);
		}catch(Exception e){
			return null;
		}
	}
	
	/**
	 * 
	 */
	public void printAttributes(){
		printAttributes("",false);
	}
	
	/**
	 * 
	 * @param prepend
	 * @param newLines
	 */
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
	
	
	
	
	/**
	 * Contains information about the precursor from which this MS/MS scan was obtained
	 */
	private HashMap<String, String> precursor_attributes = new HashMap<String, String>();
	public static final String ATTRIBUTE_PRECURSOR_SCAN_NUMBER = "precursorScanNum";
	public static final String ATTRIBUTE_PRECURSOR_INTENSITY = "precursorIntensity";
	public static final String ATTRIBUTE_PRECURSOR_CHARGE = "precursorCharge";
	public static final String ATTRIBUTE_PRECURSOR_ACTIVATION_METHOD = "activationMethod";
	public static final String ATTRIBUTE_PRECURSOR_MZ = "precursorMz";
	
	public void addPrecursorAttributes(Attributes attributes){
		for(int i=0;i<attributes.getLength();i++)
			this.precursor_attributes.put(attributes.getLocalName(i), attributes.getValue(i));
	}
	public void addPrecursorMZ(String mz){
		this.precursor_attributes.put(ATTRIBUTE_PRECURSOR_MZ, mz);
	}
	public boolean hasPrecursorAttribute(String attributeID){
		return this.precursor_attributes.containsKey(attributeID);
	}
	public String getPrecursorAttribute(String attributeID){
		try{
			return this.precursor_attributes.get(attributeID);
		}catch(Exception e){
			return null;
		}
	}
	
	
	/**
	 * Contains information about the peaks themselves
	 */
	private HashMap<String, String> peaks_attributes = new HashMap<String, String>();
	public static final String ATTRIBUTE_PEAKS_COMPRESSION_TYPE = "compressionType";
	public static final String ATTRIBUTE_PEAKS_COMPRESSED_LENGTH = "compressedLen";
	public static final String ATTRIBUTE_PEAKS_PRECISION = "precision";
	public static final String ATTRIBUTE_PEAKS_BYTE_ORDER = "byteOrder";

	public void addPeaksAttributes(Attributes attributes){
		for(int i=0;i<attributes.getLength();i++)
			this.peaks_attributes.put(attributes.getLocalName(i), attributes.getValue(i));
	}
	
	private HashMap<Float, Float> peaks = new HashMap<Float, Float>();
	public void addPeaks(String peaks){
		Base64 decode = new Base64();
		ByteBuffer bytes = ByteBuffer.wrap(decode.decode(peaks)).order(ByteOrder.BIG_ENDIAN);

		FloatBuffer floatbuffer = bytes.asFloatBuffer();
		float[] floatarray = new float[floatbuffer.remaining()];
		floatbuffer.get(floatarray);
		
		for(int i=0;i<floatarray.length;i+=2){
			this.peaks.put(floatarray[i], floatarray[i+1]);	
		}
	}
	public HashMap<Float, Float> getPeaks(){
		return this.peaks;
	}
	public float getPeakIntensity(float mz){
		return this.peaks.get(mz);
	}
	
	
	public String toString(){
		String toReturn = "MS2 scan number: "+getAttribute(MS1_scan.ATTRIBUTE_SCAN_NUMBER)+
						"\nMS1 charge: "+getPrecursorAttribute(MS1_scan.ATTRIBUTE_PRECURSOR_CHARGE)+
						"\nMS1 m/z: "+getPrecursorAttribute(MS1_scan.ATTRIBUTE_PRECURSOR_MZ)+
						"\nMS1 intensity: "+getPrecursorAttribute(MS1_scan.ATTRIBUTE_PRECURSOR_INTENSITY);
		return toReturn;
	}
	
}
