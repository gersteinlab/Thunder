package objects;

import java.util.List;

import net.sf.samtools.Cigar;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecord.SAMTagAndValue;

public class SAMRecordReduced {

	private int _samFlag = -1;
	private int _alignmentStart, _alignmentEnd;
	private String _readName, _referenceName, _readString;
	private boolean _readIsNegativeStrand;
	private Cigar _CIGAR = null;
	private List<SAMTagAndValue> _SAMAttributes;
	
	
	/**
	 * A reduced memory SAMRecord object (keeps only alignment start/stop, reference name, read name, read sequence, strand, CIGAR)
	 * @param originalSAMRecord
	 */
	public SAMRecordReduced(SAMRecord originalSAMRecord){
		_alignmentStart = originalSAMRecord.getAlignmentStart();
		_alignmentEnd = originalSAMRecord.getAlignmentEnd();
		_referenceName = originalSAMRecord.getReferenceName();
		_readIsNegativeStrand = originalSAMRecord.getReadNegativeStrandFlag();
		_readName = originalSAMRecord.getReadName();
		_readString = originalSAMRecord.getReadString();
		_CIGAR = originalSAMRecord.getCigar();
		_samFlag = originalSAMRecord.getFlags();
	}
	
	public boolean isPrimaryAlignment(){
		//if(_samFlag < 0)
		//	return false;
		//else 
		if(_samFlag == 0  ||  _samFlag == 16)
			return true;
		else
			return false;
	}
	
	/**
	 * A reduced memory SAMRecord object
	 * recommended (for lower memory footprint): new SAMRecordReduced(originalRecord, true, true, false, false);
	 * 
	 * @param originalSAMRecord
	 * @param keepReadName
	 * @param keepCIGAR
	 * @param keepReadSequence
	 * @param keepSAMAttributes
	 */
	public SAMRecordReduced(SAMRecord originalSAMRecord, boolean keepReadName, boolean keepCIGAR, boolean keepReadSequence, boolean keepSAMAttributes){
		_alignmentStart = originalSAMRecord.getAlignmentStart();
		_alignmentEnd = originalSAMRecord.getAlignmentEnd();
		_referenceName = originalSAMRecord.getReferenceName();
		_readIsNegativeStrand = originalSAMRecord.getReadNegativeStrandFlag();
		
		if(keepReadName)
			_readName = originalSAMRecord.getReadName();
		
		if(keepCIGAR)
			_CIGAR = originalSAMRecord.getCigar();
		
		if(keepReadSequence)
			_readString = originalSAMRecord.getReadString();
		
		if(keepSAMAttributes)
			_SAMAttributes = originalSAMRecord.getAttributes();
	}
	
	
	public int getAlignmentStart(){ return _alignmentStart; }
	public int getAlignmentEnd(){ return _alignmentEnd; }
	public String getReferenceName(){ return _referenceName; }
	public String getReadName(){ return _readName; }
	public boolean getReadNegativeStrandFlag(){ return _readIsNegativeStrand; }
	public Cigar getCigar(){ return _CIGAR; }
	public String getReadSequence(){ return _readString; }
	public List<SAMTagAndValue> getAttributes(){ return _SAMAttributes; }
	
}
