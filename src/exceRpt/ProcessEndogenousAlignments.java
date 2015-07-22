package exceRpt;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;
import java.util.Map.Entry;

import main.Thunder;
import net.sf.samtools.Cigar;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileReader.ValidationStringency;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.apache.commons.cli.ParseException;

public class ProcessEndogenousAlignments {

	private BufferedWriter _optimalAlignmentWriter;
	HashMap <String, Hairpin> miRNAprecursors = new HashMap <String, Hairpin>();


	public ProcessEndogenousAlignments(String pathForOptimalAlignments) throws IOException{
		_optimalAlignmentWriter = new BufferedWriter(new FileWriter(pathForOptimalAlignments+"/endogenousAlignments_Accepted.txt"));
	}


	/*
	 * 
	 */
	public void tidyUp() throws IOException{
		_optimalAlignmentWriter.flush();
		_optimalAlignmentWriter.close();
	}



	/**
	 * 
	 * @param thisRead
	 */
	public void assignRead(HashMap<SAMRecord, String> thisRead){
		Iterator<Entry<SAMRecord,String>> it = thisRead.entrySet().iterator();
		Entry<SAMRecord,String> tmp;

		// store the alignment records and the libraries to which they map
		HashMap<String, ArrayList<SAMRecord>> readsByLibrary = new HashMap<String, ArrayList<SAMRecord>>();  
		ArrayList<SAMRecord> keepAlignments = new ArrayList<SAMRecord>();
		//SAMRecord[] multimaps = new SAMRecord[thisRead.size()];
		//String[] libraries = new String[thisRead.size()];

		// set up counters for each library
		//HashMap<String, Integer> counters = new HashMap<String, Integer>(); 
		boolean readHasSenseAlignment = false;


		/*
		 * Loop through alignments and count libraries that capture alignments
		 */
		while(it.hasNext()){
			tmp = it.next();
			boolean readNegativeStrand = tmp.getKey().getReadNegativeStrandFlag();
			if(!readNegativeStrand  &&  !readHasSenseAlignment)
				readHasSenseAlignment = true;

			if(tmp.getValue().equals("miRNA")){
				if(!readNegativeStrand)
					readsByLibrary = addReadAssignment(readsByLibrary, "miRNA_sense", tmp.getKey());
				else
					readsByLibrary = addReadAssignment(readsByLibrary, "miRNA_antisense", tmp.getKey());
			}else if(tmp.getValue().equals("tRNA")){
				if(!readNegativeStrand)
					readsByLibrary = addReadAssignment(readsByLibrary, "tRNA_sense", tmp.getKey());
				else
					readsByLibrary = addReadAssignment(readsByLibrary, "tRNA_antisense", tmp.getKey());
			}else if(tmp.getValue().equals("piRNA")){
				if(!readNegativeStrand)
					readsByLibrary = addReadAssignment(readsByLibrary, "piRNA_sense", tmp.getKey());
				else
					readsByLibrary = addReadAssignment(readsByLibrary, "piRNA_antisense", tmp.getKey());
			}else if(tmp.getValue().equals("circRNA")){
				if(!readNegativeStrand)
					readsByLibrary = addReadAssignment(readsByLibrary, "circRNA_sense", tmp.getKey());
				else
					readsByLibrary = addReadAssignment(readsByLibrary, "circRNA_antisense", tmp.getKey());
			}else if(tmp.getValue().equals("gencode")){
				if(!readNegativeStrand)
					readsByLibrary = addReadAssignment(readsByLibrary, "gencode_sense", tmp.getKey());
				else
					readsByLibrary = addReadAssignment(readsByLibrary, "gencode_antisense", tmp.getKey());
			}

			//System.out.println(tmp.getKey().getReadName()+"\t"+tmp.getKey().getReferenceName()+"\t"+tmp.getValue()+"\t"+(!readNegativeStrand));
		}


		/*
		 * Logic block to prioritise alignments:
		 * 
		 * 1- prefer sense alignments over antisense alignments
		 * 2- miRNA
		 * 3- tRNA
		 * 4- piRNA
		 * 5- gencode
		 * 6- circRNA
		 */
		String keptLibrary = "none";
		if(readHasSenseAlignment){
			if(readsByLibrary.containsKey("miRNA_sense")){
				keepAlignments = readsByLibrary.get("miRNA_sense");
				keptLibrary = "miRNA_sense";
			}else if(readsByLibrary.containsKey("tRNA_sense")){
				keepAlignments = readsByLibrary.get("tRNA_sense");
				keptLibrary = "tRNA_sense";
			}else if(readsByLibrary.containsKey("piRNA_sense")){
				keepAlignments = readsByLibrary.get("piRNA_sense");
				keptLibrary = "piRNA_sense";
			}else if(readsByLibrary.containsKey("gencode_sense")){
				keepAlignments = readsByLibrary.get("gencode_sense");
				keptLibrary = "gencode_sense";
			}else if(readsByLibrary.containsKey("circRNA_sense")){
				keepAlignments = readsByLibrary.get("circRNA_sense");
				keptLibrary = "circRNA_sense";
			}
		}else{
			if(readsByLibrary.containsKey("miRNA_antisense")){
				keepAlignments = readsByLibrary.get("miRNA_antisense");
				keptLibrary = "miRNA_antisense";
			}else if(readsByLibrary.containsKey("tRNA_antisense")){
				keepAlignments = readsByLibrary.get("tRNA_antisense");
				keptLibrary = "tRNA_antisense";
			}else if(readsByLibrary.containsKey("piRNA_antisense")){
				keepAlignments = readsByLibrary.get("piRNA_antisense");
				keptLibrary = "piRNA_antisense";
			}else if(readsByLibrary.containsKey("gencode_antisense")){
				keepAlignments = readsByLibrary.get("gencode_antisense");
				keptLibrary = "gencode_antisense";
			}else if(readsByLibrary.containsKey("circRNA_antisense")){
				keepAlignments = readsByLibrary.get("circRNA_antisense");
				keptLibrary = "circRNA_antisense";
			}
		}


		/*
		 * Write 'accepted' alignments and proportionally assign reads
		 */
		int nMultimaps = keepAlignments.size();
		ArrayList<String> matureAlignments = new ArrayList<String>(); 
		//System.out.println("nMaps = "+nMultimaps);
		Iterator<SAMRecord> it2 = keepAlignments.iterator();
		SAMRecord tmp2;
		while(it2.hasNext()){
			
			tmp2 = it2.next();
			//System.out.println(tmp2.getReferenceName());
			String[] refIDbits = tmp2.getReferenceName().split(":");
			String isGenomeMapped = refIDbits[0];
			String mapsTo = refIDbits[2];
			for(int i=3;i<refIDbits.length;i++)
				mapsTo = mapsTo.concat(":"+refIDbits[i]);

			if(keptLibrary.equals("miRNA_sense")){
				miRNAprecursors.get(mapsTo).addRead(tmp2.getReadName(), tmp2.getAlignmentStart(), tmp2.getAlignmentEnd(), false);
			}else if(keptLibrary.equals("miRNA_antisense")){
				miRNAprecursors.get(mapsTo).addRead(tmp2.getReadName(), tmp2.getAlignmentStart(), tmp2.getAlignmentEnd(), true);
			}
			else if(keptLibrary.equals("tRNA_sense")){
				if(!readCounts_tRNA_sense.containsKey(mapsTo))
					readCounts_tRNA_sense.put(mapsTo, 0.0);
				readCounts_tRNA_sense.put(mapsTo, readCounts_tRNA_sense.get(mapsTo) + 1.0/(nMultimaps+0.0));
			}else if(keptLibrary.equals("tRNA_antisense")){
				if(!readCounts_tRNA_antisense.containsKey(mapsTo))
					readCounts_tRNA_antisense.put(mapsTo, 0.0);
				readCounts_tRNA_antisense.put(mapsTo, readCounts_tRNA_antisense.get(mapsTo) + 1.0/(nMultimaps+0.0));
			}
			else if(keptLibrary.equals("piRNA_sense")){
				if(!readCounts_piRNA_sense.containsKey(mapsTo))
					readCounts_piRNA_sense.put(mapsTo, 0.0);
				readCounts_piRNA_sense.put(mapsTo, readCounts_piRNA_sense.get(mapsTo) + 1.0/(nMultimaps+0.0));
			}else if(keptLibrary.equals("piRNA_antisense")){
				if(!readCounts_piRNA_antisense.containsKey(mapsTo))
					readCounts_piRNA_antisense.put(mapsTo, 0.0);
				readCounts_piRNA_antisense.put(mapsTo, readCounts_piRNA_antisense.get(mapsTo) + 1.0/(nMultimaps+0.0));
			}
			else if(keptLibrary.equals("circRNA_sense")){
				if(!readCounts_circRNA_sense.containsKey(mapsTo))
					readCounts_circRNA_sense.put(mapsTo, 0.0);
				readCounts_circRNA_sense.put(mapsTo, readCounts_circRNA_sense.get(mapsTo) + 1.0/(nMultimaps+0.0));
			}else if(keptLibrary.equals("circRNA_antisense")){
				if(!readCounts_circRNA_antisense.containsKey(mapsTo))
					readCounts_circRNA_antisense.put(mapsTo, 0.0);
				readCounts_circRNA_antisense.put(mapsTo, readCounts_circRNA_antisense.get(mapsTo) + 1.0/(nMultimaps+0.0));
			}
			else if(keptLibrary.equals("gencode_sense")){
				if(!readCounts_gencode_sense.containsKey(mapsTo))
					readCounts_gencode_sense.put(mapsTo, 0.0);
				readCounts_gencode_sense.put(mapsTo, readCounts_gencode_sense.get(mapsTo) + 1.0/(nMultimaps+0.0));
			}else if(keptLibrary.equals("gencode_antisense")){
				if(!readCounts_gencode_antisense.containsKey(mapsTo))
					readCounts_gencode_antisense.put(mapsTo, 0.0);
				readCounts_gencode_antisense.put(mapsTo, readCounts_gencode_antisense.get(mapsTo) + 1.0/(nMultimaps+0.0));
			}
			
			
			
			/*
			 * Write chosen alignments
			 */
			//System.out.print(tmp2.getReadName()+"\t"+mapsTo+"\t"+keptLibrary+"\t"+tmp2.getAlignmentStart()+"\t"+tmp2.getAlignmentEnd()+"\tNM:i:"+tmp2.getIntegerAttribute("NM")+"\tMD:Z:"+tmp2.getStringAttribute("MD"));
			try {
				_optimalAlignmentWriter.write(tmp2.getReadName()+"\t"+keptLibrary+"\t"+isGenomeMapped+"\t"+mapsTo+"\t"+tmp2.getAlignmentStart()+"\t"+tmp2.getAlignmentEnd()+"\tNM:i:"+tmp2.getIntegerAttribute("NM")+"\tMD:Z:"+tmp2.getStringAttribute("MD"));
			} catch (IOException e) {
				// TODO Auto-generated catch block
				e.printStackTrace();
			}

			
			
			/*
			 * If the read aligns to a precursor miRNA, try to find a mature sequences that it overlaps
			 */
			if(keptLibrary.equals("miRNA_sense")  ||  keptLibrary.equals("miRNA_antisense")){
				//miRNAprecursors.get(mapsTo).getMatureOverlapForRead(tmp2.getReadName());
				String maps2mature = miRNAprecursors.get(mapsTo).getMatureOverlapForRead(tmp2.getReadName());
				//System.out.println("\t"+maps2mature);
				try { _optimalAlignmentWriter.write("\t"+maps2mature+"\n");
				} catch (IOException e) { e.printStackTrace(); }
				if(maps2mature.length() > 0){
					if(!matureAlignments.contains(maps2mature))
						matureAlignments.add(maps2mature);
				}else{
					// read does not hit annotated mature miRNA, assign to hairpin instead
					if(keptLibrary.equals("miRNA_sense")){
						if(!readCounts_precursorMiRNA_sense.containsKey(mapsTo))
							readCounts_precursorMiRNA_sense.put(mapsTo, 0.0);
						readCounts_precursorMiRNA_sense.put(mapsTo, readCounts_precursorMiRNA_sense.get(mapsTo)+(1.0/(nMultimaps+0.0)));
					}else{
						if(!readCounts_precursorMiRNA_antisense.containsKey(mapsTo))
							readCounts_precursorMiRNA_antisense.put(mapsTo, 0.0);
						readCounts_precursorMiRNA_antisense.put(mapsTo, readCounts_precursorMiRNA_antisense.get(mapsTo)+(1.0/(nMultimaps+0.0)));
					}
				}
			}else{
				try { _optimalAlignmentWriter.write("\tNA\n");
				} catch (IOException e) { e.printStackTrace(); }
				//System.out.println("\tNA");
			}
		}


		/*
		 * For mature miRNA alignments, collapse replicate mature IDs (from mature sequences that align to multiple precursors) 
		 */
		if(keptLibrary.equals("miRNA_sense")  ||  keptLibrary.equals("miRNA_antisense")){
			//System.out.println(matureAlignments.size());
			Iterator<String> it3 = matureAlignments.iterator();
			String tmpID;
			while(it3.hasNext()){
				tmpID = it3.next();
				if(keptLibrary.equals("miRNA_sense")){
					if(!readCounts_matureMiRNA_sense.containsKey(tmpID))
						readCounts_matureMiRNA_sense.put(tmpID, 0.0);
					readCounts_matureMiRNA_sense.put(tmpID, readCounts_matureMiRNA_sense.get(tmpID)+(1.0/(matureAlignments.size()+0.0)));
				}else{
					if(!readCounts_matureMiRNA_antisense.containsKey(tmpID))
						readCounts_matureMiRNA_antisense.put(tmpID, 0.0);
					readCounts_matureMiRNA_antisense.put(tmpID, readCounts_matureMiRNA_antisense.get(tmpID)+(1.0/(matureAlignments.size()+0.0)));
				}
			}
		}else{

		}
	}
	
	
	

	private HashMap<String, Double> readCounts_matureMiRNA_sense = new HashMap<String, Double>();
	private HashMap<String, Double> readCounts_precursorMiRNA_sense = new HashMap<String, Double>();
	private HashMap<String, Double> readCounts_tRNA_sense = new HashMap<String, Double>();
	private HashMap<String, Double> readCounts_piRNA_sense = new HashMap<String, Double>();
	private HashMap<String, Double> readCounts_circRNA_sense = new HashMap<String, Double>();
	private HashMap<String, Double> readCounts_gencode_sense = new HashMap<String, Double>();
	
	private HashMap<String, Double> readCounts_matureMiRNA_antisense = new HashMap<String, Double>();
	private HashMap<String, Double> readCounts_precursorMiRNA_antisense = new HashMap<String, Double>();
	private HashMap<String, Double> readCounts_tRNA_antisense = new HashMap<String, Double>();
	private HashMap<String, Double> readCounts_piRNA_antisense = new HashMap<String, Double>();
	private HashMap<String, Double> readCounts_circRNA_antisense = new HashMap<String, Double>();
	private HashMap<String, Double> readCounts_gencode_antisense = new HashMap<String, Double>();

	

	/**
	 * @throws IOException 
	 * 
	 */
	public void writeCounts_miRNA(String basePath) throws IOException{
		//
		// Do the mature miRNAs
		//
		BufferedWriter out = new BufferedWriter(new FileWriter(basePath+"/mature_sense.tmp"));
		Iterator<String> it = readCounts_matureMiRNA_sense.keySet().iterator();
		String tmpID = "";
		while(it.hasNext()){
			tmpID = it.next();
			out.write(tmpID+"\t"+readCounts_matureMiRNA_sense.get(tmpID)+"\n");
		}
		out.flush();
		out.close();
		out = new BufferedWriter(new FileWriter(basePath+"/mature_antisense.tmp"));
		it = readCounts_matureMiRNA_antisense.keySet().iterator();
		while(it.hasNext()){
			tmpID = it.next();
			out.write(tmpID+"\t"+readCounts_matureMiRNA_antisense.get(tmpID)+"\n");
		}
		out.flush();
		out.close();

		//
		// Do the precursor miRNAs
		//
		out = new BufferedWriter(new FileWriter(basePath+"/hairpin_sense.tmp"));
		it = readCounts_precursorMiRNA_sense.keySet().iterator();
		while(it.hasNext()){
			tmpID = it.next();
			out.write(tmpID+"\t"+readCounts_precursorMiRNA_sense.get(tmpID)+"\n");
		}
		out.flush();
		out.close();
		out = new BufferedWriter(new FileWriter(basePath+"/hairpin_antisense.tmp"));
		it = readCounts_precursorMiRNA_antisense.keySet().iterator();
		while(it.hasNext()){
			tmpID = it.next();
			out.write(tmpID+"\t"+readCounts_precursorMiRNA_antisense.get(tmpID)+"\n");
		}
		out.flush();
		out.close();
	}


	
	/**
	 * 
	 * @param basePath
	 * @throws IOException
	 */
	public void writeCounts_otherLibraries(String basePath) throws IOException{
		// Do the tRNAs
		if(readCounts_tRNA_sense.size() > 0){
			BufferedWriter out = new BufferedWriter(new FileWriter(basePath+"/tRNA_sense.tmp"));
			Iterator<String> it = readCounts_tRNA_sense.keySet().iterator();
			String tmpID = "";
			while(it.hasNext()){
				tmpID = it.next();
				out.write(tmpID+"\t"+readCounts_tRNA_sense.get(tmpID)+"\n");
			}
			out.flush();
			out.close();
		}
		if(readCounts_tRNA_antisense.size() > 0){
			BufferedWriter out = new BufferedWriter(new FileWriter(basePath+"/tRNA_antisense.tmp"));
			Iterator<String> it = readCounts_tRNA_antisense.keySet().iterator();
			String tmpID = "";
			while(it.hasNext()){
				tmpID = it.next();
				out.write(tmpID+"\t"+readCounts_tRNA_antisense.get(tmpID)+"\n");
			}
			out.flush();
			out.close();
		}

		// Do the piRNAs
		if(readCounts_piRNA_sense.size() > 0){
			BufferedWriter out = new BufferedWriter(new FileWriter(basePath+"/piRNA_sense.tmp"));
			Iterator<String> it = readCounts_piRNA_sense.keySet().iterator();
			String tmpID = "";
			while(it.hasNext()){
				tmpID = it.next();
				//System.out.println(tmpID+"\t"+readCounts_piRNA_sense.get(tmpID));
				out.write(tmpID+"\t"+readCounts_piRNA_sense.get(tmpID)+"\n");
			}
			out.flush();
			out.close();
		}
		if(readCounts_piRNA_antisense.size() > 0){
			BufferedWriter out = new BufferedWriter(new FileWriter(basePath+"/piRNA_antisense.tmp"));
			Iterator<String> it = readCounts_piRNA_antisense.keySet().iterator();
			String tmpID = "";
			while(it.hasNext()){
				tmpID = it.next();
				//System.out.println(tmpID+"\t"+readCounts_piRNA_sense.get(tmpID));
				out.write(tmpID+"\t"+readCounts_piRNA_antisense.get(tmpID)+"\n");
			}
			out.flush();
			out.close();
		}

		// Do the circularRNAs
		if(readCounts_circRNA_sense.size() > 0){
			BufferedWriter out = new BufferedWriter(new FileWriter(basePath+"/circularRNA_sense.tmp"));
			Iterator<String> it = readCounts_circRNA_sense.keySet().iterator();
			String tmpID = "";
			while(it.hasNext()){
				tmpID = it.next();
				//System.out.println(tmpID+"\t"+readCounts_piRNA_sense.get(tmpID));
				out.write(tmpID+"\t"+readCounts_circRNA_sense.get(tmpID)+"\n");
			}
			out.flush();
			out.close();
		}
		if(readCounts_circRNA_antisense.size() > 0){
			BufferedWriter out = new BufferedWriter(new FileWriter(basePath+"/circularRNA_antisense.tmp"));
			Iterator<String> it = readCounts_circRNA_antisense.keySet().iterator();
			String tmpID = "";
			while(it.hasNext()){
				tmpID = it.next();
				//System.out.println(tmpID+"\t"+readCounts_piRNA_sense.get(tmpID));
				out.write(tmpID+"\t"+readCounts_circRNA_antisense.get(tmpID)+"\n");
			}
			out.flush();
			out.close();
		}
		
		// Do the gencode alignments
		if(readCounts_gencode_sense.size() > 0){
			BufferedWriter out = new BufferedWriter(new FileWriter(basePath+"/gencode_sense.tmp"));
			Iterator<String> it = readCounts_gencode_sense.keySet().iterator();
			String tmpID = "";
			while(it.hasNext()){
				tmpID = it.next();
				//System.out.println(tmpID+"\t"+readCounts_piRNA_sense.get(tmpID));
				out.write(tmpID+"\t"+readCounts_gencode_sense.get(tmpID)+"\n");
			}
			out.flush();
			out.close();
		}
		if(readCounts_gencode_antisense.size() > 0){
			BufferedWriter out = new BufferedWriter(new FileWriter(basePath+"/gencode_antisense.tmp"));
			Iterator<String> it = readCounts_gencode_antisense.keySet().iterator();
			String tmpID = "";
			while(it.hasNext()){
				tmpID = it.next();
				//System.out.println(tmpID+"\t"+readCounts_piRNA_sense.get(tmpID));
				out.write(tmpID+"\t"+readCounts_gencode_antisense.get(tmpID)+"\n");
			}
			out.flush();
			out.close();
		}
	}

	

	/**
	 * 
	 * @param readsByLibrary
	 * @param refType
	 * @param alignment
	 * @return
	 */
	private HashMap<String, ArrayList<SAMRecord>> addReadAssignment(HashMap<String, ArrayList<SAMRecord>> readsByLibrary, String refType, SAMRecord alignment){
		if(!readsByLibrary.containsKey(refType))
			readsByLibrary.put(refType, new ArrayList<SAMRecord>());

		readsByLibrary.get(refType).add(alignment);
		return readsByLibrary;
	}



	/**
	 * 
	 * @param path_readAlignments
	 */
	public void read_Reads(File path_readAlignments){
		/*
		 * Read the mature sequences
		 */
		SAMFileReader inputSam = new SAMFileReader(path_readAlignments);
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		//inputSam.setValidationStringency(ValidationStringency.LENIENT);
		SAMRecord thisRecord;
		HashMap<SAMRecord, String> thisRead = new HashMap<SAMRecord, String>(); 

		SAMRecordIterator it = inputSam.iterator();
		String lastReadID = null;
		while(it.hasNext()){
			thisRecord = it.next();
			//System.out.println(thisRecord.getReferenceName());
			if(!thisRecord.getReadName().equals(lastReadID)  &&  lastReadID != null){
				// new, non first
				assignRead(thisRead);
				thisRead = new HashMap<SAMRecord, String>();
			}
			thisRead.put(thisRecord, thisRecord.getReferenceName().split(":")[1]);
			lastReadID = thisRecord.getReadName();
			//System.out.println(thisRecord.getReferenceName().split(":")[0]+"\t"+thisRecord.getReferenceName()+"\t"+thisRecord.getReadName());

		}
		// assign the final read!
		assignRead(thisRead);
		
		inputSam.close();
	}


	/**
	 * 
	 * @param path_hairpin2genome
	 * @param path_mature2hairpin
	 */
	public void read_miRNAinfo(File path_hairpin2genome, File path_mature2hairpin){
		/*
		 * Read the hairpin alignments to the genome
		 */
		SAMFileReader inputSam = new SAMFileReader(path_hairpin2genome);
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		SAMRecord thisRecord;
		SAMRecordIterator it = inputSam.iterator();
		while(it.hasNext()){
			thisRecord = it.next();
			//System.out.println(thisRecord.getReferenceName()+"\t"+thisRecord.getReadName()+"\t"+thisRecord.getReadString());
			if(!miRNAprecursors.containsKey(thisRecord.getReadName()))
				miRNAprecursors.put(thisRecord.getReadName(), new Hairpin(thisRecord.getReadName()));

			if(!thisRecord.getReadUnmappedFlag())
				miRNAprecursors.get(thisRecord.getReadName()).addGenomicAlignment(new Alignment(thisRecord.getReferenceName(), thisRecord.getAlignmentStart(), thisRecord.getAlignmentEnd(), thisRecord.getReadNegativeStrandFlag(), thisRecord.getCigar(), thisRecord.getReadString()));
		}
		inputSam.close();

		/*Iterator<String> it2 = miRNAprecursors.keySet().iterator();
		while(it2.hasNext()){
			Hairpin tmp = miRNAprecursors.get(it2.next());
			System.out.println(tmp.getID()+"\t"+tmp.getNumberOfGenomicAlignments());
		}*/


		/*
		 * Read the mature sequences
		 */
		inputSam = new SAMFileReader(path_mature2hairpin);
		inputSam.setValidationStringency(ValidationStringency.SILENT);
		it = inputSam.iterator();
		while(it.hasNext()){
			thisRecord = it.next();

			//System.out.println(thisRecord.getReferenceName()+"\t"+thisRecord.getReadName());

			if(!thisRecord.getReadUnmappedFlag()  &&  miRNAprecursors.containsKey(thisRecord.getReferenceName()))
				miRNAprecursors.get(thisRecord.getReferenceName()).addMatureMiRNA(new Alignment(thisRecord.getReadName(), thisRecord.getAlignmentStart(), thisRecord.getAlignmentEnd(), thisRecord.getReadNegativeStrandFlag(), thisRecord.getCigar(), thisRecord.getReadString()));


		}
		inputSam.close();
	}

	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("SAM/BAM").hasArg().withDescription("Path to HAIRPIN alignments to the genome").create("hairpin2genome"));
		options.addOption(OptionBuilder.withArgName("SAM/BAM").hasArg().withDescription("Path to MATURE alignments to the hairpins").create("mature2hairpin"));
		options.addOption(OptionBuilder.withArgName("SAM/BAM").hasArg().withDescription("Path to READ alignments to the hairpins").create("reads2all"));
		options.addOption(OptionBuilder.withArgName("directory").hasArg().withDescription("Base path to write the results into").create("outputPath"));
		return options;
	}

	public static void main(String[] args) throws ParseException, IOException {
		/*String hairpin2genome = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/MICRO_RNA/miRBase_v21_hairpin_hsa_hg19_aligned.sam";
		String mature2hairpin = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/MICRO_RNA/miRBase_v21_mature_hairpin_hsa_aligned.sam";
		String reads_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/endogenousAlignments_ALL.sam";
		//String reads_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants/test.sam";
		String output_path = "/Users/robk/WORK/YALE_offline/exRNA/TESTING/newEndogenousQuants";
		args = new String[]{"ProcessEndogenousAlignments",
				"--hairpin2genome",hairpin2genome,
				"--mature2hairpin",mature2hairpin,
				"--reads2all",reads_path,
				"--outputPath",output_path
		};*/

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption("hairpin2genome") && cmdArgs.hasOption("mature2hairpin") && cmdArgs.hasOption("reads2all") && cmdArgs.hasOption("outputPath")){
			//ProcessEndogenousAlignments engine = new ProcessEndogenousAlignments(new File(cmdArgs.getOptionValue("hairpin2genome")), new File(cmdArgs.getOptionValue("mature2hairpin")));
			ProcessEndogenousAlignments engine = new ProcessEndogenousAlignments(cmdArgs.getOptionValue("outputPath"));

			// Read hairpin alignments to the genome and mature alignments to the hairpins
			Thunder.printLineErr("Reading miRNA annotation info");
			engine.read_miRNAinfo(new File(cmdArgs.getOptionValue("hairpin2genome")), new File(cmdArgs.getOptionValue("mature2hairpin")));

			// read and process the RNA-seq read alignments
			Thunder.printLineErr("Processing RNA-seq alignments");
			engine.read_Reads(new File(cmdArgs.getOptionValue("reads2all")));

			// write the miRNA counts
			Thunder.printLineErr("Writing read counts and tidying up");
			engine.writeCounts_miRNA(cmdArgs.getOptionValue("outputPath"));
			// write the counts from other libraries
			engine.writeCounts_otherLibraries(cmdArgs.getOptionValue("outputPath"));
			
			// close global buffered writer(s)
			engine.tidyUp();
			
			Thunder.printLineErr("Done!");

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ProcessEndogenousAlignments", getCmdLineOptions());
			System.err.println();
		}



	}

}


class Hairpin{
	private String _id; 
	public Hairpin(String id){
		_id = id;
	}

	public ArrayList<Alignment> _genomicAlignments = new ArrayList<Alignment>();
	public HashMap<String, Alignment> _matureMiRNAs = new HashMap<String, Alignment>();
	public HashMap<String, Alignment> _reads = new HashMap<String, Alignment>();

	//public void addGenomicAlignment(String chr, int start, int stop, boolean isNegativeStrand){
	public void addGenomicAlignment(Alignment newAlignment){ _genomicAlignments.add(newAlignment); }
	public void addMatureMiRNA(Alignment newAlignment){
		if(_matureMiRNAs.containsKey(newAlignment.getReferenceID())){
			if(_matureMiRNAs.get(newAlignment.getReferenceID()).isNegativeStrand()  &&  !newAlignment.isNegativeStrand())
				_matureMiRNAs.put(newAlignment.getReferenceID(), newAlignment); //if the existing alignment is antisense to the hairpin, replace it with a sense alignment				
		}else{
			_matureMiRNAs.put(newAlignment.getReferenceID(), newAlignment);
		}
		 
	}
	public void addRead(String readID, int start, int stop, boolean isNegativeStrand){ 
		_reads.put(readID, new Alignment(_id, start, stop, isNegativeStrand)); 
	}


	public String getID(){ return _id; }
	public int getNumberOfGenomicAlignments(){ return _genomicAlignments.size(); }
	public ArrayList<Alignment> getGenomicAlignments(){ return _genomicAlignments; }


	public String getMatureOverlapForRead(String readID){
		double minFractionForOverlap = 0.8;
		String result = "";
		String sep = "";
		Iterator<String> matures = _matureMiRNAs.keySet().iterator();
		//System.out.println("\n>>"+readID+"\t"+_reads.get(readID).getReferenceID()+"\t"+_reads.get(readID).getStart()+"\t"+_reads.get(readID).getStop()+"\t"+_reads.get(readID).isNegativeStrand());
		while(matures.hasNext()){
			String tmpID = matures.next();
			_matureMiRNAs.get(tmpID).setReferenceID(_id);
			//System.out.println(">"+tmpID+"\t"+_matureMiRNAs.get(tmpID).getReferenceID()+"\t"+_matureMiRNAs.get(tmpID).getStart()+"\t"+_matureMiRNAs.get(tmpID).getStop()+"\t"+_matureMiRNAs.get(tmpID).isNegativeStrand());

			if(_reads.get(readID).overlaps(_matureMiRNAs.get(tmpID), minFractionForOverlap)){
				result = result.concat(sep+tmpID);
				sep = "|";
			}
		}
		return result;
	}

}

class Alignment{
	private String _referenceID, _querySequence;
	private int _start, _stop;
	private boolean _isNegativeStrand;
	private Cigar _cigar;
	public Alignment(String referenceID, int start, int stop, boolean isNegativeStrand, Cigar cigar, String querySequence){
		_referenceID = referenceID;
		_start = start;
		_stop = stop;
		_isNegativeStrand = isNegativeStrand;
		_cigar = cigar;
		_querySequence = querySequence;
	}
	public Alignment(String referenceID, int start, int stop, boolean isNegativeStrand){
		_referenceID = referenceID;
		_start = start;
		_stop = stop;
		_isNegativeStrand = isNegativeStrand;
	}

	public void setReferenceID(String newID){ _referenceID = newID; }
	public String getReferenceID(){ return _referenceID; }
	public int getStart(){ return _start; }
	public int getStop(){ return _stop; }
	public boolean isNegativeStrand(){ return _isNegativeStrand; }
	public Cigar getCIGAR(){ return _cigar; }
	public String getSequence(){ return _querySequence; }

	public boolean overlaps(Alignment otherAlignment, double minFracOverlap){
		boolean result = false;
		if(_referenceID.equals(otherAlignment.getReferenceID())){
			if(_isNegativeStrand == otherAlignment.isNegativeStrand()){
				if((_start >= otherAlignment.getStart()  &&  _start <= otherAlignment.getStop())  ||  (_stop >= otherAlignment.getStart()  &&  _stop <= otherAlignment.getStop())){
					int diff = Math.abs(_start - otherAlignment.getStart()) + Math.abs(_stop - otherAlignment.getStop()); 
					int max = Math.max(_stop, otherAlignment.getStop()) - Math.min(_start, otherAlignment.getStart());
					//System.out.print(1.0 - (diff / (max+0.0))+"\t"+minFracOverlap+"\t"+(1.0 - (diff / (max+0.0)) >= minFracOverlap));
					if((1.0 - (diff / (max+0.0)) >= minFracOverlap))
						result = true;
				}
			}
		}
		return result;
	}
}




