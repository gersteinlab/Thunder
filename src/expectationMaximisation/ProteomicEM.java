package expectationMaximisation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Iterator;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import main.Thunder;
import objects.Alignment;
import objects.MS1_scan;
import objects.Peptide;
import objects.Spectra;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import proteome.ProcessPepDigest;
import proteome.ReadXTandem;
import proteome.SAXHandler_mzXML;
import proteome.SpectraAlignmentEngine;

public class ProteomicEM {
	public boolean verbose = false;

	private String _outputPrefix;
	public ProteomicEM(String outputPrefix){
		_outputPrefix = outputPrefix;
	}


	/*private HashMap<String, HashMap<Integer,Integer>> _isoformPeptideCounts;
	public void readDigestInfo(File input) throws IOException{
		System.err.println(this.getTime()+" Reading peptide digest info: "+input.getAbsolutePath());
		_isoformPeptideCounts = ProcessPepDigest.processEMBOSSoutput(input);
		System.err.println(this.getTime()+" Done.");
	}*/
	private HashMap<String, Integer> _isoformPeptideCounts = new HashMap<String, Integer>();
	public void readDigestInfo(File input) throws IOException{
		Thunder.printLineErr("Performing in-silico peptide digest: "+input.getAbsolutePath());
		ArrayList<String> enzymes = new ArrayList<String>();
		enzymes.add(ProcessPepDigest.ENZYME_TRYPSIN);
		enzymes.add(ProcessPepDigest.ENZYME_LYSC);
		//_isoformPeptideCounts = ProcessPepDigest.digestProteinsFromFasta(input, enzymes);

		HashMap<String, Integer> tmpPeptides = ProcessPepDigest.digestProteinsFromFasta(input, enzymes);

		Iterator<String> it = tmpPeptides.keySet().iterator();
		String tmpID;
		String[] isoformFrameGene;
		while(it.hasNext()){
			tmpID = it.next();
			isoformFrameGene = SpectraAlignmentEngine.parseIsoformID(tmpID);
			_isoformPeptideCounts.put(isoformFrameGene[0]+"_"+isoformFrameGene[1], tmpPeptides.get(tmpID));
		}

		Thunder.printLineErr("Done.");
	}


	private HashMap<String, MS1_scan> _ms1Intensities = new HashMap<String, MS1_scan>();
	/**
	 * 
	 * @param input
	 */
	public void readMZXML(File input){
		Thunder.printLineErr("Reading spectra intensities (MS1): "+input.getAbsolutePath());

		try {
			SAXParserFactory factory = SAXParserFactory.newInstance();
			SAXParser saxParser = factory.newSAXParser();
			SAXHandler_mzXML handler = new SAXHandler_mzXML();

			saxParser.parse(input, handler);
			_ms1Intensities = handler.getMS1Intensities();
			//System.err.println();

			Thunder.printLineErr("N spectra: "+_ms1Intensities.size());

		}catch (Exception e) {
			e.printStackTrace();
		}finally{
		}

		/*Iterator<String> it = ms1Intensities.keySet().iterator();
		for(int i=0;i<10;i++){
			System.out.println(it.next());
		}*/

		Thunder.printLineErr("Done.");
	}



	/**
	 * 
	 * @param input
	 * @throws Exception
	 */
	public void readSpectra(File tandemResults, double fdrMax, boolean restrictToCDS) throws Exception{
		Thunder.printLineErr("Reading spectra alignments: "+tandemResults.getAbsolutePath());

		// Read spectra
		SpectraAlignmentEngine spectraAlignments = ReadXTandem.readTandemXML(tandemResults, fdrMax, true);

		String thisGeneID = "";
		String thisTranscriptID = "";
		String thisSpectraID = "";
		String thisFrame = "";
		ArrayList<String> tmp_transcripts = new ArrayList<String>();

		Alignment currentAlignment = null;
		Iterator<String> allSpectra = spectraAlignments.getAlignments().keySet().iterator();
		Iterator<Alignment> currentSpectraAlignments;

		ArrayList<String> tmp_genes;

		String[] isoformFrameGene;
		ArrayList<String> unknownIsoformIDs = new ArrayList<String>();
		while(allSpectra.hasNext()){
			currentSpectraAlignments = spectraAlignments.getAlignments().get(allSpectra.next()).iterator();

			int spectraIntensity = 1;
			boolean hasMS1 = false;

			tmp_genes = new ArrayList<String>();
			tmp_transcripts = new ArrayList<String>();

			while(currentSpectraAlignments.hasNext()){
				currentAlignment = currentSpectraAlignments.next();

				isoformFrameGene = SpectraAlignmentEngine.parseIsoformID(currentAlignment.getProtein().getIsoformID());
				thisTranscriptID = isoformFrameGene[0];
				thisGeneID = isoformFrameGene[2];
				thisFrame = isoformFrameGene[1];
				//thisSpectraID = currentAlignment.getSpectra().getAttribute(Spectra.ATTRIBUTE_ID);
				thisSpectraID = currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_SEQ);


				if(transcriptID_2_geneID.containsKey(thisTranscriptID+"_"+thisFrame)){					
					thisGeneID = transcriptID_2_geneID.get(thisTranscriptID+"_"+thisFrame);

					if(!tmp_genes.contains(thisGeneID))
						tmp_genes.add(thisGeneID);

					/*if(thisGeneID.equals("ENSG00000181652.16")){
						System.out.println("Gene has CDS: "+this.genes.get(thisGeneID).hasCDS());
						System.out.println(thisTranscriptID+" biotype: "+this.genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getBiotype());
						System.out.println(thisTranscriptID+" has CDS: "+this.genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).hasCDS());
					}*/

					if(restrictToCDS  &&  this.genes.get(thisGeneID).hasCDS()){
						// If this gene has a CDS, and we only want CDS peptides, only allow alignments to coding transcripts
						int peptideStart = Integer.valueOf(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_START)).intValue();
						int peptideEnd = Integer.valueOf(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_END)).intValue();
						if(this.genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).hasCDS()  &&
								this.genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getBiotype().equals("protein_coding")  &&
								peptideStart >= this.genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getCDSStart()  &&
								peptideEnd <= this.genes.get(thisGeneID).getTranscript(thisTranscriptID+"_"+thisFrame).getCDSStop()){
							tmp_transcripts.add(thisTranscriptID+"_"+thisFrame);
						}
					}else{
						// if this gene has no annotated CDS, add all peptides
						tmp_transcripts.add(thisTranscriptID+"_"+thisFrame);
					}


				}else{
					// this transcript is not in the GTF! 
					if(!unknownIsoformIDs.contains(currentAlignment.getProtein().getIsoformID())){	
						unknownIsoformIDs.add(currentAlignment.getProtein().getIsoformID());
						System.err.println("WARNING: Unable to find annotation information for "+currentAlignment.getProtein().getIsoformID());
					}
				}

				if(!hasMS1){
					if(_ms1Intensities.containsKey(currentAlignment.getSpectra().getAttribute(Spectra.ATTRIBUTE_ID))){
						hasMS1 = true;
						spectraIntensity = (int)Math.round(Double.valueOf(_ms1Intensities.get(currentAlignment.getSpectra().getAttribute(Spectra.ATTRIBUTE_ID)).getPrecursorAttribute(MS1_scan.ATTRIBUTE_PRECURSOR_INTENSITY)).doubleValue());
						//System.out.println(currentAlignment.getSpectra().getAttribute(Spectra.ATTRIBUTE_ID)+": "+spectraIntensity);
					}
				}
			}

			// add peptides to gene if there is only one gene:
			if(tmp_genes.size() == 1){
				if(hasMS1)
					this.genes.get(thisGeneID).addRead(thisSpectraID, tmp_transcripts, spectraIntensity);
				else
					this.genes.get(thisGeneID).addRead(thisSpectraID, tmp_transcripts);
			}

		}
		Thunder.printLineErr("Done.");
	}





	/**
	 * 
	 * @param input
	 * @throws Exception
	 */
	public void readGTF(File input) throws Exception{
		Thunder.printLineErr("Reading GTF: "+input.getAbsolutePath());
		BufferedReader in = new BufferedReader(new FileReader(input));
		String line = "";

		while((line=in.readLine()) != null){
			if(!line.startsWith("##")){
				parseLine(line);
			}
		}
		in.close();
		Thunder.printLineErr("Done.");
	}

	private HashMap<String, EM_Gene> genes = new HashMap<String, EM_Gene>();
	private HashMap<String, String> transcriptID_2_geneID = new HashMap<String, String>();



	/*
	 * Read GTF entries to memory - IMPORTANT: read all transcripts in THREE FRAMES!
	 */
	private void parseLine(String line){
		String[] bits = line.split(" |\t");;

		String featureType = bits[2].trim(); 
		if(featureType.equals("exon")  ||  featureType.equals("CDS")){

			String geneID = trimAttribute(bits[9].trim());
			String transcriptID = trimAttribute(bits[11].trim());

//			if(geneID.equals("ENSG00000181652.16")){
//				System.out.println(featureType);
//			}

			// Create gene
			if(!genes.containsKey(geneID))
				genes.put(geneID, new EM_Gene(geneID));

			// Placeholder variables
			int tmp_peptideCount_1 = 0;
			int tmp_peptideCount_2 = 0;
			int tmp_peptideCount_3 = 0;

			// Add transcript in all three frames to this gene
			if(!genes.get(geneID).containsTranscript(transcriptID+"_1")){
				genes.get(geneID).addTranscript(transcriptID+"_1");
				genes.get(geneID).addTranscript(transcriptID+"_2");
				genes.get(geneID).addTranscript(transcriptID+"_3");

				// Set transcript biotype
				String tmpAtt = trimAttribute(bits[19].trim());
				genes.get(geneID).getTranscript(transcriptID+"_1").setBiotype(tmpAtt);
				genes.get(geneID).getTranscript(transcriptID+"_2").setBiotype(tmpAtt);
				genes.get(geneID).getTranscript(transcriptID+"_3").setBiotype(tmpAtt);

				//System.out.println(transcriptID+"_1");
				if(_isoformPeptideCounts.containsKey(transcriptID+"_1"))
					tmp_peptideCount_1 = _isoformPeptideCounts.get(transcriptID+"_1");
				if(_isoformPeptideCounts.containsKey(transcriptID+"_2"))
					tmp_peptideCount_2 = _isoformPeptideCounts.get(transcriptID+"_2");
				if(_isoformPeptideCounts.containsKey(transcriptID+"_3"))
					tmp_peptideCount_3 = _isoformPeptideCounts.get(transcriptID+"_3");

				// Add gene-transcript map
				if(!transcriptID_2_geneID.containsKey(transcriptID+"_1")){
					transcriptID_2_geneID.put(transcriptID+"_1", geneID);
					transcriptID_2_geneID.put(transcriptID+"_2", geneID);
					transcriptID_2_geneID.put(transcriptID+"_3", geneID);
				}
			}

			// Add transcript or CDS length in nucleotides
			if(featureType.equals("CDS")){
				genes.get(geneID).getTranscript(transcriptID+"_1").addCDS(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				genes.get(geneID).getTranscript(transcriptID+"_2").addCDS(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				genes.get(geneID).getTranscript(transcriptID+"_3").addCDS(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				genes.get(geneID).setProteinCoding(true);
			}else if(featureType.equals("exon")){
				genes.get(geneID).getTranscript(transcriptID+"_1").addExon(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				genes.get(geneID).getTranscript(transcriptID+"_2").addExon(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				genes.get(geneID).getTranscript(transcriptID+"_3").addExon(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
			}

			
			
			// Reset exon lengths based on in-silico digest
			genes.get(geneID).getTranscript(transcriptID+"_1").setExonLength(tmp_peptideCount_1);
			genes.get(geneID).getTranscript(transcriptID+"_2").setExonLength(tmp_peptideCount_2);
			genes.get(geneID).getTranscript(transcriptID+"_3").setExonLength(tmp_peptideCount_3);
			/*
			 * TODO: calculate peptides from the CDS sequence, rather than from the whole transcript
			 */
			genes.get(geneID).getTranscript(transcriptID+"_1").setCodingExonLength(tmp_peptideCount_1);
			genes.get(geneID).getTranscript(transcriptID+"_2").setCodingExonLength(tmp_peptideCount_2);
			genes.get(geneID).getTranscript(transcriptID+"_3").setCodingExonLength(tmp_peptideCount_3);
			






			// Add feature length to the transcript/CDS length
			//int length = Integer.valueOf(bits[4]).intValue()-Integer.valueOf(bits[3]).intValue();
			//if(bits[2].equals("CDS"))
			//	genes.get(geneID).getTranscript(transcriptID).addCodingExonLength(length);
			//else
			//	genes.get(geneID).getTranscript(transcriptID).addExonLength(length);
		}
	}
	public static String trimAttribute(String in){
		return in.replaceAll("^\"|\";$", "");
	}




	/**
	 * 
	 * @param inputFile
	 * @throws Exception
	 */
	public void readTranscriptExpressions(File inputFile) throws IOException{
		Thunder.printLineErr("Reading expression data for priors: "+inputFile.getAbsolutePath());
		BufferedReader in = new BufferedReader(new FileReader(inputFile));
		String line = "";
		String[] bits;


		// determine from the file header if this is an RNA-seq (eXpress) output or a Thunder output form the footprint EM
		String inputDiscriminator = in.readLine().split("\t")[0];

		while((line=in.readLine()) != null){
			bits = line.split("\t");

			//transcriptQuants.put(bits[1].trim(), Double.valueOf(bits[14].trim()).doubleValue());
			String transcriptID = bits[1].trim();
			if(this.transcriptID_2_geneID.containsKey(transcriptID+"_1")){
				double tmpTPM = 0.0;
				if(inputDiscriminator.equals("bundle_id")){		// this is an eXpress file
					tmpTPM = Double.valueOf(bits[14].trim()).doubleValue();
				}else if(inputDiscriminator.equals("geneID")){		// this is a Thunder footprintEM file
					// "geneID\ttranscriptID\ttranscriptLength_exons\ttranscriptLength_CDS\tnFootprintsMappedToThisGene\tnReadsMappedToThisTranscript\thasFlatPrior\tprior\tflag\tnIterations\teffectiveReadCount\treadDensity_transcript\treadDensity_CDS\ttranscriptFractionAfterFootprints\ttranscriptBiotype");
					// use the transcript fraction as the TPM, this will get translated to the same prior
					tmpTPM = Double.valueOf(bits[13].trim()).doubleValue();
				}else{
					// Dunno!
				}
				this.genes.get(this.transcriptID_2_geneID.get(transcriptID+"_1")).getTranscript(transcriptID+"_1").setTPM(tmpTPM);
				this.genes.get(this.transcriptID_2_geneID.get(transcriptID+"_2")).getTranscript(transcriptID+"_2").setTPM(tmpTPM);
				this.genes.get(this.transcriptID_2_geneID.get(transcriptID+"_3")).getTranscript(transcriptID+"_3").setTPM(tmpTPM);

			}
		}
		in.close();

		Thunder.printLineErr("Done.");
	}



	/**
	 * 
	 */
	public void initialisePriors(){
		//System.out.println(this.getTime()+" Initialising priors...");

		Iterator<String> iterator = genes.keySet().iterator();
		String geneID = "";

		while(iterator.hasNext()){
			geneID = iterator.next();
			genes.get(geneID).setPriors();
		}

		//System.out.println(this.getTime()+" Done.");
	}


	/**
	 * 
	 * @param maxIterations
	 * @param converganceDistance
	 */
	public void doEM(int maxIterations, double converganceDistance, boolean outputAll) throws IOException{
		Thunder.printLineErr("Running EM (maxIterations: "+maxIterations+"), writing to: "+this._outputPrefix+".exprs");
		Iterator<String> iterator = genes.keySet().iterator();
		int totalGenes = genes.size();
		int count = 0;
		int percent = 0;
		String geneID = "";
		EM_Core_optimised emCore;

		PrintWriter out = new PrintWriter(new FileWriter(this._outputPrefix+".exprs"));

		if(_ms1Intensities.size() == 0)
			out.println("geneID\ttranscriptID\tframe\tnPeptidesMappedToThisGene\tnumberObservedPeptides\tnumberObservablePeptides\thasFlatPrior\tprior\tflag\tnIterations\teffectivePeptideCount\ttranscriptFractionAfterPeptides\ttranscriptBiotype");
		else
			out.println("geneID\ttranscriptID\tframe\tnPeptidesMappedToThisGene\tnumberObservedPeptides\tnumberObservablePeptides\thasFlatPrior\tprior\tflag\tnIterations\taveragePrecursorIntensity\ttranscriptFractionAfterPeptides\ttranscriptBiotype");

		while(iterator.hasNext()){
			geneID = iterator.next();

			// Change the effective exon length based on the number of mapped peptides (use # theoretical peptides when the in silico digest works well enough...)
			if(genes.get(geneID).getReadCounts().size() > 0){
				for(int i=0;i<genes.get(geneID).getTranscriptIDs().size();i++){
					genes.get(geneID).getTranscript(genes.get(geneID).getTranscriptIDs().get(i)).setExonLength(genes.get(geneID).getReadCounts().size());
					genes.get(geneID).getTranscript(genes.get(geneID).getTranscriptIDs().get(i)).setCodingExonLength(genes.get(geneID).getReadCounts().size());
				}
			}

			// Run the EM
			emCore = new EM_Core_optimised(genes.get(geneID));
			emCore.setReadLength(0);
			EM_Result res = emCore.runEM(maxIterations, converganceDistance);//, out);

			// Format and print EM results
			if(outputAll  ||  res.getGene().getReadCounts().size() > 0){

				for(int i=0;i<res.getGene().getTranscriptIDs().size();i++){
					String thisTranscriptID = res.getGene().getTranscriptIDs().get(i);
					int nPeptides = res.getGene().getTranscriptReadCounts(thisTranscriptID);
					double tmp_avgPrecursorIntensity = 0.0;
					if(nPeptides > 0){
						if(_ms1Intensities.size() == 0)
							tmp_avgPrecursorIntensity = res.getFinalEffectiveReadCount().get(thisTranscriptID).doubleValue()+0.0;
						else
							tmp_avgPrecursorIntensity = (res.getFinalEffectiveReadCount().get(thisTranscriptID).doubleValue()+0.0) / nPeptides;
					}
					int tmp_coverage = res.getGene().getLengths().get(thisTranscriptID);
					double tmp_finalLikelihood = res.getFinalTranscriptLikelihoods().get(thisTranscriptID).doubleValue();
					//out.printf(geneID+"\t"+thisTranscriptID+"\t"+res.getGene().hasFlatPrior()+"\t%e\t"+res.getFlag()+"\t"+res.getNIterations()+"\t%f"+nPeptides+"\t%e\t%e\n", res.getGene().getPriors().get(thisTranscriptID), tmp_avgPrecursorIntensity, tmp_coverage, tmp_finalLikelihood);
					String[] transcriptIDbits = thisTranscriptID.split("_");
					out.printf(geneID+"\t"+transcriptIDbits[0]+"\t"+transcriptIDbits[1]+"\t"+res.getGene().getReadCounts().size()+"\t"+nPeptides+"\t"+tmp_coverage+"\t"+res.getGene().hasFlatPrior()+"\t%e\t"+res.getFlag()+"\t"+res.getNIterations()+"\t%f\t%e\t"+res.getGene().getTranscript(thisTranscriptID).getBiotype()+"\n", res.getGene().getPriors().get(thisTranscriptID), tmp_avgPrecursorIntensity, tmp_finalLikelihood);
				}
			}

			// Update the status bar
			count++;
			if(Math.round((count*100.0)/totalGenes) > percent){
				percent = (int)Math.round((count*100.0)/totalGenes);
				Thunder.printProgressBar(percent);
			}
		}
		System.err.println("");
		out.flush();
		out.close();
		Thunder.printLineErr("Done.");
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withLongOpt("annotation").withArgName("path").hasArg().withDescription("Path to the GTF file containing the transcript/gene relationship").create(Thunder.OPT_PATH_ANNOTATION));
		options.addOption(OptionBuilder.withLongOpt("fasta").withArgName("path").hasArg().withDescription("Path to the fasta file containing [AA] isoform sequences so that we can predict expected peptide products").create("s"));
		options.addOption(OptionBuilder.withLongOpt("priors").withArgName("path").hasArg().withDescription("[optional] Text file containing RNA-seq transcript expression quantifications from eXpress").create("e"));
		options.addOption(OptionBuilder.withLongOpt("spectra_mapped").withArgName("path").hasArg().withDescription("X!Tandem output XML file containing spectra alignments").create("f"));		
		options.addOption(OptionBuilder.withLongOpt("spectra_raw").withArgName("path").hasArg().withDescription("[optional] mzXML file containing containing raw spectra peaks (required for quantification)").create("mzxml"));		
		options.addOption(OptionBuilder.withLongOpt("outputPrefix").withArgName("outputPath").hasArg().withDescription("Path to which to output results").create("o"));
		options.addOption(OptionBuilder.withDescription("Count only peptides in the CDS (where available)").create("cds"));
		options.addOption(OptionBuilder.withArgName("int").hasArg().withDescription("Maximum number of EM iterations [default: 500]").create("N"));
		options.addOption(OptionBuilder.withArgName("double").hasArg().withDescription("Criteria for EM convergence [default: 0.0001]").create("c"));
		options.addOption(OptionBuilder.withArgName("double").hasArg().withDescription("[default: 0.05] Based on hits to the reverse decoy DB, define an acceptable false discovery rate (FDR)").create("fdr"));
		options.addOption(OptionBuilder.withDescription("Write all genes/isoforms, even those with no observed peptide spectra").create("v"));
		return options;
	}




	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {

		/*
		args = new String[]{"IsoformEM_Proteomics", 
				//"-f", "/Users/robk/Desktop/EM_TEST/ProteoData/allSpectra_C8.txt",
				"-f", "/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Proteomics/HEK.C8.Orbi.gencode21.TandemOutput.xml",
				"-mzxml", "/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Proteomics/EOT12-0395.HEK.C8.Orbi.mzxml",
				"-a", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation_noSelenocysteine.gtf",
				//"-e", "/Users/robk/Box Sync/Work/HEK293_RNAseq/STAR_eXpress/A1-Total_2x75/results.xprs",
				"-e", "/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/A1_totalRNA/results.xprs",
				"-s", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation.protein.fa",
				"--outputPrefix","/Users/robk/WORK/YALE_offline/My_Proteomics/HEK/RNA-seq_NEW/Proteomics/ProteomicEM/TEST_TEST_TEST"};
		 */
		//"-o", "/Users/robk/Desktop/EM_TEST/TEST_NEW",
		//"-DEV"};


		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption(Thunder.OPT_PATH_ANNOTATION)  &&  cmdArgs.hasOption("f")  &&  cmdArgs.hasOption("s") && cmdArgs.hasOption("o")){
			System.err.println();

			ProteomicEM engine = new ProteomicEM(cmdArgs.getOptionValue("o"));

			// Clean up
			//if(cmdArgs.hasOption("v")){
			//engine.removeFile(new File(tmp_path));
			//	engine.verbose = true;
			//}

			//String output_path = cmdArgs.getOptionValue("o");
			//File output_priors = new File(output_path+"/priors.txt");
			//File output_alignments = new File(tmp_path+"/alignments.txt");


			// Create the output directory and temp dir if they don't already exist
			//engine.makeOutputDirectories(output_path);

			// read the number of expected peptides per isoform per frame
			engine.readDigestInfo(new File(cmdArgs.getOptionValue("s")));

			// Read the transcripts in the GTF
			engine.readGTF(new File(cmdArgs.getOptionValue(Thunder.OPT_PATH_ANNOTATION)));

			// Read the eXpress transcript quantifications or choose flat priors
			if(cmdArgs.hasOption("e")){
				engine.readTranscriptExpressions(new File(cmdArgs.getOptionValue("e")));
			}

			// Convert the transcript expressions to priors and output
			engine.initialisePriors();

			// Set desired FDR for spectra
			double fdrMax = 0.05;
			if(cmdArgs.hasOption("fdr"))
				fdrMax = Double.valueOf(cmdArgs.getOptionValue("fdr")).doubleValue();

			// If provided an mzXML file containing precursor intensities, read it now
			if(cmdArgs.hasOption("mzxml"))
				engine.readMZXML(new File(cmdArgs.getOptionValue("mzxml")));

			// Parse the MS/MS alignments
			engine.readSpectra(new File(cmdArgs.getOptionValue("f")), fdrMax, cmdArgs.hasOption("cds"));

			// set the number of iterations to cap the EM
			int maxIterations = 500;
			if(cmdArgs.hasOption("N"))
				maxIterations = Integer.valueOf(cmdArgs.getOptionValue("N")).intValue();

			// set the convergence criteria for the EM
			double convergenceDistance = 1.0/10000.0;
			if(cmdArgs.hasOption("c"))
				convergenceDistance = Double.valueOf(cmdArgs.getOptionValue("c")).doubleValue();

			boolean outputAll = false;
			if(cmdArgs.hasOption("v"))
				outputAll = true;

			// Run the EM!
			engine.doEM(maxIterations, convergenceDistance, outputAll);

			Thunder.printLineErr("All Done!");

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.setWidth(200);
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" IsoformEM_Proteomics", getCmdLineOptions());
			System.err.println();
		}


	}

	/*public String getTime(){
		return((new SimpleDateFormat("yyyy/MM/dd HH:mm:ss")).format(Calendar.getInstance().getTime()));
	}*/




	public void makeOutputDirectories(String output_path){
		File tmp = new File(output_path);
		if(!tmp.exists())
			tmp.mkdir();
		//tmp = new File(output_path+"/tmp");
		//if(!tmp.exists())
		//	tmp.mkdir();
	}

	public void removeFile(File file){
		if(file.isDirectory()){
			File[] children = file.listFiles();
			for(int i=0;i<children.length;i++)
				removeFile(children[i]);
		}
		file.delete();
		//System.out.println("Deleted: "+file.getAbsolutePath());
	}

}



