package expectationMaximisation;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
import java.text.SimpleDateFormat;
import java.util.ArrayList;
import java.util.Calendar;
import java.util.HashMap;
import java.util.Iterator;

import main.Thunder;
import net.sf.samtools.SAMFileHeader.SortOrder;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMFileWriter;
import net.sf.samtools.SAMFileWriterFactory;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

public class FootprintEM {
	public boolean verbose = false;
	private String _outputPrefix;
	public FootprintEM(String outputPrefix){
		_outputPrefix = outputPrefix;
	}

	public String recordSummary(SAMRecord record){
		return record.getReadName()+"\t"+record.getReferenceName()+"\t"+record.getIntegerAttribute("NH")+"\t"+record.getIntegerAttribute("HI");
	}


	public void readSAM(File input, boolean restrictToCDS) throws Exception{
		System.err.println(this.getTime()+" Reading alignments: "+input.getAbsolutePath());
		SAMFileReader inputSam = new SAMFileReader(input);
		SAMRecord thisRecord;
		SAMRecordIterator it = inputSam.iterator();

		if(verbose)
			System.err.println("BAM sorted by: "+inputSam.getFileHeader().getSortOrder().toString());

		if(inputSam.getFileHeader().getSortOrder().compareTo(SortOrder.queryname) != 0){
			System.err.println("BAM must be sorted by readID ('queryname')\nThis BAM is sorted by: "+inputSam.getFileHeader().getSortOrder().toString());
			System.exit(0);
		}

		String thisGeneID = "";
		String thisTranscriptID = "";
		String lastReadID = "";
		String thisReadID = "";
		boolean multimapper = false;
		ArrayList<String> tmp_transcripts = new ArrayList<String>();

		while(it.hasNext()){
			thisRecord = it.next();
			thisTranscriptID = thisRecord.getReferenceName();
			thisReadID = thisRecord.getReadName();

			if(transcriptID_2_geneID.containsKey(thisTranscriptID)){

				if(!thisReadID.equals(lastReadID)){  		// new read
					if(lastReadID.length() > 0){     		// not first read
						if(!multimapper){            		// maps to only one gene
							if(tmp_transcripts.size() > 0){ // read aligns to at least one transcript (may not happen for CDS only filtering)
								// add read to gene:
								this.genes.get(thisGeneID).addRead(lastReadID, tmp_transcripts);
							}
						}
					}
					
					lastReadID = thisReadID;
					thisGeneID = transcriptID_2_geneID.get(thisTranscriptID);
					multimapper = false;
					tmp_transcripts = new ArrayList<String>();
					
					
					//tmp_transcripts.add(thisTranscriptID);
					if(restrictToCDS  &&  this.genes.get(thisGeneID).hasCDS()){
						// If this gene has a CDS, and we only want CDS reads, only allow alignments to coding transcripts
						if(this.genes.get(thisGeneID).getTranscript(thisTranscriptID).hasCDS()  &&
								this.genes.get(thisGeneID).getTranscript(thisTranscriptID).getBiotype().equals("protein_coding")  &&
								thisRecord.getAlignmentStart() >= this.genes.get(thisGeneID).getTranscript(thisTranscriptID).getCDSStart()  &&
								thisRecord.getAlignmentEnd() <= this.genes.get(thisGeneID).getTranscript(thisTranscriptID).getCDSStop()){
							tmp_transcripts.add(thisTranscriptID);
						}
					}else{
						// if this gene has no annotated CDS, add all reads 
						tmp_transcripts.add(thisTranscriptID);
					}
			
					
				}else{
					if( ! transcriptID_2_geneID.get(thisTranscriptID).equals(thisGeneID))
						multimapper = true;

					if(!multimapper){
						if(restrictToCDS  &&  this.genes.get(thisGeneID).hasCDS()){
							// If this gene has a CDS, and we only want CDS reads, only allow alignments to coding transcripts
							if(this.genes.get(thisGeneID).getTranscript(thisTranscriptID).hasCDS()  &&
									this.genes.get(thisGeneID).getTranscript(thisTranscriptID).getBiotype().equals("protein_coding")  &&
									thisRecord.getAlignmentStart() >= this.genes.get(thisGeneID).getTranscript(thisTranscriptID).getCDSStart()  &&
									thisRecord.getAlignmentEnd() <= this.genes.get(thisGeneID).getTranscript(thisTranscriptID).getCDSStop()){
								tmp_transcripts.add(thisTranscriptID);
							}
						}else{
							// if this gene has no annotated CDS, add all reads 
							tmp_transcripts.add(thisTranscriptID);
						}
					}
				}

				/*if(thisGeneID.equals("ENSG00000023572.4")){
				System.out.println("thisReadID = "+thisReadID);
				System.out.println("thisTranscriptID = "+thisTranscriptID);
				System.out.println("multimapper = "+multimapper);
			}*/
			}

		}

		// for the last read!
		if(lastReadID.length() > 0){     		// not first read
			if(!multimapper){            		// maps to only one gene
				if(tmp_transcripts.size() > 0){ // read aligns to at least one transcript (may not happen for CDS only filtering)
					//System.out.println("LAST: thisReadID = "+thisReadID);
					//System.out.println("LAST: thisTranscriptID = "+thisTranscriptID);
					//System.out.println("LAST: multimapper = "+multimapper);
					
					// add read to gene:
					this.genes.get(thisGeneID).addRead(thisReadID, tmp_transcripts);
				}
			}
		}

		inputSam.close();
		System.err.println(this.getTime()+" Done.");
	}






	public void readGTF(File input) throws Exception{
		System.err.println(this.getTime()+" Reading GTF: "+input.getAbsolutePath());
		BufferedReader in = new BufferedReader(new FileReader(input));
		String line = "";

		while((line=in.readLine()) != null){
			if(!line.startsWith("##")){
				parseLine(line);
			}
		}
		in.close();
		System.err.println(this.getTime()+" Done.");
	}

	private HashMap<String, EM_Gene> genes = new HashMap<String, EM_Gene>();
	private HashMap<String, String> transcriptID_2_geneID = new HashMap<String, String>();


	private void parseLine(String line){
		String[] bits = line.split(" |\t");;

		if(bits[2].equals("exon")  ||  bits[2].equals("CDS")){
			//if(bits[2].equals("CDS")){

			String geneID = trimAttribute(bits[9].trim());
			String transcriptID = trimAttribute(bits[11].trim());

			// Create gene
			if(!genes.containsKey(geneID))
				genes.put(geneID, new EM_Gene(geneID));

			// Add transcript to gene
			if(!genes.get(geneID).containsTranscript(transcriptID))
				genes.get(geneID).addTranscript(transcriptID);

			// Add gene-transcript map
			if(!transcriptID_2_geneID.containsKey(transcriptID))
				transcriptID_2_geneID.put(transcriptID, geneID);

			// Set transcript biotype
			genes.get(geneID).getTranscript(transcriptID).setBiotype(trimAttribute(bits[19].trim()));

			// Add feature length to the transcript/CDS length
			//int length = Integer.valueOf(bits[4]).intValue()-Integer.valueOf(bits[3]).intValue();
			if(bits[2].equals("CDS")){
				//genes.get(geneID).getTranscript(transcriptID).addCodingExonLength(length);
				genes.get(geneID).getTranscript(transcriptID).addCDS(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
				genes.get(geneID).setProteinCoding(true);
			}else{
				//genes.get(geneID).getTranscript(transcriptID).addExonLength(length);
				genes.get(geneID).getTranscript(transcriptID).addExon(Integer.valueOf(bits[3]).intValue(), Integer.valueOf(bits[4]).intValue());
			}
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
		System.err.println(this.getTime()+" Reading eXpress expression data: "+inputFile.getAbsolutePath());
		BufferedReader in = new BufferedReader(new FileReader(inputFile));
		String line = "";
		String[] bits;

		// burn the first line
		in.readLine();

		while((line=in.readLine()) != null){
			bits = line.split("\t");

			//transcriptQuants.put(bits[1].trim(), Double.valueOf(bits[14].trim()).doubleValue());
			String transcriptID = bits[1].trim();
			if(this.transcriptID_2_geneID.containsKey(transcriptID))
				this.genes.get(this.transcriptID_2_geneID.get(transcriptID)).getTranscript(transcriptID).setTPM(Double.valueOf(bits[14].trim()).doubleValue());;
		}
		in.close();

		System.err.println(this.getTime()+" Done.");
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

	private HashMap<String, String> _readMap = new HashMap<String, String>();
	public void doEM(int maxIterations, double converganceDistance, boolean outputAll) throws IOException{
		//System.err.println(getTime()+" Running EM (maxIterations: "+maxIterations+"), writing to: "+this._outputPrefix+".exprs");
		Thunder.printLineErr("Running EM (maxIterations: "+maxIterations+"), writing to: "+this._outputPrefix+".exprs");
		Iterator<String> iterator = genes.keySet().iterator();
		int totalGenes = genes.size();
		int count = 0;
		int percent = 0;
		String geneID = "";
		EM_Core_optimised emCore;

		PrintWriter out = new PrintWriter(new FileWriter(this._outputPrefix+".exprs"));
		//System.out.println("geneID\ttranscriptID\thasFlatPrior\tprior\tflag\tnIterations\teffectiveReadCount\treadDensity\ttranscriptFractionAfterFootprints");
		out.println("geneID\ttranscriptID\ttranscriptLength_exons\ttranscriptLength_CDS\tnFootprintsMappedToThisGene\tnReadsMappedToThisTranscript\thasFlatPrior\tprior\tflag\tnIterations\teffectiveReadCount\treadDensity_transcript\treadDensity_CDS\ttranscriptFractionAfterFootprints\ttranscriptBiotype");

		while(iterator.hasNext()){
			geneID = iterator.next();
			emCore = new EM_Core_optimised(genes.get(geneID));
			EM_Result res = emCore.runEM(maxIterations, converganceDistance);//, out);
			_readMap = emCore.sampleReadAssignments(_readMap);

			// Format and print EM results
			if(outputAll  ||  res.getGene().getReadCounts().size() > 0){
				for(int i=0;i<res.getGene().getTranscriptIDs().size();i++){
					String thisTranscriptID = res.getGene().getTranscriptIDs().get(i);
					double tmp_effectiveCount = res.getFinalEffectiveReadCount().get(thisTranscriptID).doubleValue()+0.0;
					//double tmp_coverage = tmp_effectiveCount / res.getGene().getLengths().get(thisTranscriptID);
					double tmp_finalLikelihood = res.getFinalTranscriptLikelihoods().get(thisTranscriptID).doubleValue();
					int exonlength = res.getGene().getTranscript(thisTranscriptID).getExonLength();
					int cdslength = res.getGene().getTranscript(thisTranscriptID).getCodingExonLength();
					double density_exon = tmp_effectiveCount/(exonlength+0.0);
					double density_cds = tmp_effectiveCount/(cdslength+0.0);
					if(tmp_effectiveCount == 0.0)
						density_exon = density_cds = 0.0;
					else if(cdslength == 0)
						density_cds = 0.0;
					out.printf(geneID+"\t"+thisTranscriptID+"\t"+exonlength+"\t"+cdslength+"\t"+res.getGene().getReadCounts().size()+"\t"+res.getGene().getTranscriptReadCounts(thisTranscriptID)+"\t"+res.getGene().hasFlatPrior()+"\t%e\t"+res.getFlag()+"\t"+res.getNIterations()+"\t%f\t%e\t%e\t%e\t"+res.getGene().getTranscript(thisTranscriptID).getBiotype()+"\n", res.getGene().getPriors().get(thisTranscriptID), tmp_effectiveCount, density_exon, density_cds, tmp_finalLikelihood);
				}
			}			
			count++;
			if(Math.round((count*100.0)/totalGenes) > percent)
				percent = (int)Math.round((count*100.0)/totalGenes);
			Thunder.printProgressBar(percent);

		}
		//System.err.println("readMap.size() = "+_readMap.size());
		//System.err.println("\r                                                                                                               ");
		System.err.println("");

		out.flush();
		out.close();
		//System.err.println(this.getTime()+" Done.");
		Thunder.printLineErr("Done.");
	}


	/**
	 * 
	 * @param readMap
	 * @throws IOException
	 */
	public void outputSampleReadAssignments(File inputSAM) throws IOException{
		//PrintWriter out = new PrintWriter(new FileWriter(this._outputPrefix+".bam"));
		Thunder.printLineErr("Writing sampled read alignments to: "+this._outputPrefix+".sample.bam");

		SAMFileReader inputSam = new SAMFileReader(inputSAM);
		SAMFileWriterFactory fact=new SAMFileWriterFactory();
		fact.setCreateIndex(true);
		SAMFileWriter outputSam = fact.makeBAMWriter(inputSam.getFileHeader(), false, new File(this._outputPrefix+".sample.bam"));

		SAMRecord thisRecord;
		SAMRecordIterator it = inputSam.iterator();
		String thisTranscriptID,thisReadID; 
		while(it.hasNext()){
			thisRecord = it.next();
			thisTranscriptID = thisRecord.getReferenceName();
			thisReadID = thisRecord.getReadName();

			if(_readMap.containsKey(thisReadID)){
				if(thisTranscriptID.equals(_readMap.get(thisReadID))){
					outputSam.addAlignment(thisRecord);
					_readMap.remove(thisReadID);
				}
			}
		}
		outputSam.close();
		inputSam.close();
		Thunder.printLineErr("Done.");
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("annotationPath").hasArg().withDescription("Path to the GTF file containing the transcript/gene relationship").create(Thunder.OPT_PATH_ANNOTATION));
		options.addOption(OptionBuilder.withArgName("results.xprs").hasArg().withDescription("[optional] Text file containing RNA-seq transcript expression quantifications from eXpress").create("e"));
		options.addOption(OptionBuilder.withArgName("footprint alignments").hasArg().withDescription("SAM/BAM file containing alignments against the transcriptome").create("f"));
		options.addOption(OptionBuilder.withLongOpt("outputPrefix").withArgName("outputPath").hasArg().withDescription("Output prefix for results").create("o"));
		options.addOption(OptionBuilder.withLongOpt("sample").withDescription("Output sampling of read alignments based on EM transcript likelihoods").create("s"));
		options.addOption(OptionBuilder.withArgName("maxIterations_EM").hasArg().withDescription("Maximum number of EM iterations [default: 500]").create("N"));
		options.addOption(OptionBuilder.withArgName("double").hasArg().withDescription("Criteria for EM convergence [default: 0.0001]").create("c"));
		options.addOption(OptionBuilder.withDescription("Write all genes/transcripts, even those with no observed footprints").create("v"));
		options.addOption(OptionBuilder.withDescription("Count only footprints in the CDS (where available)").create("cds"));
		options.addOption(OptionBuilder.withDescription("Keep intermediate files").create("DEV"));
		return options;
	}


	/**
	 * 
	 * @param args
	 * @throws Exception
	 */
	public static void main(String[] args) throws Exception {

		/*args = new String[]{"IsoformEM_Footprints", 
						//"-f", "/Users/robk/Desktop/EM_TEST/TEST_NEW/D1-Footprint_1x50_gencode.v18.annotation_mappedReads.sorted.bam",
						"-f", "/Users/robk/Desktop/EM_TEST/TEST_NEW/ALL-Footprint_1x50_gencode.v18.annotation_mappedReads.sorted.bam",
						"-a", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v18.annotation.filtered.gtf",
						"-e", "/Users/robk/Box Sync/Work/HEK293_RNAseq/STAR_eXpress/A1-Total_2x75/results.xprs",
						"-o", "/Users/robk/Desktop/EM_TEST/TEST_NEW",
						"-DEV"};
		 */

		/*args = new String[]{"IsoformEM_Footprints",
				"-a","/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v21.annotation_noSelenocysteine.gtf",
				"-f","/Users/robk/WORK/YALE_offline/My Proteomics/HEK/RNA-seq_NEW/Footprinting/FP_original_1_Genome_Aligned.toTranscriptome.sorted.bam",
				"-N","1000",
				"--sample",
				"-o","/Users/robk/WORK/YALE_offline/My Proteomics/HEK/RNA-seq_NEW/Footprinting/FootprintEM/prior_none-EM_full-FP_D1-NEWNEW"};
		 */	

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption(Thunder.OPT_PATH_ANNOTATION) && cmdArgs.hasOption("f") && cmdArgs.hasOption("o")){
			System.err.println();


			FootprintEM engine = new FootprintEM(cmdArgs.getOptionValue("o"));

			//String output_path = cmdArgs.getOptionValue("o");
			//File output_priors = new File(output_path+"/priors.txt");
			//File output_alignments = new File(tmp_path+"/alignments.txt");


			// Create the output directory and temp dir if they don't already exist
			//engine.makeOutputDirectories(output_path);

			Thunder.printLineErr("Restrict to CDS aligned reads: "+cmdArgs.hasOption("cds"));
			
			// Read the transcripts in the GTF
			engine.readGTF(new File(cmdArgs.getOptionValue(Thunder.OPT_PATH_ANNOTATION)));

			// Read the eXpress transcript quantifications or choose flat priors
			if(cmdArgs.hasOption("e")){
				engine.readTranscriptExpressions(new File(cmdArgs.getOptionValue("e")));
			}

			// Convert the transcript expressions to priors and output
			engine.initialisePriors();

			// Parse the footprint alignments
			engine.readSAM(new File(cmdArgs.getOptionValue("f")), cmdArgs.hasOption("cds"));

			// Make temporary copies of the perl EM executables:
			//engine.makePerlExecutables(tmp_path);

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
			//engine.runEM("perl -I "+tmp_path+" "+tmp_path+"/main.foot.print.pl "+maxIterations+" "+output_priors+" "+output_alignments+" "+output_path);
			engine.doEM(maxIterations, convergenceDistance, outputAll);

			// Sample the reads?
			if(cmdArgs.hasOption("s")){
				engine.outputSampleReadAssignments(new File(cmdArgs.getOptionValue("f")));
			}

			// Clean up
			if(cmdArgs.hasOption("DEV")){
				//engine.removeFile(new File(tmp_path));
				engine.verbose = true;
			}

			System.err.println(engine.getTime()+" All Done!");

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" IsoformEM_Footprints", getCmdLineOptions());
			System.err.println();
		}


	}

	public String getTime(){
		return((new SimpleDateFormat("yyyy/MM/dd HH:mm:ss")).format(Calendar.getInstance().getTime()));
	}




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



