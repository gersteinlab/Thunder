package proteome;

import java.io.File;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import utils.IO_utils;

public class ReadXTandem {

	private static SpectraAlignmentEngine createSpectraEngine(File tandemResults) throws Exception{
		return ReadXTandem.createSpectraEngine(tandemResults, "");
	}

	private static SpectraAlignmentEngine createSpectraEngine(File tandemResults, String spetraIDprefix) throws Exception{
		SpectraAlignmentEngine spectraEngine = new SpectraAlignmentEngine(spetraIDprefix);

		try {
			SAXParserFactory factory = SAXParserFactory.newInstance();
			SAXParser saxParser = factory.newSAXParser();
			SAXHandler_XTandem handler = new SAXHandler_XTandem();

			handler.attachSpectraAlignmentEngine(spectraEngine);
			//Thunder.printLineErr("Reading X!Tandem results: ");
			saxParser.parse(tandemResults, handler);

		}catch (Exception e) {
			e.printStackTrace();
		}finally{

		}

		System.err.println();
		IO_utils.printLineErr("  N spectra:\t"+spectraEngine.countSpectra());
		IO_utils.printLineErr("  N alignments:\t"+spectraEngine.getNumberOfAlignments());
		//System.out.println("N peptides:\t"+spectraEngine.countPeptides());
		//System.out.println("N proteins:\t"+spectraEngine.countProteins());
		//System.out.println();

		return(spectraEngine);
	}


	/**
	 * 
	 * @param tandemResults
	 * @param fdrMax
	 * @param removeCRAPome
	 * @return
	 * @throws Exception
	 */
	public static SpectraAlignmentEngine processTandemXML(File tandemResults, double fdrMax, boolean removeCRAPome) throws Exception{
		IO_utils.printLineErr("Reading spectra alignments: "+tandemResults.getAbsolutePath());
		SpectraAlignmentEngine spectraEngine = ReadXTandem.createSpectraEngine(tandemResults);
		spectraEngine.processAlignments(fdrMax, new File(tandemResults.getAbsolutePath()+".processedAlignments.txt"), removeCRAPome);
		return(spectraEngine);
	}

	/**
	 * 
	 * @param tandemResults
	 * @param fdrMax
	 * @param removeCRAPome
	 * @return
	 * @throws Exception
	 */
	public static SpectraAlignmentEngine processTandemXML(File[] tandemResults, String outputPrefix, double fdrMax, boolean removeCRAPome) throws Exception{
		if(tandemResults.length >= 1){
			IO_utils.printLineErr("Reading spectra alignments: "+tandemResults[0].getAbsolutePath());
			SpectraAlignmentEngine spectraEngineMAIN = ReadXTandem.createSpectraEngine(tandemResults[0], tandemResults[0].getName()+"_");

			SpectraAlignmentEngine[] spectraEngines = new SpectraAlignmentEngine[tandemResults.length-1];
			for(int i=1;i<tandemResults.length;i++){
				IO_utils.printLineErr("Reading spectra alignments: "+tandemResults[i].getAbsolutePath());
				spectraEngines[i-1] = ReadXTandem.createSpectraEngine(tandemResults[i], tandemResults[i].getName()+"_");
				
				// merge alignments
				spectraEngineMAIN.mergeAlignments(spectraEngines[i-1].getAlignments());
			}

			// NO FDR
			//String tmpFileName = tandemResults[0].getParent()+"/MERGED.processedAlignments_FDR{none}.txt";
			String tmpFileName = outputPrefix+".processedAlignments_FDR{none}.txt";
			IO_utils.printLineErr("Writing merged, unfiltered, spectra alignments to: "+tmpFileName);
			spectraEngineMAIN.writeAlignments(new File(tmpFileName));
			
			// process FDR etc and return
			//tmpFileName = tandemResults[0].getParent()+"/MERGED.processedAlignments_FDR{"+fdrMax+"}.txt";
			tmpFileName = outputPrefix+".processedAlignments_FDR{"+fdrMax+"}.txt";
			IO_utils.printLineErr("Writing merged FDR filtered spectra alignments to: "+tmpFileName);
			spectraEngineMAIN.processAlignments(fdrMax, new File(tmpFileName), removeCRAPome);
			return(spectraEngineMAIN);
		}
		else{
			return null;
		}
	}



	public static String parseGeneID(String idString){
		String[] tmpParsed;
		tmpParsed = idString.split(" ");
		for(int i=0;i<tmpParsed.length;i++){
			if(tmpParsed[i].startsWith("GN=")){
				idString = tmpParsed[i].substring(3);
				break;
			}
		}
		return idString;
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("tandemPath").hasArg().withDescription("path to the X!Tandem XML file containing the spectra identifications").create(Thunder.OPT_PATH_SPECTRA_TANDEM));
		//options.addOption(OptionBuilder.withArgName("databasePath").hasArg().withDescription("add spectra to the sample database at this path\n(e.g. /path/to/db)").create(Thunder.OPT_PATH_DB_SAMPLE));
		//options.addOption(OptionBuilder.withArgName("tableName").hasArg().withDescription("name of the table to which to add the annotation data").create(Thunder.OPT_DB_TABLE_NAME));
		//options.addOption(new Option(Thunder.OPT_CHOICE_DB_FORCE_REFRESH, "reads and adds the annotation to the database, regardless of whether the data already exists"));
		//options.addOption(OptionBuilder.withLongOpt("IDcuffmerge").withDescription("gene ID format used to map the MS/MS spectra is derived from cuffmerge [default]").create(Thunder.OPT_FORMAT_CUFFMERGE));
		//options.addOption(OptionBuilder.withLongOpt("IDswissprot").withDescription("gene ID format used to map the MS/MS spectra is derived from SwissProt").create(Thunder.OPT_FORMAT_SWISSPROT));
		//options.addOption(OptionBuilder.withLongOpt("IDgencode").withDescription("gene ID format used to map the MS/MS spectra is derived from Gencode").create(Thunder.OPT_FORMAT_GENCODE));
		//options.addOption(OptionBuilder.withLongOpt("IDtrinity").withDescription("gene ID format used to map the MS/MS spectra is derived from De-novo RNA-seq assembly using Trinity").create(Thunder.OPT_FORMAT_TRINITY));
		options.addOption(OptionBuilder.withArgName("double").hasArg().withDescription("[default: 0.05] Based on hits to the reverse decoy DB, define an acceptable false discovery rate (FDR)").create("fdr"));
		options.addOption(OptionBuilder.withLongOpt("swissprot").withDescription("these alignments are against the swissprot reference").create("s"));
		return options;
	}


	public static void main(String[] args) throws Exception {

		//args = new String[]{"-"+Thunder.OPT_PATH_SPECTRA_TANDEM, "/Users/robk/WORK/YALE_offline/My Proteomics/HEK/RNA-seq_NEW/Proteomics/HEK.C8.Orbi.gencode21.TandemOutput.xml"};
		args = new String[]{"-s", "-"+Thunder.OPT_PATH_SPECTRA_TANDEM, "/Users/robk/WORK/YALE_offline/My Proteomics/HEK/RNA-seq_NEW/Proteomics/HEK.C8.Orbi.swissprot.TandemOutput.xml"};
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		if(cmdArgs.hasOption(Thunder.OPT_PATH_SPECTRA_TANDEM)){

			//
			double fdrMax = 0.05;
			if(cmdArgs.hasOption("fdr"))
				fdrMax = Double.valueOf(cmdArgs.getOptionValue("fdr")).doubleValue();

			boolean removeCRAPome = true;
			if(cmdArgs.hasOption("s"))
				removeCRAPome = false;

			processTandemXML(new File(cmdArgs.getOptionValue(Thunder.OPT_PATH_SPECTRA_TANDEM)), fdrMax, removeCRAPome);	

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ReadTandem", getCmdLineOptions());
			System.out.println();
		}
	}


	//public static void main(String[] args) throws Exception {
	//	main(Thunder.parseArgs(args, getCmdLineOptions()));
	//}

}
