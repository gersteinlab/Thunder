package proteome;

import java.io.File;
import java.sql.SQLException;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import main.Thunder;
import objects.Isoform;
import objects.Peptide;
import objects.Spectra;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;
import org.xml.sax.Attributes;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.DefaultHandler;

import database.DBConnect;
import database.DBConnect_SQLite;

public class ReadXTandem_toDB {


	public static void readTandemXML(File tandemResults, int referenceIDFormat, String dbPath, String tableName, boolean forceNewDB) throws Exception{
		DBConnect db = new DBConnect_SQLite(dbPath);
		try {

			ProteomicsDataToDB init = new ProteomicsDataToDB(db, tableName, referenceIDFormat, forceNewDB);

			if(init.isNewDB()==true || forceNewDB){

				SAXParserFactory factory = SAXParserFactory.newInstance();
				SAXParser saxParser = factory.newSAXParser();
				SAXHandler handler = new SAXHandler();
				db.setAutoCommit(false);
				handler.attachDB(init);
				saxParser.parse(tandemResults, handler);
				db.commit();
				db.setAutoCommit(true);
				handler.getDataset().closeConnection();

			}else{
				//				this.dataset = init;
			}
		}catch (Exception e) {
			e.printStackTrace();
		}finally{
			db.closeConnection();
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
		options.addOption(OptionBuilder.withArgName("databasePath").hasArg().withDescription("add spectra to the sample database at this path\n(e.g. /path/to/db)").create(Thunder.OPT_PATH_DB_SAMPLE));
		options.addOption(OptionBuilder.withArgName("tableName").hasArg().withDescription("name of the table to which to add the annotation data").create(Thunder.OPT_DB_TABLE_NAME));
		options.addOption(new Option(Thunder.OPT_CHOICE_DB_FORCE_REFRESH, "reads and adds the annotation to the database, regardless of whether the data already exists"));
		options.addOption(OptionBuilder.withLongOpt("IDcuffmerge").withDescription("gene ID format used to map the MS/MS spectra is derived from cuffmerge [default]").create(Thunder.OPT_FORMAT_CUFFMERGE));
		options.addOption(OptionBuilder.withLongOpt("IDswissprot").withDescription("gene ID format used to map the MS/MS spectra is derived from SwissProt").create(Thunder.OPT_FORMAT_SWISSPROT));
		options.addOption(OptionBuilder.withLongOpt("IDgencode").withDescription("gene ID format used to map the MS/MS spectra is derived from Gencode").create(Thunder.OPT_FORMAT_GENCODE));
		options.addOption(OptionBuilder.withLongOpt("IDtrinity").withDescription("gene ID format used to map the MS/MS spectra is derived from De-novo RNA-seq assembly using Trinity").create(Thunder.OPT_FORMAT_TRINITY));
		return options;
	}


	public static void main(CommandLine cmdArgs) throws Exception {
		/* Force refresh this table? */
		boolean forceNewDB = false;
		if(cmdArgs.hasOption(Thunder.OPT_CHOICE_DB_FORCE_REFRESH)){
			forceNewDB = true;
		}

		/* Get reference geneID format */
		//int referenceType = ProteomicsDataToDB.REFERENCE_FORMAT_RNASEQ;
		int referenceType = ProteomicsDataToDB.REFERENCE_FORMAT_GENCODE;
		if(cmdArgs.hasOption(Thunder.OPT_FORMAT_SWISSPROT)){
			referenceType = ProteomicsDataToDB.REFERENCE_FORMAT_SWISSPROT;
		}else if(cmdArgs.hasOption(Thunder.OPT_FORMAT_GENCODE)){
			referenceType = ProteomicsDataToDB.REFERENCE_FORMAT_GENCODE;
		}else if(cmdArgs.hasOption(Thunder.OPT_FORMAT_TRINITY)){
			referenceType = ProteomicsDataToDB.REFERENCE_FORMAT_TRINITY;
		}


		if(cmdArgs.hasOption(Thunder.OPT_PATH_DB_SAMPLE) && cmdArgs.hasOption(Thunder.OPT_PATH_SPECTRA_TANDEM)){

			File tandemOutput = new File(cmdArgs.getOptionValue(Thunder.OPT_PATH_SPECTRA_TANDEM));
			String tableName = tandemOutput.getName().replaceAll("\\.|-", "_");
			if(cmdArgs.hasOption(Thunder.OPT_DB_TABLE_NAME)){
				tableName = cmdArgs.getOptionValue(Thunder.OPT_DB_TABLE_NAME);
			}

			readTandemXML(tandemOutput, referenceType, cmdArgs.getOptionValue(Thunder.OPT_PATH_DB_SAMPLE), tableName, forceNewDB);	
		
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ReadTandem", getCmdLineOptions());
			System.out.println();
		}
	}

	
	public static void main(String[] args) throws Exception {
		main(Thunder.parseArgs(args, getCmdLineOptions()));
	}

}



/**
 * Defines how to deal with each line of XML from the X!Tandem output file
 * 
 */
class SAXHandler extends DefaultHandler{

	private ProteomicsDataToDB dataset;
	public void attachDB(ProteomicsDataToDB db) throws Exception{
		this.dataset = db;
	}

	private Spectra thisSpectra;
	private Peptide thisPeptide;
	private Isoform thisIsoform;
	//private boolean parsingPeptide = false;
	//private boolean parsingNote = false;

	private int spectraCount = 0;
	private void updateSpectraCount(){
		this.spectraCount ++;
		if(this.spectraCount % 100000 == 0){ System.out.println(this.spectraCount); }
		else if(this.spectraCount % 10000 == 0){ System.out.print(this.spectraCount); }
		else if(this.spectraCount % 1000 == 0){ System.out.print("."); }
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
			//parsingNote = true;
		}else if(qName.equals("peptide")){
			//parsingPeptide = true;
		}else if(qName.equals("domain")){
			//parsingPeptide = false;	
			this.thisPeptide = new Peptide(attributes);

		}else if(qName.equals("aa")){
			thisPeptide.addModification(attributes);
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
				dataset.addSpectra(thisSpectra, thisPeptide, thisIsoform);
			} catch (NumberFormatException e) {
				e.printStackTrace();
			} catch (SQLException e) {
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
			//thisIsoform.setIsoformID(new String(ch, start, length));
		//}
		//else if(parsingPeptide){
		//	thisIsoform.appendIsoformSequence((new String(ch, start, length)).replaceAll(" |\t|\n", "").trim());
		//}
	}

	/**
	 * Return the dataset object holding the parsed XML data
	 * @return
	 */
	public ProteomicsDataToDB getDataset(){
		return this.dataset;
	}

}


