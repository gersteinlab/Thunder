package proteome;

import java.io.File;
import java.sql.SQLException;

import javax.xml.parsers.SAXParser;
import javax.xml.parsers.SAXParserFactory;

import main.Thunder;

import objects.MS1_scan;

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

public class ReadmzXML {


	/**
	 * Reads the contents of an mzXML file to the database  
	 * @param filePath
	 * @param dbURL
	 * @param tableName
	 * @param forceNewDB
	 */
	public static void readFile(String filePath, String dbPath, String tableName, boolean forceNewDB){

		try {
			DBConnect db = new DBConnect_SQLite(dbPath);

			ProteomicsDataToDB init = new ProteomicsDataToDB(db, tableName, forceNewDB);
			if(init.isNewDB()==true || forceNewDB){
				SAXParserFactory factory = SAXParserFactory.newInstance();
				SAXParser saxParser = factory.newSAXParser();
				MySAXHandler_mzXML handler = new MySAXHandler_mzXML();
				db.setAutoCommit(false);
				handler.attachDB(init);
				saxParser.parse(filePath, handler);
				db.commit();
				db.setAutoCommit(true);
				handler.getDataset().closeConnection();
			}
			db.closeConnection();
		} catch (Exception e) {
			e.printStackTrace();
		}
	}



	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("spectraPath").hasArg().withDescription("path to the mzXML file containing the spectra").create(Thunder.OPT_PATH_SPECTRA_MZXML));
		options.addOption(OptionBuilder.withArgName("databasePath").hasArg().withDescription("add spectra to the sample database at this path\n(e.g. /path/to/db)").create(Thunder.OPT_PATH_DB_SAMPLE));
		options.addOption(OptionBuilder.withArgName("tableName").hasArg().withDescription("name of the table to which to add the annotation data").create(Thunder.OPT_DB_TABLE_NAME));
		options.addOption(new Option(Thunder.OPT_CHOICE_DB_FORCE_REFRESH, "reads and adds the annotation to the database, regardless of whether the data already exists"));
		return options;
	}



	public static void main(String[] args) throws Exception {
		
		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		boolean forceNewDB = false;
		if(cmdArgs.hasOption(Thunder.OPT_CHOICE_DB_FORCE_REFRESH)){
			forceNewDB = true;
		}
		
		if(cmdArgs.hasOption(Thunder.OPT_PATH_SPECTRA_MZXML) && cmdArgs.hasOption(Thunder.OPT_PATH_DB_SAMPLE)){
			File inputFile = new File(cmdArgs.getOptionValue(Thunder.OPT_PATH_SPECTRA_MZXML));
			String tableName = inputFile.getName().replaceAll("\\.|-", "_");
			if(cmdArgs.hasOption(Thunder.OPT_DB_TABLE_NAME)){
				tableName = cmdArgs.getOptionValue(Thunder.OPT_DB_TABLE_NAME);
			}
			
			readFile(inputFile.getAbsolutePath(), cmdArgs.getOptionValue(Thunder.OPT_PATH_DB_SAMPLE), tableName, forceNewDB);
			
		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ReadMzXML", getCmdLineOptions());
			System.out.println();
		}
		
	}

}




/**
 * Defines how to deal with each line of XML from the X!Tandem output file
 * 
 */
class MySAXHandler_mzXML extends DefaultHandler{

	private ProteomicsDataToDB dataset;
	public void attachDB(ProteomicsDataToDB db) throws Exception{
		this.dataset = db;
	}

	private MS1_scan thisScan;
	private boolean parsingPrecursorInfo = false;
	private boolean parsingPeaks = false;

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
					this.dataset.addSpectra_mzXML(thisScan);
				} catch (NumberFormatException e) {
					e.printStackTrace();
				} catch (SQLException e) {
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

	/**
	 * Return the dataset object holding the parsed XML data
	 * @return
	 */
	public ProteomicsDataToDB getDataset(){
		return this.dataset;
	}
}



