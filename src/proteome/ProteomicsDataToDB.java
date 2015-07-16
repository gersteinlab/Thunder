package proteome;

import java.sql.PreparedStatement;
import java.sql.SQLException;
import java.sql.Statement;

import objects.Isoform;
import objects.MS1_scan;
import objects.Peptide;
import objects.Spectra;
import database.DBConnect;

public class ProteomicsDataToDB {

	private Statement statement;
	private PreparedStatement preparedstatement;




	/**
	 * 
	 * @param tableName
	 * @param forceNewDB
	 * @throws SQLException
	 */
	public void createTables_mzXML(String tableName, boolean forceNewDB) throws SQLException{
		if(this.db.containsTable(tableName)==false || forceNewDB){
			statement.execute("DROP TABLE IF EXISTS "+tableName);
			statement.execute("CREATE TABLE "+tableName+"(" +
					"spectraID INT not null, " +
					"retentionTime VARCHAR(33) not null, " +
					"charge INT not null, " +
					"mz DOUBLE not null, " +
					"intensity DOUBLE not null, " +
					"primary key (spectraID)" +
					")");
		}
	}



	//TODO: Add modifications to spectra2Peptide!
	/**
	 * 
	 * @param tableName
	 * @param forceNewDB
	 * @throws SQLException
	 */
	public void createTable_XTandemResult(String tableName, boolean forceNewDB) throws SQLException{
		if(this.db.containsTable(tableName)==false || forceNewDB){
			statement.execute("DROP TABLE IF EXISTS "+tableName);
			statement.execute("CREATE TABLE "+tableName+"(" +
					"peptide_sequence varchar(100) not null, " +
					"spectraID int not null, " +
					"geneID varchar(30) not null, " +
					"isoformID varchar(30) not null, " +
					"frame int not null, " +
					"isReverse boolean not null, " +
					"start int not null, " +
					"end int not null, " +
					"expectation double not null, " +
					"mh double not null, " +
					"delta double not null, " +
					"hyperscore double not null, " +
					"nextscore double not null, " +
					"y_score double not null, " +
					"y_ions int not null, " +
					"b_score double not null, " +
					"b_ions int not null, " +
					"missed_cleavages int not null, " +
					"modifications varchar(200) null" +
					")");
		}
	}



	/**
	 * 
	 * @param thisScan
	 * @throws NumberFormatException
	 * @throws SQLException
	 */
	public void addSpectra_mzXML(MS1_scan thisScan) throws NumberFormatException, SQLException{
		this.preparedstatement.setInt(1, Integer.valueOf(thisScan.getAttribute(MS1_scan.ATTRIBUTE_SCAN_NUMBER)).intValue());
		this.preparedstatement.setString(2, thisScan.getAttribute(MS1_scan.ATTRIBUTE_RETENTION_TIME));
		this.preparedstatement.setInt(3, Integer.valueOf(thisScan.getPrecursorAttribute(MS1_scan.ATTRIBUTE_PRECURSOR_CHARGE)).intValue());
		this.preparedstatement.setDouble(4, Double.valueOf(thisScan.getPrecursorAttribute(MS1_scan.ATTRIBUTE_PRECURSOR_MZ)).doubleValue());
		this.preparedstatement.setDouble(5, Double.valueOf(thisScan.getPrecursorAttribute(MS1_scan.ATTRIBUTE_PRECURSOR_INTENSITY)).doubleValue());
		this.preparedstatement.addBatch();
		this.preparedstatement.executeBatch();
	}


	public void addSpectra(Spectra spectra, Peptide peptide, Isoform isoform) throws NumberFormatException, SQLException{

		
//		if(Integer.valueOf(spectra.getAttribute(Spectra.ATTRIBUTE_ID)).intValue() == 18156){
//			System.out.println("ADDING TO DB...");
//			for(int i=0;i<peptide.getModifications().size();i++){
//				for(int j=0;j<peptide.getModifications().get(i).getLength();j++){
//					System.out.println(peptide.getAttribute(Peptide.ATTRIBUTE_ID)+": "+peptide.getModifications().get(i).getLocalName(j)+"="+peptide.getModifications().get(i).getValue(j));
//				}
//				
//			}
//			System.out.println(peptide.getModificationsString());		
//		}
		
		//		String[] parsedID = parseIsoformID(isoform.getAttribute(Isoform.ATTRIBUTE_LABEL));
		String[] parsedID = parseIsoformID(isoform.getIsoformID());

		boolean isReverse = false;
		if(parsedID.length == 4)
			isReverse = true;

		this.preparedstatement.setString(1, peptide.getAttribute(Peptide.ATTRIBUTE_SEQ));
		this.preparedstatement.setInt(2, Integer.valueOf(spectra.getAttribute(Spectra.ATTRIBUTE_ID)).intValue());
		this.preparedstatement.setString(3, parsedID[2]);
		this.preparedstatement.setString(4, parsedID[0]);
		this.preparedstatement.setInt(5, Integer.valueOf(parsedID[1]).intValue());
		this.preparedstatement.setBoolean(6, isReverse);
		this.preparedstatement.setInt(7, Integer.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_START)).intValue());
		this.preparedstatement.setInt(8, Integer.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_END)).intValue());
		this.preparedstatement.setDouble(9, Double.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_EXPECT)).doubleValue());
		this.preparedstatement.setDouble(10, Double.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_MH)).doubleValue());
		this.preparedstatement.setDouble(11, Double.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_DELTA)).doubleValue());
		this.preparedstatement.setDouble(12, Double.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_HYPERSCORE)).doubleValue());
		this.preparedstatement.setDouble(13, Double.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_NEXTSCORE)).doubleValue());
		this.preparedstatement.setDouble(14, Double.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_Y_SCORE)).doubleValue());
		this.preparedstatement.setInt(15, Integer.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_Y_IONS)).intValue());
		this.preparedstatement.setDouble(16, Double.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_B_SCORE)).doubleValue());
		this.preparedstatement.setInt(17, Integer.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_B_IONS)).intValue());
		this.preparedstatement.setInt(18, Integer.valueOf(peptide.getAttribute(Peptide.ATTRIBUTE_MISSED_CLEAVAGES)).intValue());
		this.preparedstatement.setString(19, peptide.getModificationsString());
		this.preparedstatement.addBatch();
		this.preparedstatement.executeBatch();
	}



	/**
	 * Returns an id String[] containing {IsoformID, frame, GeneID, [reverse -- optional]}
	 * @param idString
	 * @return
	 */
	public String[] parseIsoformID(String idString){

		String[] tmpParsed = idString.split(" ");
		String frame = "1";
		String geneID = "";
		String transcriptID = "";
		//int transcriptIndex = 0;

		//		System.out.println("this.thisReferenceIDFormat = "+this.thisReferenceIDFormat);
		//		System.out.println(tmpParsed[0]);

		if(this.thisReferenceIDFormat == ProteomicsDataToDB.REFERENCE_FORMAT_RNASEQ  ||  this.thisReferenceIDFormat == ProteomicsDataToDB.REFERENCE_FORMAT_GENCODE  ||  
				this.thisReferenceIDFormat == ProteomicsDataToDB.REFERENCE_FORMAT_SWISSPROT){
			//		System.out.println(idString);
			// Find the part of the ID string containing the geneID
			for(int i=0;i<tmpParsed.length;i++){
				if(tmpParsed[i].startsWith("GN=")){
					geneID = tmpParsed[i].substring(3);
					break;
				}
			}
			if(tmpParsed.length > 1){
				if(this.thisReferenceIDFormat == ProteomicsDataToDB.REFERENCE_FORMAT_RNASEQ  ||  this.thisReferenceIDFormat == ProteomicsDataToDB.REFERENCE_FORMAT_GENCODE){
//					frame = tmpParsed[0].split("_")[2];
					frame = tmpParsed[0].substring(tmpParsed[0].length()-1);
					transcriptID = tmpParsed[0].substring(0, tmpParsed[0].length()-2); 
				}
//				else if(this.thisReferenceIDFormat == ProteomicsDataToDB.REFERENCE_FORMAT_GENCODE){
//					frame = tmpParsed[0].split("_")[1];
//					transcriptID = tmpParsed[0].substring(0, tmpParsed[0].length()-2); 
//				}
				else if(this.thisReferenceIDFormat == ProteomicsDataToDB.REFERENCE_FORMAT_SWISSPROT){
					transcriptID = tmpParsed[0];
					if(geneID.endsWith(":reversed"))
						geneID = geneID.replaceAll(":reversed", "");
				}
			}else{
				// For spectra mapped to the same reference as eXpress
				tmpParsed = idString.split(":");
				if(tmpParsed[0].split("\\|").length > 1){
					// this is a contaminant entry from cRAP
					frame = "1";
					transcriptID = tmpParsed[0];
				}else{
					frame = tmpParsed[0].substring(tmpParsed[0].length()-1);
					transcriptID = tmpParsed[0].substring(0, tmpParsed[0].length()-2);
				}
			}
		}
		else if(this.thisReferenceIDFormat == ProteomicsDataToDB.REFERENCE_FORMAT_TRINITY){
			if(tmpParsed.length > 1){
				tmpParsed = tmpParsed[0].split("_");
				geneID = tmpParsed[0];
				transcriptID = tmpParsed[2];
				frame = tmpParsed[3];
			}else{
				geneID = tmpParsed[0];
				transcriptID = geneID;
			}
		}

		//System.out.println(idString+" -- geneID:"+geneID+"\ttranscriptID:"+transcriptID+"\tframe:"+frame);

		if(idString.endsWith(":reversed"))
			tmpParsed = new String[]{transcriptID, frame, geneID, "reverse"};
		else{
			tmpParsed = new String[]{transcriptID, frame, geneID};
		}

		return tmpParsed;
	}



	/**
	 * 
	 * @throws SQLException
	 */
	public void closeConnection() throws SQLException{
		this.statement.close();
		this.preparedstatement.close();
		this.db.closeConnection();
	}



	/**
	 * 
	 * @return
	 */
	public Statement getDBStatement(){
		return this.statement;
	}


	public static final int REFERENCE_FORMAT_SWISSPROT = 1;
	public static final int REFERENCE_FORMAT_RNASEQ = 2;
	public static final int REFERENCE_FORMAT_GENCODE = 3;
	public static final int REFERENCE_FORMAT_TRINITY = 4;
	private int thisReferenceIDFormat;

	private DBConnect db;
	private boolean isNewDB;

	public boolean isNewDB(){ return isNewDB; }


	/**
	 * For reading mapped XTandem data to the DB
	 * 
	 * @param dbURL
	 * @param tableName
	 * @param referenceIDFormat
	 * @param forceNewDB
	 * @throws Exception
	 */
	public ProteomicsDataToDB(DBConnect db, String tableName, int referenceIDFormat, boolean forceNewDB) throws Exception{
		//		this.connection = DBConnect.getConnection_H2(dbURL, "sa", "");
		this.db = db;
		this.statement = db.createStatement();

		this.isNewDB = false;
		if(this.db.containsTable(tableName) == false){ this.isNewDB = true; }

		createTable_XTandemResult(tableName, forceNewDB);
		this.preparedstatement = db.getPreparedStatement("INSERT INTO "+tableName+" VALUES (?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?, ?)");
		this.thisReferenceIDFormat = referenceIDFormat;
	}



	/**
	 * For reading unmapped MZXML data to the DB
	 * 
	 * @param dbURL
	 * @param tableName
	 * @param forceNewDB
	 * @throws Exception
	 */
	public ProteomicsDataToDB(DBConnect db, String tableName, boolean forceNewDB) throws Exception{
		this.db = db;
		this.statement = db.createStatement();

		this.isNewDB = false;
		if(this.db.containsTable(tableName) == false){ this.isNewDB = true; }
		createTables_mzXML(tableName, forceNewDB);
		this.preparedstatement = db.getPreparedStatement("INSERT INTO "+tableName+" VALUES (?, ?, ?, ?, ?)");
	}


}
