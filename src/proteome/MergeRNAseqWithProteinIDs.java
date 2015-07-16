package proteome;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.Writer;
import java.sql.ResultSet;
import java.sql.Statement;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.OptionBuilder;
import org.apache.commons.cli.Options;

import database.DBConnect;
import database.DBConnect_SQLite;

public class MergeRNAseqWithProteinIDs {

	private int referenceType = 0;


	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		options.addOption(OptionBuilder.withArgName("annotationPath").withLongOpt("pathToAnnotationDB").hasArg().withDescription("Path of the database containing the annotation used to assemble transcripts(e.g. /path/to/db)").create(Thunder.OPT_PATH_DB_ANNOTATION));
		options.addOption(OptionBuilder.withArgName("sampleDBPath").withLongOpt("pathToSampleDB").hasArg().withDescription("add spectra to the sample database at this path\n(e.g. /path/to/db)").create(Thunder.OPT_PATH_DB_SAMPLE));
		options.addOption(OptionBuilder.withLongOpt("spectraFile").hasArg().withDescription("MGF file containing MS/MS spectra").create(Thunder.OPT_PATH_SPECTRA_TANDEM));
		options.addOption(OptionBuilder.withLongOpt("tableName").hasArg().withDescription("name of the table to which to add the annotation data").create(Thunder.OPT_DB_TABLE_NAME));
		options.addOption(OptionBuilder.withLongOpt("forceRefreshDB").withDescription("reads and adds the spectra data in this XML to the database, regardless of whether data from a file of the same name already exists").create(Thunder.OPT_CHOICE_DB_FORCE_REFRESH));
		options.addOption(OptionBuilder.withLongOpt("IDcuffmerge").withDescription("gene ID format used to map the MS/MS spectra is derived from cuffmerge [default]").create(Thunder.OPT_FORMAT_CUFFMERGE));
		options.addOption(OptionBuilder.withLongOpt("IDswissprot").withDescription("gene ID format used to map the MS/MS spectra is derived from SwissProt").create(Thunder.OPT_FORMAT_SWISSPROT));
		options.addOption(OptionBuilder.withLongOpt("IDgencode").withDescription("gene ID format used to map the MS/MS spectra is derived from Gencode").create(Thunder.OPT_FORMAT_GENCODE));
		//		options.addOption(new Option(Thunder.OPT_FORMAT_CUFFMERGE, "gene ID format used to map the MS/MS spectra is derived from cuffmerge [default]"));
		//		options.addOption(new Option(Thunder.OPT_FORMAT_SWISSPROT, "gene ID format used to map the MS/MS spectra is derived from SwissProt"));
		//		options.addOption(new Option(Thunder.OPT_FORMAT_GENCODE, "gene ID format used to map the MS/MS spectra is derived from Gencode"));
		return options;
	}




	/**
	 * Constructor, reads data
	 * @param paths Array containing locations of the 3 required data files {annotationGTF, cufflinksGTF, MSMSspectra}
	 * @throws Exception 
	 */
	public MergeRNAseqWithProteinIDs(String[] args) throws Exception{

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());

		/* Force refresh this table? */
		boolean forceNewDB = false;
		if(cmdArgs.hasOption(Thunder.OPT_CHOICE_DB_FORCE_REFRESH)){
			forceNewDB = true;
		}
		
		/* Get reference geneID format */
		this.referenceType = ProteomicsDataToDB.REFERENCE_FORMAT_RNASEQ;
		if(cmdArgs.hasOption(Thunder.OPT_FORMAT_SWISSPROT)){
			this.referenceType = ProteomicsDataToDB.REFERENCE_FORMAT_SWISSPROT;
		}else if(cmdArgs.hasOption(Thunder.OPT_FORMAT_GENCODE)){
			this.referenceType = ProteomicsDataToDB.REFERENCE_FORMAT_GENCODE;
		}

		if(cmdArgs.hasOption(Thunder.OPT_PATH_DB_SAMPLE) && cmdArgs.hasOption(Thunder.OPT_PATH_SPECTRA_TANDEM) && cmdArgs.hasOption(Thunder.OPT_PATH_DB_ANNOTATION)){
			
			// read the X!Tandem results to the DB
			ReadXTandem_toDB.main(cmdArgs);

			// 
			File tandemOutput = new File(cmdArgs.getOptionValue(Thunder.OPT_PATH_SPECTRA_TANDEM));
			String tableName = tandemOutput.getName().replaceAll("\\.|-", "_");
			if(cmdArgs.hasOption(Thunder.OPT_DB_TABLE_NAME)){
				tableName = cmdArgs.getOptionValue(Thunder.OPT_DB_TABLE_NAME);
			}

			// do the intersections with the annotation
			processSpectra(tandemOutput, cmdArgs.getOptionValue(Thunder.OPT_PATH_DB_SAMPLE), tableName, this.referenceType, forceNewDB);

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(Thunder.THUNDER_EXE_COMMAND+" ReadTandem", getCmdLineOptions());
			System.out.println();
		}
	}



	public void processSpectra(File spectra, String dbPath, String tableName, int geneIDformat, boolean forceNewDB) throws Exception{

		DBConnect db = new DBConnect_SQLite(dbPath);
		Statement st = db.createStatement();

		//
		// Keep only spectra / gene assignments for simplicity
		//
		System.out.print("Collapsing spectra by gene assignment...");
		st.executeUpdate("DROP TABLE IF EXISTS tmp_spectraGeneAssignments");
		st.executeUpdate("CREATE TABLE tmp_spectraGeneAssignments " +
				"AS SELECT DISTINCT peptide_Sequence,spectraID,geneID,isReverse,expectation,COUNT(isoformID) " +
				"FROM "+tableName+" " +
				"GROUP BY spectraID,peptide_Sequence,geneID,isReverse,expectation " +
				"ORDER BY spectraID");
		System.out.print("Indexing spectra...");
		st.executeUpdate("CREATE INDEX tmpIndex_spectraGeneAssignments_geneID ON tmp_spectraGeneAssignments (geneID)");
		System.out.println("Done");

		ResultSet rs;

		if(this.referenceType == ProteomicsDataToDB.REFERENCE_FORMAT_SWISSPROT){

			//
			// Get annotated spectra mapping to KNOWN proteome
			//
			System.out.print("Annotating spectra assigned to SwissProt proteins...");
			rs = st.executeQuery("SELECT * FROM tmp_spectraGeneAssignments");

			// Write the annotated spectra
			System.out.print("Writing annotated spectra...");
			Writer out = new BufferedWriter(new FileWriter(spectra.getAbsolutePath()+".peptideSummary.txt"));

			int nCols = rs.getMetaData().getColumnCount();
			String colNames = "";
			for(int i=1;i<=nCols;i++)
				colNames += rs.getMetaData().getColumnName(i)+"\t";
			out.write(colNames+"\n");

			String toOut = "";
			Object tmp;
			while(rs.next()){
				toOut = "";
				for(int i=1;i<=nCols;i++){
					tmp = rs.getObject(i);
					if(tmp != null)
						toOut += tmp.toString()+"\t";
					else
						toOut += "null\t";
				}
				out.write(toOut+"\n");
			}
			out.flush();
			System.out.println("Done");
			out.close();
			rs.close();

		}else if(this.referenceType == ProteomicsDataToDB.REFERENCE_FORMAT_GENCODE){

			//
			// Get annotated spectra mapping to KNOWN proteome
			//
			System.out.print("Annotating spectra assigned to GENCODE transcripts...");
			rs = st.executeQuery("SELECT * FROM tmp_spectraGeneAssignments");

			// Write the annotated spectra
			System.out.print("Writing annotated spectra...");
			Writer out = new BufferedWriter(new FileWriter(spectra.getAbsolutePath()+".peptideSummary.txt"));

			int nCols = rs.getMetaData().getColumnCount();
			String colNames = "";
			for(int i=1;i<=nCols;i++)
				colNames += rs.getMetaData().getColumnName(i)+"\t";
			out.write(colNames+"\n");

			String toOut = "";
			Object tmp;
			while(rs.next()){
				toOut = "";
				for(int i=1;i<=nCols;i++){
					tmp = rs.getObject(i);
					if(tmp != null)
						toOut += tmp.toString()+"\t";
					else
						toOut += "null\t";
				}
				out.write(toOut+"\n");
			}
			out.flush();
			System.out.println("Done");
			out.close();
			rs.close();

		}else if(this.referenceType == ProteomicsDataToDB.REFERENCE_FORMAT_RNASEQ){
			//
			// Create a merged annotation from the gencode data and the cufflinks build
			//
			System.out.print("Merging gene annotation and transcript builds...");
			st.executeUpdate("DROP TABLE IF EXISTS mergedAnnotation_GeneLevel");
			st.executeUpdate("CREATE TABLE mergedAnnotation_GeneLevel(GeneID_cuff varchar(30),chromosome varchar(20),start int,stop int," +
					"strand char(1),geneName varchar(30),transcriptID_cuffbest varchar(30),transcriptID_gencode varchar(30),chromosome_gencode varchar(20),start_gencode int," +
					"stop_gencode int,strand_gencode char(1),geneID_ensembl varchar(30),geneName_ensembl varchar(30),geneType varchar(30),geneStatus varchar(30)," +
					"transcriptName varchar(30),transcriptType varchar(50),transcriptStatus varchar(30),transcriptID_Havana varchar(30),geneID_Havana varchar(30)," +
					"primary key(GeneID_cuff)) " +
					"AS SELECT * FROM CUFFMERGEANNOTATION_GENELEVEL " +
					"JOIN GENCODE14ANNOTATION " +
					"ON CUFFMERGEANNOTATION_GENELEVEL.transcriptID_ensembl=GENCODE14ANNOTATION.TRANSCRIPTID " +
					"WHERE CUFFMERGEANNOTATION_GENELEVEL.transcriptID_ensembl LIKE 'ENST%'");
			// Some cufflinks genes have no 'real' annotated transcripts (such as small non-coding RNAs etc) so re-do search based on ensembl gene ID:
			st.executeUpdate("MERGE INTO mergedAnnotation_GeneLevel KEY(GeneID_cuff) " +
					"SELECT * FROM CUFFMERGEANNOTATION_GENELEVEL " +
					"JOIN GENCODE14ANNOTATION ON CUFFMERGEANNOTATION_GENELEVEL.transcriptID_ensembl=GENCODE14ANNOTATION.GENEID " +
					"WHERE CUFFMERGEANNOTATION_GENELEVEL.transcriptID_ensembl LIKE 'ENSG%'");
			System.out.print("Indexing merged annotation...");
			st.executeUpdate("CREATE INDEX mergedIndexCuffGenes ON mergedAnnotation_GeneLevel (geneID_cuff)");
			System.out.println("Done");


			//
			// Get annotated spectra mapping to KNOWN transcriptome
			//
			System.out.print("Annotating spectra assigned to KNOWN transcripts...");
			rs = st.executeQuery("SELECT * FROM tmp_spectraGeneAssignments " +
					"JOIN mergedAnnotation_GeneLevel " +
					"ON tmp_spectraGeneAssignments.geneID=mergedAnnotation_GeneLevel.geneID_cuff");

			// Write the annotated spectra
			System.out.print("Writing annotated spectra...");
			Writer out = new BufferedWriter(new FileWriter(spectra.getAbsolutePath()+".peptideSummary.txt"));

			int nCols = rs.getMetaData().getColumnCount();
			String colNames = "";
			for(int i=1;i<=nCols;i++)
				colNames += rs.getMetaData().getColumnName(i)+"\t";
			out.write(colNames+"\n");

			String toOut = "";
			Object tmp;
			while(rs.next()){
				toOut = "";
				for(int i=1;i<=nCols;i++){
					tmp = rs.getObject(i);
					if(tmp != null)
						toOut += tmp.toString()+"\t";
					else
						toOut += "null\t";
				}
				out.write(toOut+"\n");
			}
			out.flush();
			System.out.println("Done");



			//
			// Get annotated spectra mapping to NOVEL transcriptome
			//
			System.out.print("Annotating spectra assigned to NOVEL transcripts...");
			rs = st.executeQuery("SELECT * FROM tmp_spectraGeneAssignments " +
					"JOIN CUFFMERGEANNOTATION_GENELEVEL " +
					"ON tmp_spectraGeneAssignments.geneID=CUFFMERGEANNOTATION_GENELEVEL.geneID " +
					"WHERE CUFFMERGEANNOTATION_GENELEVEL.TRANSCRIPTID_ENSEMBL LIKE 'CUFF%' " +
					"AND CUFFMERGEANNOTATION_GENELEVEL.GENENAME IS NULL");
			// Write spectra mapped to novel sequences
			System.out.print("Writing annotated spectra...");
			int nCols_novel = rs.getMetaData().getColumnCount();
			//colNames_novel = "";
			//for(int i=1;i<=nCols_novel;i++)
			//	colNames_novel += rs.getMetaData().getColumnName(i)+"\t";
			//System.out.println("\n\n"+colNames_novel+"\n");
			while(rs.next()){
				toOut = "";
				for(int i=1;i<=nCols_novel;i++){
					tmp = rs.getObject(i);
					if(tmp != null)
						toOut += tmp.toString()+"\t";
					else
						toOut += "null\t";
				}
				out.write(toOut+"\t\t\t\t\t\t\t\t\t\t\t\t\t\n");
			}
			out.flush();
			System.out.println("Done");



			//
			// Get annotated spectra mapping to DECOY transcriptome
			//
			System.out.print("Annotating spectra assigned to DECOY transcripts...");
			rs = st.executeQuery("SELECT * FROM tmp_spectraGeneAssignments WHERE isReverse = TRUE");
			// Write spectra mapped to the reverse sequences
			System.out.print("Writing annotated spectra...");
			int nCols_decoy = rs.getMetaData().getColumnCount();
			//		colNames = "";
			//		for(int i=1;i<=nCols_decoy;i++)
			//			colNames += rs.getMetaData().getColumnName(i)+"\t";
			//		System.out.println("\n\n"+colNames+"\n");
			while(rs.next()){
				toOut = "";
				for(int i=1;i<=nCols_decoy;i++){
					tmp = rs.getObject(i);
					if(tmp != null)
						toOut += tmp.toString()+"\t";
					else
						toOut += "null\t";
				}
				out.write(toOut+"\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\t\n");
			}
			out.flush();
			System.out.println("Done");


			rs.close();
			out.close();
		}
		st.close();
		db.closeConnection();
		System.out.println("All Done!");
	}



	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {

		if(args.length == 0){

			// spectra IDs from X!Tandem:
			//String idFormat = MergeRNAseqWithProteinIDs.OPT_FORMAT_CUFFMERGE;
			//String spectraPath = "/Users/robk/Desktop/MergeProteomicsRNAseq/DetectedTranscripts_Total_Fpkm0.1.justCuff.TandemOut";
			//spectraPath = "/Users/robk/Desktop/MergeProteomicsRNAseq/DetectedTranscripts_L10a_Fpkm0.1.TandemOut.all.withReverse";
			//spectraPath = "/Users/robk/Desktop/MergeProteomicsRNAseq/HEK_C8_5600_RNAseq_Total_Fpkm0.1.TandemOut";
			//			spectraPath = "/Users/robk/Desktop/MergeProteomicsRNAseq/HEK_C8_Orbi_CEXmerge_DetectedTranscripts_L10a_Fpkm0.1.TandemOut.all";
			//			spectraPath = "/Users/robk/Desktop/MergeProteomicsRNAseq/F456743_HeLa_DetectedTranscripts_Total_Fpkm0.1.TandemOut.all.withReverse";


			//			idFormat = MergeRNAseqWithProteinIDs.OPT_FORMAT_GENCODE;
			//			spectraPath = "/Users/robk/Desktop/MergeProteomicsRNAseq/F456743_HeLa_StockGencode14ProteinCoding.TandemOut.all.withReverse";
			//			spectraPath = "/Users/robk/Desktop/MergeProteomicsRNAseq/HEK_C8_5600_StockGencode14ProteinCoding.TandemOut.all.withReverse";


			//			idFormat = MergeRNAseqWithProteinIDs.OPT_FORMAT_SWISSPROT;
			//			spectraPath = "/Users/robk/Desktop/MergeProteomicsRNAseq/HEK.C8.5600TT.SwissprotIncIso.TandemOutput.all.xml";



			//String annotationDB = "jdbc:h2:~/Desktop/MergeProteomicsRNAseq/SEXTRACTOR_DB";
			//			String cuffmergeDB = "jdbc:h2:~/Desktop/MergeProteomicsRNAseq/ANNOTATION";
			//			args = new String[]{"-"+MergeRNAseqWithProteinIDs.OPT_URL_ANNOTATION, annotationDB, "-"+MergeRNAseqWithProteinIDs.OPT_URL_CUFFMERGE, cuffmergeDB, "-"+MergeRNAseqWithProteinIDs.OPT_PATH_SPECTRA, spectraPath};

			//MergeRNAseqWithProteinIDs.OPT_CHOICE_OVERWRITESPECTRADB
			//args = new String[]{"-"+MergeRNAseqWithProteinIDs.OPT_CHOICE_OVERWRITESPECTRADB, "-"+idFormat, "-"+MergeRNAseqWithProteinIDs.OPT_URL_ANNOTATION_ALL, annotationDB, "-"+MergeRNAseqWithProteinIDs.OPT_PATH_SPECTRA, spectraPath};
		}

//		args = new String[]{"ReadTandem","--IDgencode","-i","HEK_C8_5600_RNAseq_L10a_Fpkm0.1.TandemOut.oneSeq","-S","SampleDB_L10a.sqlite","--pathToAnnotationDB","/Users/rob/Work/YALE_OFFLINE/gencode14.sqlite"};
		new MergeRNAseqWithProteinIDs(args);

	}


}
