package transcriptome;

import genome.Variant;

import java.io.File;
import java.sql.ResultSet;
import java.sql.SQLException;
import java.sql.Statement;
import java.util.Iterator;
import java.util.Map.Entry;
import java.util.TreeMap;

import net.sf.samtools.Cigar;
import net.sf.samtools.CigarElement;
import net.sf.samtools.CigarOperator;
import net.sf.samtools.SAMFileReader;
import net.sf.samtools.SAMRecord;
import net.sf.samtools.SAMRecordIterator;
import annotation.Chromosome;
import database.DBConnect;
import database.DBConnect_SQLite;

public class ParseAlignmentsCoveringVariants {

	private static boolean verbose = true;
	private Entry<Integer,Variant> currentVariant;
	private Integer currentVariantPosition;
	private TreeMap<Integer, Variant> variants;

	public void readSAM(File input, String dbPath) throws Exception{
		String tableName = "genomicVariants";
		DBConnect db = new DBConnect_SQLite(dbPath);
		Statement st = db.createStatement();

		if(db.containsTable(tableName)){

			SAMFileReader inputSam = new SAMFileReader(input);

			String lastChromosome = "";

			SAMRecord thisRecord;
			SAMRecordIterator it = inputSam.iterator();
			
			Iterator<Integer> it_var;
			int thisVarPos;
			while(it.hasNext()){
				thisRecord = it.next();
				if(thisRecord.getIntegerAttribute("NH").equals(1)){
					if(verbose) System.out.print(thisRecord.getReferenceName()+":"+thisRecord.getAlignmentStart()+"-"+thisRecord.getAlignmentEnd()+"   "+thisRecord.getReadString());
					if(verbose) System.out.println("   NM:"+thisRecord.getIntegerAttribute("NM")+"   NH:"+thisRecord.getIntegerAttribute("NH")+"   "+thisRecord.getCigarString());

					if( ! thisRecord.getReferenceName().equals(lastChromosome) ){ //new chromosome
						//						variants = resultSet2HashMap(st.executeQuery("SELECT * FROM "+tableName+" WHERE CHROMOSOME = '"+thisRecord.getReferenceName()+"' ORDER BY position ASC LIMIT 5"));
						variants = resultSet2HashMap(st.executeQuery("SELECT * FROM "+tableName+" WHERE CHROMOSOME = '"+thisRecord.getReferenceName()+"' AND position > 35400000 AND position < 35500000 AND genotype = '1/1' ORDER BY position ASC LIMIT 10"));
						lastChromosome = thisRecord.getReferenceName(); 
						thisVarPos = 0;
					}


					/*currentVariantPosition = 14501;
					//					currentVariantPosition = 14502;
					currentVariantPosition = 15005;
					currentVariantPosition = 15078;
					currentVariantPosition = 15835;*/

					it_var = variants.keySet().iterator();
					while(it_var.hasNext()){
						thisVarPos = it_var.next();
						
						if(thisVarPos < thisRecord.getAlignmentStart()){
							it_var.remove();
							variants.remove(thisVarPos);
						}else if(thisVarPos <= thisRecord.getAlignmentEnd()){
							// variant within read span
							int readPos = readVariantPosition(thisRecord.getAlignmentStart(), thisRecord.getAlignmentEnd(), thisRecord.getCigar());
							if(readPos >= 0){ // no variant within this read, remove variant
								System.out.println(" -- readCoversSNP: "+readPos+" base = "+thisRecord.getReadString().charAt(readPos-1));
							}
						}
						
					}
					


				}
			}

			inputSam.close();

		}else{
			//quit
		}

		st.close();
		db.closeConnection();
	}



	/**
	 * If the read maps over the variant, return the position in the read sequence that corresponds to the variant base
	 * @param readStart
	 * @param readEnd
	 * @param cigar
	 * @return
	 */
	public int readVariantPosition(int readStart, int readEnd, Cigar cigar){
		int ret = -1;
		CigarElement el;
		CigarOperator op;

		if(currentVariantPosition >= readStart  &&  currentVariantPosition <= readEnd){
			int var_offset = currentVariantPosition - readStart;
			int tmp = 0;
			int cigarSum_all = 0;
			int cigarSum_ref = 0;
			int cigarSum_read = 0;


			if(verbose) System.out.println("SNP at: "+currentVariantPosition+" (+"+var_offset+" from read-start)");

			for(int i=0;i<cigar.numCigarElements();i++){
				//get the cigar operator for this element
				el = cigar.getCigarElement(i);
				op = el.getOperator();
				tmp = el.getLength();


				if(verbose)System.out.print("cigar: "+op.name()+" "+tmp);

				if(op.consumesReferenceBases()  &&  op.consumesReadBases()){
					cigarSum_all += tmp;
					cigarSum_ref += tmp;
					cigarSum_read += tmp;
					if(verbose) System.out.print(" [ref="+tmp+"] [read="+tmp+"]");
				}else if(op.consumesReferenceBases()){ // if the cigar operator refers to reference bases 
					cigarSum_ref += tmp;
					cigarSum_all += tmp;
					if(verbose) System.out.print(" [ref="+tmp+"]");
				}else if(op.consumesReadBases()){ // if the cigar operator refers to read bases 
					cigarSum_read += tmp;
					if(verbose) System.out.print(" [read="+tmp+"]");
				}

				if(verbose) System.out.print(" [sum_ref="+cigarSum_ref+"]");
				if(verbose) System.out.print(" [sum_read="+cigarSum_read+"]");
				if(verbose) System.out.print(" [sum_all="+cigarSum_all+"]");

				if(cigarSum_ref > var_offset){ // this cigar element contains the variant
					if(op.consumesReadBases()){ // if this cigar element is a gap in the read, ignore it!
						ret = cigarSum_read - (cigarSum_all - var_offset) + 1;
						if(verbose) System.out.println(" [ret = "+ret+"]");
					}
					break;
				}
				if(verbose) System.out.println();
			}
		}
		return ret;
	}



	/**
	 * Load variants from the DB returned in the given ResultSet object 
	 * @param rs
	 * @return
	 * @throws SQLException
	 */
	public TreeMap<Integer, Variant> resultSet2HashMap(ResultSet rs) throws SQLException{
		TreeMap<Integer, Variant> hm = new TreeMap<Integer, Variant>();
		while(rs.next())
			hm.put(rs.getInt(2), new Variant(Chromosome.getChromIndex(rs.getString(1)),rs.getInt(2),rs.getString(4).toUpperCase(),new String[]{rs.getString(5).toUpperCase()},rs.getString(7)));
		currentVariant = hm.firstEntry();
		currentVariantPosition = hm.firstKey(); 
		return hm;
	}



	/**
	 * @param args
	 * @throws Exception 
	 */
	public static void main(String[] args) throws Exception {
		ParseAlignmentsCoveringVariants engine = new ParseAlignmentsCoveringVariants();
		engine.readSAM(new File("/Users/robk/Desktop/tmp.sam"), "/Users/robk/Desktop/SAMPLEDB_HSB123.sqlitedb");
	}

}
