package genome;

import java.io.BufferedReader;
import java.io.FileReader;
import java.util.HashMap;

import annotation.Chromosome;

public class VCFparser {
	private FileReader     _fr   = null;
	private BufferedReader _br   = null;
	private String         _line = "";
	private int _format_ind = -1, _id_ind = -1, _info_ind = -1;
	private String _sampleID = "";

	public VCFparser(String fileName,String id)
	{
		try {
			_fr = new FileReader(fileName);
			_br = new BufferedReader(_fr);
			readLine();
		} catch (Exception ex) {
			System.err.println("Can't open file " + fileName);
			System.err.println(ex.toString());
		}
		if (id != null) _sampleID = id;
	}

	public boolean hasMoreInput()
	{
		return (_line != null);
	}

	public Variant parseLine()
	{
		String line = _line;
		readLine();

		if (line == null || line.length() == 0) return null;

		//		StringTokenizer toks = new StringTokenizer(line);
		String[] bits = line.split("\t");
		if (line.startsWith("#")) {
			if (line.startsWith("#CHROM")) {
				for(int i=0;i<bits.length;i++){
					if (bits[i].equals("INFO"))    _info_ind   = i;
					if (bits[i].equals("FORMAT"))    _format_ind   = i;
					if (bits[i].equals(_sampleID))    _id_ind   = i;
				}
			}
			return null;
		}


		int chr = -1,pos = -1;
		String REF = "", ALT = "", phase = "0/0", id = "";//,info = "";
		// Parsing chromosome
		chr = Chromosome.getChromIndex(bits[0]);
		// Parsing position
		pos = Integer.parseInt(bits[1]);
		// Parsing id
		id = bits[2];
		// Parsing reference allele
		REF = bits[3];
		// Parsing alternative allele
		ALT = bits[4];
		// Info field
		//info = bits[_info_ind];
		// Output format
		
		if(_format_ind == -1){
			
		}else{
			setAttributes(bits[_format_ind], bits[_id_ind]);
			if (atts.containsKey("GT")){
				phase = atts.get("GT");
			}
		}

		if (ALT.equals("<DEL>")) return null; // Imprecise SV

		// Upper casing
		REF = REF.toUpperCase();
		ALT = ALT.toUpperCase();    

		// Splitting
		String[] alts = ALT.split(",");
		int n = alts.length;

		// Check 
		for (int i = 0;i < n;i++) 
			if      (REF.length() == 1 && alts[i].length() == 1) ; // SNP
			else if (REF.length() == 0 || alts[i].length() == 0) {
				System.err.println("Skipping invalid record:");
				System.err.println(" " + line);
				return null;
			}

		// Adjustment of first base
		if (REF.length() > 0) {
			boolean same = true;
			for (int i = 0;i < n;i++) 
				if (alts[i].length() == 0 || REF.charAt(0) != alts[i].charAt(0)) {
					same = false;
					break;
				}
			if (same) {
				pos++;
				REF = REF.substring(1);
				for (int i = 0;i < n;i++) alts[i] = alts[i].substring(1);
			}
		}

		// Adjustment of last
		if (REF.length() > 0) {
			boolean same = true;
			int indREF = REF.length() - 1;
			for (int i = 0;i < n;i++) {
				int len = alts[i].length();
				if (len == 0 || 
						REF.charAt(indREF) != alts[i].charAt(len - 1)) {
					same = false;
					break;
				}
			}
			if (same) {
				REF = REF.substring(0,indREF);
				for (int i = 0;i < n;i++)
					alts[i] = alts[i].substring(0,alts[i].length() - 1);
			}

			// System.out.print(REF);
			// for (int i = 0;i < n;i++)
			//     System.out.print(" " + alts[i]);
			// System.out.println(pos);
		}

		return new Variant(chr,pos,REF,alts,phase,atts,id);
	}

	private void readLine()
	{
		_line = null;
		try {
			_line = _br.readLine();
		} catch (Exception ex) { }

		if (_line == null) {
			try {
				_br.close();
				_fr.close();
			} catch (Exception ex) { }
			_br = null;
			_fr = null;
		}
	}


	private HashMap<String, String> atts = new HashMap<String, String>();
	private void setAttributes(String format, String values){
		String[] keys = format.split(":");
		String[] vals = values.split(":");

		if(keys.length == vals.length){
			for(int i=0;i<keys.length;i++){
				atts.put(keys[i], vals[i]);
			}
		}
	}


}
