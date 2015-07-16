package annotation;

public class Chromosome {
	
	public static int getChromIndex(String chr){
		int ret = -1;

		if (chr.length() > 3 && chr.substring(0,3).equalsIgnoreCase("chr"))
			chr = chr.substring(3);
		if (chr.equalsIgnoreCase("X"))
			ret = 23;
		else if (chr.equalsIgnoreCase("Y"))
			ret = 24;
		else if (chr.equalsIgnoreCase("M") || chr.equalsIgnoreCase("MT")) 
			ret = 25;
		else
			try {
				ret = Integer.parseInt(chr);
			} catch (Exception e) {
				//System.err.println("Unknown chromosome " + chr + ".");
			}

		return ret;
	}
	
	public static String getChromID(int chrIndex){
		return getChromID(chrIndex, true);
	}
	
	public static String getChromID(int chrIndex, boolean addCHR){
		String ret = "";
		
		if (chrIndex==23)
			ret = "X";
		else if (chrIndex==24)
			ret = "Y";
		else if (chrIndex==25)
			ret = "M";
		else
			ret = chrIndex+"";
		
		if(addCHR)
			ret = "chr"+ret;
		
		return ret;
	}
	
}
