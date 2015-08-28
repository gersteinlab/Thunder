package sequenceTools;

public class Sequence {

	
	/**
	 * Calculates the reverse complement of the given sequence
	 * 
	 * @param seq
	 * @return
	 */
	public static String reverseComplement(String seq){
		final StringBuilder out = new StringBuilder(seq.length());
		for (int i = seq.length(); i > 0; i--)
			try{
				out.append(replace(seq.charAt(i-1)));
			}catch(IllegalArgumentException e){
				out.append("N");
			}
		return out.toString();
	}

	/**
	 * Specifies the reverse complement of each nucleotide 
	 * 
	 * @param in
	 * @return
	 */
	private static char replace(char in) {
		switch (in) {
		case 'A': return 'T';
		case 'T': return 'A';
		case 'C': return 'G';
		case 'G': return 'C';
		case 'N': return 'N';
		case 'a': return 't';
		case 't': return 'a';
		case 'c': return 'g';
		case 'g': return 'c';
		case 'n': return 'n';
		default: throw new IllegalArgumentException();
		}
	}

}
