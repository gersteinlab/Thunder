package sequenceTools;


/**
 * Smith-Waterman local alignment algorithm.
 */

import java.io.IOException;
import java.util.ArrayList;
import java.util.List;

import sequenceTools.SimpleChaining.Match;

/** 
 * Design Note: this class implements AminoAcids interface: a simple fix customized to amino acids, since that is all we deal with in this class
 * Supporting both DNA and Aminoacids, will require a more general design.
 */

public class Aligner_SmithWaterman {

	private final double scoreThreshold = 19.9;

	/**
	 * The first input string
	 */
	private String _query_sequence;

	/**
	 * The second input String
	 */
	private String _reference_sequence;

	/**
	 * The lengths of the input strings
	 */
	private int _query_length, _reference_length;

	/**
	 * The score matrix.
	 * The true scores should be divided by the normalization factor.
	 */
	private double[][] score;

	/**
	 * The normalization factor.
	 * To get the true score, divide the integer score used in computation
	 * by the normalization factor.
	 */
	static final double NORM_FACTOR = 1.0;

	/**
	 * The similarity function constants.
	 * They are amplified by the normalization factor to be integers.
	 */
	//static final int MATCH_SCORE = 10;
	//static final int MISMATCH_SCORE = -8;
	//static final int INDEL_SCORE = -9;
	static final int MATCH_SCORE = 1;
	static final int MISMATCH_SCORE = -1;
	static final int INDEL_SCORE = -100;

	/**
	 * Constants of directions.
	 * Multiple directions are stored by bits.
	 * The zero direction is the starting point.
	 */
	static final int DR_LEFT = 1; // 0001
	static final int DR_UP = 2;   // 0010
	static final int DR_DIAG = 4; // 0100
	static final int DR_ZERO = 8; // 1000

	/**
	 * The directions pointing to the cells that
	 * give the maximum score at the current cell.
	 * The first index is the column index.
	 * The second index is the row index.
	 */
	private int[][] prevCells;


	public Aligner_SmithWaterman(String query_sequence, String reference_sequence) {
		_query_sequence = reverse(query_sequence); 
		_reference_sequence = reverse(reference_sequence);
		//this.str1 = str1;
		//this.str2 = str2;
		_query_length = query_sequence.length();
		_reference_length = reference_sequence.length();

		score = new double[_query_length+1][_reference_length+1];
		prevCells = new int[_query_length+1][_reference_length+1];

		buildMatrix();
	}

	/**
	 * Reverses the input String
	 * @param in
	 * @return
	 */
	public static String reverse(String in){ return new StringBuilder(in).reverse().toString(); }; 


	/**
	 * Compute the similarity score of substitution: use a substitution matrix if the cost model
	 * The position of the first character is 1.
	 * A position of 0 represents a gap.
	 * @param i Position of the character in str1
	 * @param j Position of the character in str2
	 * @return Cost of substitution of the character in str1 by the one in str2
	 */
	private double similarity(int i, int j) {
		if (i == 0 || j == 0) {
			// it's a gap (indel)
			return INDEL_SCORE;
		}

		return (_query_sequence.charAt(i - 1) == _reference_sequence.charAt(j  - 1)) ? MATCH_SCORE : MISMATCH_SCORE;
		//return Blosum.getDistance(str1.charAt(i-1), str2.charAt(j-1)); 
	}

	/**
	 * Build the score matrix using dynamic programming.
	 * Note: The indel scores must be negative. Otherwise, the
	 * part handling the first row and column has to be
	 * modified.
	 */
	@SuppressWarnings("unused")
	private void buildMatrix() {
		if (INDEL_SCORE >= 0) {
			throw new Error("Indel score must be negative");
		}

		//		if (isDistanceMatrixNull()) {
		//			throw new Error ("Distance Matrix is NULL");
		//		}

		int i; // length of prefix substring of str1
		int j; // length of prefix substring of str2

		// base case
		score[0][0] = 0;
		prevCells[0][0] = DR_ZERO; // starting point

		// the first row
		for (i = 1; i <= _query_length; i++) {
			score[i][0] = 0;
			prevCells[i][0] = DR_ZERO;
		}

		// the first column
		for (j = 1; j <= _reference_length; j++) {
			score[0][j] = 0;
			prevCells[0][j] = DR_ZERO;
		}

		// the rest of the matrix
		for (i = 1; i <= _query_length; i++) {
			for (j = 1; j <= _reference_length; j++) {
				double diagScore = score[i - 1][j - 1] + similarity(i, j);
				double upScore = score[i][j - 1] + similarity(0, j);
				double leftScore = score[i - 1][j] + similarity(i, 0);

				score[i][j] = Math.max(diagScore, Math.max(upScore,
						Math.max(leftScore, 0)));
				prevCells[i][j] = 0;

				// find the directions that give the maximum scores.
				// the bitwise OR operator is used to record multiple
				// directions.
				if (diagScore == score[i][j]) {
					prevCells[i][j] |= DR_DIAG;
				}
				if (leftScore == score[i][j]) {
					prevCells[i][j] |= DR_LEFT;
				}
				if (upScore == score[i][j]) {
					prevCells[i][j] |= DR_UP;
				}
				if (0 == score[i][j]) {
					prevCells[i][j] |= DR_ZERO;
				}
			}
		}
	}

	/**
	 * Get the maximum value in the score matrix.
	 */
	private double getMaxScore() {
		double maxScore = 0;

		// skip the first row and column
		for (int i = 1; i <= _query_length; i++) {
			for (int j = 1; j <= _reference_length; j++) {
				if (score[i][j] > maxScore) {
					maxScore = score[i][j];
				}
			}
		}

		return maxScore;
	}

	/**
	 * Get the alignment score between the two input strings.
	 */
	public double getAlignmentScore() {
		return getMaxScore() / NORM_FACTOR;
	}

	private int basesMatched = 0;
	private int basesMismatched = 0;

	/**
	 * Output the local alignments ending in the (i, j) cell.
	 * aligned1 and aligned2 are suffixes of final aligned strings
	 * found in backtracking before calling this function.
	 * Note: the strings are replicated at each recursive call.
	 * Use buffers or stacks to improve efficiency.
	 */
	private void printAlignments(int i, int j, String aligned1, String aligned2, boolean print) {
		// we've reached the starting point, so print the alignments	

		if ((prevCells[i][j] & DR_ZERO) > 0) {
			if(print){
				System.out.println(reverse(aligned1));
				System.out.println(reverse(aligned2));
				System.out.println("");
			}
			// Note: we could check other directions for longer alignments
			// with the same score. we don't do it here.
			return;
		}

		// find out which directions to backtrack
		if ((prevCells[i][j] & DR_LEFT) > 0) {
			printAlignments(i-1, j, _query_sequence.charAt(i-1) + aligned1, "_" + aligned2, print);
			basesMismatched ++;
		}
		if ((prevCells[i][j] & DR_UP) > 0) {
			printAlignments(i, j-1, "_" + aligned1, _reference_sequence.charAt(j-1) + aligned2, print);
			basesMismatched ++;
		}
		if ((prevCells[i][j] & DR_DIAG) > 0) {
			printAlignments(i-1, j-1, _query_sequence.charAt(i-1) + aligned1, _reference_sequence.charAt(j-1) + aligned2, print);
			basesMatched += 1;
		}
	}

	/** 
	 * given the bottom right corner point trace back  the top left corner.
	 *  at entry: i, j hold bottom right (end of Aligment coords)
	 *  at return:  hold top left (start of Alignment coords)
	 */
	private int [] traceback(int i, int j) {

		// find out which directions to backtrack
		while (true)
		{
			if ((prevCells[i][j] & DR_LEFT) > 0) {
				if (score[i-1][j]>0) i--;
				else	break;			    
			}
			if ((prevCells[i][j] & DR_UP) > 0) {
				//		    return traceback(i, j-1);
				if (score[i][j-1]>0) j--;
				else	break;			    
			}
			if ((prevCells[i][j] & DR_DIAG) > 0) {
				//		    return traceback(i-1, j-1);
				if (score[i-1][j-1]>0) {i--;j--;}
				else	 break;			    
			}
		}		
		int [] m ={i, j};
		return m;
	}

	private int alignmentStart_query = 0, alignmentStart_reference = 0;
	public int getAlignmentStart_query(){ return alignmentStart_query; }
	public int getAlignmentStart_reference(){ return alignmentStart_reference; }
	public int getNumberOfMatches(){ return basesMatched; }
	public int getNumberOfMismatches(){ return basesMismatched; }

	/**
	 * Find the local alignments with the maximum score.
	 */
	public void findBestAlignments(boolean print){
		double maxScore = getMaxScore();
		// skip the first row and column
		for (int i = 1; i <= _query_length; i++) {
			for (int j = 1; j <= _reference_length; j++) {
				if (score[i][j] == maxScore) {
					//System.out.println("i="+i+"  j="+j);
					//System.out.println("_query_length-i+1="+(_query_length-i+1)+"  _reference_length-j+1="+(_reference_length-j+1));
					alignmentStart_query = _query_length-i+1;
					alignmentStart_reference = _reference_length-j+1;
					printAlignments(i, j, "", "", print);
					
					i = _query_length;
					j = _reference_length;
				}
			}
		}
	}


	/**
	 * print the dynamic programming matrix
	 */
	public void printDPMatrix(){
		System.out.print("   ");
		for (int j=1; j<=_reference_length;j++)
			System.out.print ("   "+_reference_sequence.charAt(j-1));
		System.out.println();

		for (int i=0; i<=_query_length; i++){
			if (i>0)
				System.out.print(_query_sequence.charAt(i-1)+" ");
			else 
				System.out.print("  ");

			for (int j=0; j<=_reference_length; j++){
				System.out.print(score[i][j]/NORM_FACTOR+" ");
			}
			System.out.println();
		}
	}




	/**
	 *  Return a set of Matches idenfied in Dynamic programming matrix. 
	 * A match is a pair of subsequences whose score is higher than the 
	 * preset scoreThreshold
	 **/
	public List<Match> getMatches()
	{
		ArrayList<Match> matchList = new ArrayList<Match>();
		int fA=0, fB=0;
		//	skip the first row and column, find the next maxScore after prevmaxScore 
		for (int i = 1; i <= _query_length; i++) {
			for (int j = 1; j <= _reference_length; j++) {
				if (score[i][j] > scoreThreshold && score[i][j]>score[i-1][j-1]
						&& score[i][j]>score[i-1][j] && score[i][j]>score[i][j-1])
				{
					if (i==_query_length || j==_reference_length ||  score[i][j]>score[i+1][j+1])
					{
						// should be lesser than prev maxScore					    	
						fA = i; 
						fB = j;
						int [] f=traceback(fA, fB); // sets the x, y to startAlignment coordinates 
						matchList.add(new SimpleChaining.Match(f[0], i, f[1], j, score[i][j]/ NORM_FACTOR));
					}
				}
			}
		}
		return matchList; // could be empty if no HSP scores are > scoreThreshold
	}


	public double getMatchFractionOfOverlap(){
		String readSeq = reverse(_query_sequence);
		//for(int i=0;i<getAlignmentStart_query()-getAlignmentStart_reference();i++){ System.out.print(" "); }
		//System.out.println(reverse(_reference_sequence));

		int lengthNoOverlap = getAlignmentStart_query()-getAlignmentStart_reference();
		int lengthOverlap = readSeq.length()-lengthNoOverlap;

		if(readSeq.length() - (getAlignmentStart_query()-getAlignmentStart_reference()) > _reference_sequence.length()){
			lengthOverlap = _reference_sequence.length();
		}

		//System.out.println("lengthNoOverlap="+lengthNoOverlap+" lengthOverlap="+lengthOverlap+" getNumberOfMatches()/lengthOverlap="+((getNumberOfMatches()+0.0)/(lengthOverlap+0.0)));

		return((getNumberOfMatches()+0.0)/(lengthOverlap+0.0));
	}


	public double getWeightedScore(){
		double nMismatches = (getNumberOfMatches()-getAlignmentScore())/2.0;
		double newScore = (getNumberOfMatches()-nMismatches-(getAlignmentStart_reference()-1)+0.0)*Math.pow(getMatchFractionOfOverlap(), 2);
		if(newScore < 0.0)
			newScore = 0.0;
		return newScore;
	}

	public String getAlignmentInfo(){
		String out = "";

		//System.out.println("Alignment start on query: "+getAlignmentStart_query());
		//System.out.println("Alignment start on reference: "+getAlignmentStart_reference());

		if(getAlignmentStart_query() >= getAlignmentStart_reference()){
			out = reverse(_query_sequence)+"\n";
			for(int i=0;i<getAlignmentStart_query()-getAlignmentStart_reference();i++){ out+=" "; }
			out += reverse(_reference_sequence);
		}else{
			for(int i=0;i<getAlignmentStart_reference()-getAlignmentStart_query();i++){ out+=" "; }
			out += reverse(_query_sequence)+"\n";
			out += reverse(_reference_sequence);
		}

		out += "\tmaxScore: "+getMaxScore();
		out += "\tnMatch:" +getNumberOfMatches();
		out += "\tnMismatch:" +getNumberOfMismatches();
		out += "\tfracMatch:" +getMatchFractionOfOverlap();
		out += "\tscore:" +getAlignmentScore();
		out += "\tweightedScore:" +getWeightedScore();



		//System.out.println("read.length() = "+_query_sequence.length());
		//System.out.println("adapter1.length() = "+_reference_sequence.length());
		/*System.out.println("N matched bases = "+getNumberOfMatches() +" -- "+(getNumberOfMatches()/(_query_sequence.length()+_reference_sequence.length()-getNumberOfMatches()+0.0)));
		System.out.println("N matched bases = "+getNumberOfMatches() +" -- "+(getNumberOfMatches()/(_reference_sequence.length()+0.0)));
		System.out.println("N mismatched bases = "+getNumberOfMismatches() +" -- "+(getNumberOfMismatches()/(_reference_sequence.length()+0.0)));
		System.out.println("Alignment start on query: "+getAlignmentStart_query());
		System.out.println("Alignment start on reference: "+getAlignmentStart_reference());
		System.out.println("getAlignmentStart_query()-getAlignmentStart_reference() = "+(getAlignmentStart_query()-getAlignmentStart_reference()));
		 */
		//getMatchFractionOfOverlap();

		return out;
	}

	public static void main(String[] args) throws IOException{

		//String read =           "GGGTGCCAATGAACTCCAGTCACCGATGT";  //29 long
		//String read = "GNCTGGTCCGATGGTAGTGGGTTATCAGAACTTGGAATTCTCGGGTGCCA"; //50 long
		//String adapter1 =                             "TGGAATTCTCGGGTGCCAAGG";

		String read = "NNNNNNNNNNGCCCCGATCGATTCGTATGCCTTCTT";
		//String read = "CCCCCCCACACAGGACACACATACAATACAACTCACACACAACATACACA";


		String adapter_Illumina_1_0_smallRNA_3p = "TCGTATGCCGTCTTCTGCTTG";
		String adapter_Illumina_1_5_smallRNA_3p = "ATCTCGTATGCCGTCTTCTGCTTGC";
		String adapter_TruSeq_smallRNA_3p = "TGGAATTCTCGGGTGCCAAGG";
		String adapter_NEB_smallRNA_3p = "AGATCGGAAGAGCACACGTCT";
		String adapter_SOLiD_smallRNA_3p = "CGCCTTGGCCGTACAGCAG";
		String adapter_TruSeq_adapter_p7 = "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

		//String str1 =            "GGGTGCCAAG GAACTCCAGTCACCGATGT";
		//String adapter2= "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";
		//String adapter3 = "TTTTGGGTGCCAAGGAACTCCAGTCACCGATGTAAAAAA";
		//String adapter4 = "TAGTGGGTTA"; //50 long
		//String adapter5 = "BLAHBLAHGNCTGGTCCGATGG"; //50 long

		Aligner_SmithWaterman sw1 = new Aligner_SmithWaterman(read, adapter_Illumina_1_0_smallRNA_3p);
		Aligner_SmithWaterman sw2 = new Aligner_SmithWaterman(read, adapter_Illumina_1_5_smallRNA_3p);
		Aligner_SmithWaterman sw3 = new Aligner_SmithWaterman(read, adapter_TruSeq_smallRNA_3p);
		Aligner_SmithWaterman sw4 = new Aligner_SmithWaterman(read, adapter_NEB_smallRNA_3p);
		Aligner_SmithWaterman sw5 = new Aligner_SmithWaterman(read, adapter_SOLiD_smallRNA_3p);
		Aligner_SmithWaterman sw6 = new Aligner_SmithWaterman(read, adapter_TruSeq_adapter_p7);

		boolean printBestAlignmentString = true;

		sw1.findBestAlignments(printBestAlignmentString);
		System.out.println("Illumina_1_0_smallRNA_3p:");
		System.out.println(sw1.getAlignmentInfo());
		//sw1.printDPMatrix();

		System.out.println("\n\n");

		sw2.findBestAlignments(printBestAlignmentString);
		System.out.println("adapter_Illumina_1_5_smallRNA_3p:");
		System.out.println(sw2.getAlignmentInfo());
		//sw2.printDPMatrix();

		System.out.println("\n\n");

		sw3.findBestAlignments(printBestAlignmentString);
		System.out.println("adapter_TruSeq_smallRNA_3p:");
		System.out.println(sw3.getAlignmentInfo());
		//sw3.printDPMatrix();

		System.out.println("\n\n");

		sw4.findBestAlignments(printBestAlignmentString);
		System.out.println("adapter_NEB_smallRNA_3p:");
		System.out.println(sw4.getAlignmentInfo());
		//sw4.printDPMatrix();

		System.out.println("\n\n");

		sw5.findBestAlignments(printBestAlignmentString);
		System.out.println("adapter_SOLiD_smallRNA_3p:");
		System.out.println(sw5.getAlignmentInfo());
		//sw5.printDPMatrix();

		System.out.println("\n\n");

		sw6.findBestAlignments(printBestAlignmentString);
		System.out.println("adapter_TruSeq_adapter_p7:");
		System.out.println(sw6.getAlignmentInfo());
		//sw6.printDPMatrix();

	}
}
