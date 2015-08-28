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
	private String str1;

	/**
	 * The second input String
	 */
	private String str2;

	/**
	 * The lengths of the input strings
	 */
	private int length1, length2;

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
	static final int INDEL_SCORE = -10;

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

	
	public Aligner_SmithWaterman(String str1, String str2) {
		this.str1 = reverse(str1); 
		this.str2 = reverse(str2);
		//this.str1 = str1;
		//this.str2 = str2;
		length1 = str1.length();
		length2 = str2.length();

		score = new double[length1+1][length2+1];
		prevCells = new int[length1+1][length2+1];

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

		return (str1.charAt(i - 1) == str2.charAt(j  - 1)) ? MATCH_SCORE : MISMATCH_SCORE;
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
		for (i = 1; i <= length1; i++) {
			score[i][0] = 0;
			prevCells[i][0] = DR_ZERO;
		}

		// the first column
		for (j = 1; j <= length2; j++) {
			score[0][j] = 0;
			prevCells[0][j] = DR_ZERO;
		}

		// the rest of the matrix
		for (i = 1; i <= length1; i++) {
			for (j = 1; j <= length2; j++) {
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
		for (int i = 1; i <= length1; i++) {
			for (int j = 1; j <= length2; j++) {
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
			printAlignments(i-1, j, str1.charAt(i-1) + aligned1, "_" + aligned2, print);
		}
		if ((prevCells[i][j] & DR_UP) > 0) {
			printAlignments(i, j-1, "_" + aligned1, str2.charAt(j-1) + aligned2, print);
		}
		if ((prevCells[i][j] & DR_DIAG) > 0) {
			printAlignments(i-1, j-1, str1.charAt(i-1) + aligned1, str2.charAt(j-1) + aligned2, print);
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


	/**
	 * Find the local alignments with the maximum score.
	 */
	public void findBestAlignments(boolean print){
		double maxScore = getMaxScore();
		// skip the first row and column
		for (int i = 1; i <= length1; i++) {
			for (int j = 1; j <= length2; j++) {
				if (score[i][j] == maxScore) {
					//System.out.println("i="+i+"  j="+j);
					//System.out.println("length1-i="+(length1-i)+"  length2-j="+(length2-j));
					alignmentStart_query = length1-i+1;
					alignmentStart_reference = length2-j+1;
					printAlignments(i, j, "", "", print);
				}
			}
		}
	}


	/**
	 * print the dynamic programming matrix
	 */
	public void printDPMatrix()
	{
		System.out.print("   ");
		for (int j=1; j<=length2;j++)
			System.out.print ("   "+str2.charAt(j-1));
		System.out.println();
		for (int i=0; i<=length1; i++)
		{
			if (i>0)
				System.out.print(str1.charAt(i-1)+" ");
			else 
				System.out.print("  ");
			for (int j=0; j<=length2; j++)
			{
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
		for (int i = 1; i <= length1; i++) {
			for (int j = 1; j <= length2; j++) {
				if (score[i][j] > scoreThreshold && score[i][j]>score[i-1][j-1]
						&& score[i][j]>score[i-1][j] && score[i][j]>score[i][j-1])
				{
					if (i==length1 || j==length2 ||  score[i][j]>score[i+1][j+1])
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

	
	public void printAlignmentInfo(){
		System.out.println("The maximum alignment score is: " +getAlignmentScore());
		System.out.println("read.length() = "+str1.length());
		System.out.println("adapter1.length() = "+str2.length());
		System.out.println("N matched bases = "+getNumberOfMatches() +" -- "+(getNumberOfMatches()/(str1.length()+str2.length()-getNumberOfMatches()+0.0)));
		System.out.println("Alignment start on query: "+getAlignmentStart_query());
		System.out.println("read:\n"+reverse(str1));
		for(int i=0;i<getAlignmentStart_query()-getAlignmentStart_reference();i++){ System.out.print(" "); }
		System.out.println(reverse(str2));
	}
	
	public static void main(String[] args) throws IOException{

		//String read =           "GGGTGCCAATGAACTCCAGTCACCGATGT";  //29 long
		String read = "GNCTGGTCCGATGGTAGTGGGTTATCAGAACTTGGAATTCTCGGGTGCCA"; //50 long
		String adapter1 =                             "TGGAATTCTCGGGTGCCAAGG";

		//String str1 =            "GGGTGCCAAG GAACTCCAGTCACCGATGT";
		String adapter2= "AGATCGGAAGAGCACACGTCTGAACTCCAGTCAC";

		String adapter3 = "TTTTGGGTGCCAAGGAACTCCAGTCACCGATGTAAAAAA";

		Aligner_SmithWaterman sw1 = new Aligner_SmithWaterman(read, adapter1);
		Aligner_SmithWaterman sw2 = new Aligner_SmithWaterman(read, adapter2);
		Aligner_SmithWaterman sw3 = new Aligner_SmithWaterman(read, adapter3);

		boolean printBestAlignmentString = true;

		sw1.findBestAlignments(printBestAlignmentString);
		sw1.printAlignmentInfo();
		//sw1.printDPMatrix();
		
		System.out.println("\n\n");

		sw2.findBestAlignments(printBestAlignmentString);
		sw2.printAlignmentInfo();
		//sw2.printDPMatrix();
		
		System.out.println("\n\n");

		sw3.findBestAlignments(printBestAlignmentString);
		sw3.printAlignmentInfo();
		//sw3.printDPMatrix();
		
		
		
	}
}
