# Thunder
Bunch of tools for RNA and protein analysis

Usage:	 java -Xmx2G -jar Thunder.jar <Command>

Available Commands: 

	GetSequenceLengths		| Get the distribution of sequence lengths in a FASTA/Q file

	FastaHeaderGrep			| Filter fasta sequences based on the sequence ID

	FilterSequencesByLength		| Filter fasta or fastq sequences based on some maximum sequence length

	RemoveHomopolymerRepeat		| Filter fasta or fastq sequences based on sequence composition

	MatchPairedEndSequence		| Match paired-end fastq sequences based on readID

	GTF2Fasta			| Extract GTF coordinates from FASTA sequence(s)

	ParseTandemOutput		| Process X!Tandem MS/MS-spectra alignments and output summary table

	IsoformEM_Footprint		| Infer most likely transcripts from ribosome footprint alignments

	IsoformEM_Proteomic		| Infer most likely isoforms from MS/MS spectra mapping

	CIGAR_2_PW			| Reads SAM alignments and converts the CIGAR strings to a position-weight matrix

	ProcessEndogenousAlignments	| Process endogenous smallRNA alignments for the exceRpt pipeline

