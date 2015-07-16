package objects;

/**
 * Very basic object to hold alignment info for a given spectra/peptide/protein
 * @author robk
 *
 */
public class Alignment{
	private Spectra spectra;
	private Peptide peptide;
	private Isoform protein;
	
	public Alignment(Spectra spec, Peptide pep, Isoform prot){
		spectra = spec;
		peptide = pep;
		protein = prot;
	}
	
	public Spectra getSpectra(){ return spectra; }
	public Peptide getPeptide(){ return peptide; }
	public Isoform getProtein(){ return protein; }
}