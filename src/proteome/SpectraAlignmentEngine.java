package proteome;

import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.ArrayList;
import java.util.Collections;
import java.util.HashMap;
import java.util.Iterator;

import objects.Alignment;
import objects.Isoform;
import objects.Peptide;
import objects.Spectra;
import utils.IO_utils;

public class SpectraAlignmentEngine {

	//private HashMap<String, Spectra> spectra = new HashMap<String, Spectra>();
	private HashMap<String, ArrayList<Alignment>> _alignments = new HashMap<String, ArrayList<Alignment>>();
	

	public HashMap<String, ArrayList<Alignment>> getAlignments(){ return _alignments; }
	
	public void mergeAlignments(HashMap<String, ArrayList<Alignment>> newAlignments){
		_alignments.putAll(newAlignments);
	}

	
	/**
	 * 
	 */
	public SpectraAlignmentEngine(){}
	public SpectraAlignmentEngine(String spectraIDprefix){
		_spectraIDprefix = spectraIDprefix;
	}
	private String _spectraIDprefix = "";

	private int alignmentCount = 0;
	public int getNumberOfAlignments(){ return alignmentCount; }


	/**
	 * 
	 * @param spec
	 * @param pep
	 * @param prot
	 */
	public void addSpectraAlignment(Spectra spec, Peptide pep, Isoform prot){
		//if(! spectra.containsKey(spec.getAttribute(Spectra.ATTRIBUTE_ID)))
		//	spectra.put(spec.getAttribute(Spectra.ATTRIBUTE_ID), spec);

		if(! _alignments.containsKey(_spectraIDprefix+spec.getAttribute(Spectra.ATTRIBUTE_ID)))
			_alignments.put(_spectraIDprefix+spec.getAttribute(Spectra.ATTRIBUTE_ID), new ArrayList<Alignment>());
		_alignments.get(_spectraIDprefix+spec.getAttribute(Spectra.ATTRIBUTE_ID)).add(new Alignment(spec, pep, prot));
		alignmentCount++;

	}


	/**
	 * 
	 * @param idString
	 * @return String[]{IsoformID, Frame, GeneID}
	 */
	public static String[] parseIsoformID(String idString){
		//System.out.println(idString);
		String[] isoformFrameGene = new String[]{"","",""};
		String[] tmpParsed;
		tmpParsed = idString.split(" ");

		for(int i=1;i<tmpParsed.length;i++){
			if(tmpParsed[i].startsWith("GN=")){
				isoformFrameGene[2] = tmpParsed[i].substring(3);
				break;
			}
		}

		//System.out.println(tmpParsed[0]);
		tmpParsed = tmpParsed[0].split("_");
		isoformFrameGene[0] = tmpParsed[0];
		isoformFrameGene[1] = tmpParsed[1];

		return isoformFrameGene;
	}



	/**
	 * 
	 */
	public void processAlignments(double fdrMax, File outputFile, boolean removeCRAPome) throws IOException{

		// TODO: really only just FDR filtering at the moment- maybe include MS1 intensities?
		filterByReverseHits(fdrMax);
		IO_utils.printLineErr("  N valid spectra:\t"+_alignments.size());

		if(removeCRAPome)
			removeCRAPomeEntries("sp|", true);

		writeAlignments(outputFile);
	}


	/*
	 * Remove spectra aligning to the CRAPome
	 */
	public void removeCRAPomeEntries(String searchString, boolean verbose){
		Alignment currentAlignment;
		String currentSpectra;
		Iterator<String> allSpectra = _alignments.keySet().iterator();
		Iterator<Alignment> currentSpectraAlignments;
		ArrayList<String> invalidSpectra = new ArrayList<String>();
		ArrayList<String> crapomeIDs = new ArrayList<String>();
		//String[] isoformFrameGene;
		while(allSpectra.hasNext()){
			currentSpectra = allSpectra.next();
			currentSpectraAlignments = _alignments.get(currentSpectra).iterator();
			while(currentSpectraAlignments.hasNext()){
				currentAlignment = currentSpectraAlignments.next();
				if(currentAlignment.getProtein().getIsoformID().contains(searchString)){
					invalidSpectra.add(currentSpectra);
					if(!crapomeIDs.contains(currentAlignment.getProtein().getIsoformID()))
						crapomeIDs.add(currentAlignment.getProtein().getIsoformID());
					break;
				}
			}
		}
		
		// remove invalid spectra
		Iterator<String> invalidIterator = invalidSpectra.iterator();
		while(invalidIterator.hasNext())
			_alignments.remove(invalidIterator.next());


		if(verbose){
			Iterator<String> it = crapomeIDs.iterator();
			IO_utils.printLineErr("Removed alignments to these CRAPome sequences: ");
			IO_utils.printErr("  ");
			while(it.hasNext())
				System.err.print(it.next()+" ");
			System.err.println("");
		}
	}



	/**
	 * 
	 * @param outputFile
	 * @throws IOException
	 */
	public void writeAlignments(File outputFile) throws IOException{
		Writer out = new BufferedWriter(new FileWriter(outputFile));
		out.write("PeptideSeq\tSpectraID\tIsoformID\tGeneID\tFrame\tStart\tStop\tExpectation_peptide\tExpectation_isoform\tMH\tCharge\tDelta\tHyperscore\tNextscore\ty_score\ty_ions\tb_score\tb_ions\tmissedCleavages\tPTMs\n");

		Alignment currentAlignment;
		Iterator<String> allSpectra = _alignments.keySet().iterator();
		Iterator<Alignment> currentSpectraAlignments;

		String[] isoformFrameGene;
		while(allSpectra.hasNext()){
			currentSpectraAlignments = _alignments.get(allSpectra.next()).iterator();
			while(currentSpectraAlignments.hasNext()){
				currentAlignment = currentSpectraAlignments.next();

				out.write(currentAlignment.getPeptide().getSequence()+"\t");
				out.write(currentAlignment.getSpectra().getAttribute(Spectra.ATTRIBUTE_ID)+"\t");

				isoformFrameGene = parseIsoformID(currentAlignment.getProtein().getIsoformID());
				out.write(isoformFrameGene[0]+"\t");
				out.write(isoformFrameGene[2]+"\t");
				out.write(isoformFrameGene[1]+"\t");

				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_START)+"\t");
				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_END)+"\t");
				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_EXPECT)+"\t");
				out.write(currentAlignment.getProtein().getAttribute(Isoform.ATTRIBUTE_EXPECT)+"\t");
				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_MH)+"\t");
				out.write(currentAlignment.getSpectra().getAttribute(Spectra.ATTRIBUTE_Z)+"\t");
				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_DELTA)+"\t");
				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_HYPERSCORE)+"\t");
				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_NEXTSCORE)+"\t");
				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_Y_SCORE)+"\t");
				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_Y_IONS)+"\t");
				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_B_SCORE)+"\t");
				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_B_IONS)+"\t");
				out.write(currentAlignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_MISSED_CLEAVAGES)+"\t");
				out.write(currentAlignment.getPeptide().getModificationsString()+"\n");

			}
			out.flush();
		}
		out.close();
	}




	/**
	 * 
	 */
	public void filterByReverseHits(double fdrMax){
		String spectraID;
		ArrayList<Alignment> tmp_alignments;
		Alignment tmp_alignment = null;

		ArrayList<Double> expectations_legit = new ArrayList<Double>();
		ArrayList<Double> expectations_decoy = new ArrayList<Double>();

		//
		// Loop through all spectra and assign expectation values to either a legit DB hit, or to the reverse decoy
		//
		Iterator<String> spectraIterator = _alignments.keySet().iterator();
		while(spectraIterator.hasNext()){
			spectraID = spectraIterator.next();

			tmp_alignments = _alignments.get(spectraID);
			Iterator<Alignment> alignmentIterator = tmp_alignments.iterator();

			int countReverseHits = 0;
			while(alignmentIterator.hasNext()){
				tmp_alignment = alignmentIterator.next();
				//System.out.println(tmp_alignment.getProtein().getAttribute(Isoform.ATTRIBUTE_LABEL));
				if(tmp_alignment.getProtein().getAttribute(Isoform.ATTRIBUTE_LABEL).endsWith(":reversed"))
					countReverseHits ++;
			}

			if(countReverseHits < tmp_alignments.size())
				expectations_legit.add(Double.valueOf(tmp_alignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_EXPECT)).doubleValue());
			else
				expectations_decoy.add(Double.valueOf(tmp_alignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_EXPECT)).doubleValue());
		}


		//
		// Calculate maximum expectation value leading to the desired FDR 
		//
		double expectationMax = calculateMaxExpectation(expectations_legit, expectations_decoy, fdrMax);
		//System.out.println("expectationMax = "+expectationMax);


		//
		// Loop through all spectra and remove those that are not valid hits
		//
		ArrayList<String> invalidSpectraIDs = new ArrayList<String>();
		spectraIterator = _alignments.keySet().iterator();
		while(spectraIterator.hasNext()){
			spectraID = spectraIterator.next();
			tmp_alignments = _alignments.get(spectraID);
			Iterator<Alignment> alignmentIterator = tmp_alignments.iterator();
			ArrayList<Alignment> invalidAlignments = new ArrayList<Alignment>();
			while(alignmentIterator.hasNext()){
				tmp_alignment = alignmentIterator.next();
				if(Double.valueOf(tmp_alignment.getPeptide().getAttribute(Peptide.ATTRIBUTE_EXPECT)).doubleValue() > expectationMax){
					// Flag spectra for removal if expectation value is too high
					invalidSpectraIDs.add(spectraID);
					break;
				}else if(tmp_alignment.getProtein().getAttribute(Isoform.ATTRIBUTE_LABEL).endsWith(":reversed")){
					// Remove alignment if it is to the reverse DB
					invalidAlignments.add(tmp_alignment);
				}
			}

			// remove invalid alignments
			Iterator<Alignment> invalidIterator = invalidAlignments.iterator();
			while(invalidIterator.hasNext())
				_alignments.get(spectraID).remove(invalidIterator.next());

			// if there are no valid alignments for this spectra, flag it for removal
			if(_alignments.get(spectraID).size() == 0){
				invalidSpectraIDs.add(spectraID);
			}
		}
		//System.out.println("N invalid spectra = "+invalidSpectraIDs.size());


		//
		// Loop through exclude list and remove bad spectra
		//
		Iterator<String> removalIterator = invalidSpectraIDs.iterator();
		while(removalIterator.hasNext()){
			_alignments.remove(removalIterator.next());
		}
	}



	/**
	 * 
	 * @param expectations_legit
	 * @param expectations_decoy
	 * @param fdrMax
	 * @return
	 */
	public double calculateMaxExpectation(ArrayList<Double> expectations_legit, ArrayList<Double> expectations_decoy, double fdrMax){
		Collections.sort(expectations_legit);
		Collections.sort(expectations_decoy);
		Object[] legitExpects = expectations_legit.toArray();
		Object[] decoyExpects = expectations_decoy.toArray();

		int countLegit = 0;
		int countDecoy = 0;
		double[][] fdrMatrix = new double[legitExpects.length][2];

		int decoyIndex = 0;
		for(int i=0;i<legitExpects.length;i++){
			if(decoyIndex < decoyExpects.length){
				//System.out.println(legitExpects[i]+"\t"+decoyExpects[decoyIndex]);

				if((Double)legitExpects[i] <= (Double)decoyExpects[decoyIndex])
					countLegit++;
				else{
					decoyIndex++;
					countDecoy++;
				}
				fdrMatrix[i][0] = (Double)legitExpects[i];
				fdrMatrix[i][1] = (countDecoy+0.0)/(countLegit+countDecoy+0.0);

			}else{
				break;
			}
		}

		double maxExpect = 1000000.0;
		//System.out.println();
		for(int i=0;i<(countLegit);i++){
			//System.out.println(fdrMatrix[i][0]+" : "+fdrMatrix[i][1]);
			if(fdrMatrix[i][1] > fdrMax){
				maxExpect = fdrMatrix[i][0]; 
				break;
			}
		}
		
		IO_utils.printLineErr("  Max expectation value for this FDR = "+maxExpect);

		return maxExpect;
	}



	public int countSpectra(){ return this._alignments.size(); }
	//public int countPeptides(){ return this.peptides.size(); }
	//public int countProteins(){ return this.proteins.size(); }

}



