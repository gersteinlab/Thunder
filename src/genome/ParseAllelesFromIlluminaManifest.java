package genome;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import annotation.Chromosome;

public class ParseAllelesFromIlluminaManifest {

	public ParseAllelesFromIlluminaManifest(){};
	
	/**
	 * Convert Illumina probe info to actual A/B alleles
	 * 
	 * @param ilmnStrand
	 * @param snp
	 * @return
	 */
	public static String[] parseGenotype(String ilmnStrand, String snp, String genomicStrand){
		String a1 = snp.substring(1, 2);
		String a2 = snp.substring(3, 4);
		
		//System.out.println(a1+"-"+a2);
		
		String alleleA = "";
		String alleleB = "";
		
		if(a1.equals("A") || a1.equals("T")){  // if A or T is the first base
			if(a2.equals("C")){
				alleleA = a1;
				alleleB = "C";
			}else if(a2.equals("G")){
				alleleA = a1;
				alleleB = "G";
			}else{
				// must be A/T or T/A
				if(ilmnStrand.equals("TOP")){
					alleleA = "A";
					alleleB = "T";
				}else{
					alleleA = "T";
					alleleB = "A";
				}
			}
		}else if(a1.equals("C") || a1.equals("G")){  // if C or G is the first base
			if(a2.equals("A")){
				alleleA = "A";
				alleleB = a1;
			}else if(a2.equals("T")){
				alleleA = "T";
				alleleB = a1;
			}else{
				// must be C/G or G/C
				if(ilmnStrand.equals("TOP")){
					alleleA = "C";
					alleleB = "G";
				}else{
					alleleA = "G";
					alleleB = "C";
				}
			}
		}
		
		if(genomicStrand.equals("-")){
			if(alleleA.equals("A"))
				alleleA = "T";
			else if(alleleA.equals("C"))
				alleleA = "G";
			else if(alleleA.equals("G"))
				alleleA = "C";
			else if(alleleA.equals("T"))
				alleleA = "A";
			
			if(alleleB.equals("A"))
				alleleB = "T";
			else if(alleleB.equals("C"))
				alleleB = "G";
			else if(alleleB.equals("G"))
				alleleB = "C";
			else if(alleleB.equals("T"))
				alleleB = "A";
		}
		
		return new String[]{alleleA,alleleB};
	}
	
	
	/**
	 * Read the Illumina array manifest
	 * 
	 * @param input
	 * @throws IOException
	 */
	public static void parseManifest(String input) throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(input));
		System.out.println("Name\tIlluminaProbeID\tAllele_A\tAllele_B\tGenomeBuild\tChromosome\tPosition");
		
		String line = "";
		String[] alleles;
		String[] bits;
		boolean header = true;
		while((line=in.readLine())!=null){
			if(header){
				if(line.startsWith("IlmnID"))
					header = false;
			}else{
				bits = line.split(",");
				if(bits.length == 21){
					alleles= parseGenotype(bits[2], bits[3], bits[20]);
					System.out.println(bits[1]+"\t"+bits[4]+"\t"+alleles[0]+"\t"+alleles[1]+"\t"+bits[8]+"\t"+bits[9]+"\t"+bits[10]);
				}
			}
		}
		in.close();
	}
	
	
	/**
	 * Read the Illumina array manifest and populate the alleles with the dbSNP variant info (i.e. determine which allele is ref and which is alt)
	 * 
	 * @param input
	 * @param variants
	 * @throws IOException
	 */
	public static HashMap<Integer, HashMap<Integer, String[]>> readManifest(String input) throws IOException{
		HashMap<Integer, HashMap<Integer, String[]>> manifestProbes = new HashMap<Integer, HashMap<Integer, String[]>>(); 
		
		System.err.print("Reading Illumina manifest...");
		BufferedReader in = new BufferedReader(new FileReader(input));
		
		int tmp_chrom;
		String line = "";
		String[] alleles;
		String[] bits;
		boolean header = true;
		while((line=in.readLine())!=null){
			if(header){
				if(line.startsWith("IlmnID"))
					header = false;
			}else{
				bits = line.split(",");
				if(bits.length == 21){
					alleles = parseGenotype(bits[2], bits[3], bits[20]);
					
					tmp_chrom = Chromosome.getChromIndex(bits[9]);
					if(!manifestProbes.containsKey(tmp_chrom)){
						manifestProbes.put(tmp_chrom, new HashMap<Integer, String[]>());
					}
					manifestProbes.get(tmp_chrom).put(Integer.valueOf(bits[10]), new String[]{bits[1],bits[4],alleles[0],alleles[1],bits[8],bits[9],bits[10],bits[20]});
					
				}
			}
		}
		in.close();
		System.err.println("Done.");
		return(manifestProbes);
	}
	
	
	/**
	 * Read the dbSNP VCF
	 * 
	 * @param args
	 * @throws IOException
	 */
	public static void parseManifest(HashMap<Integer, HashMap<Integer, String[]>> manifestProbes, String vcfFile, String sampleID) throws IOException{
		
		System.err.print("Reading variants...");
		VCFparser parser = new VCFparser(vcfFile,sampleID);
		
		System.out.println("Name\tIlluminaProbeID\tAllele_A\tAllele_B\tGenomeBuild\tChromosome\tPosition\tStrand\tdbSNP_ID\tRefAllele\tAltAllele\tAlleles_match_dbSNP");
		
		String[] tmp;
		Variant var;
		boolean refAndAltMatchAlleles;
		while (parser.hasMoreInput()) {
			var = parser.parseLine();
			
			if (var == null) continue;
			//if (var.getMaternal() == 0 && var.getPaternal() == 0) continue;
			
			if(var.getReference().length() == 1  &&  var.getAlternative().length() == 1){
				
				var.getChromosome();
				
				if(manifestProbes.containsKey(var.getChromosome())){
					if(manifestProbes.get(var.getChromosome()).containsKey(var.getPosition())){
						tmp = manifestProbes.get(var.getChromosome()).get(var.getPosition());
						refAndAltMatchAlleles = allelesMatchdbSNP(tmp[2],tmp[3],var.getReference(),var.getAlternative());
						System.out.println(tmp[0]+"\t"+tmp[1]+"\t"+tmp[2]+"\t"+tmp[3]+"\t"+tmp[4]+"\t"+tmp[5]+"\t"+tmp[6]+"\t"+tmp[7]+"\t"+var.getID()+"\t"+var.getReference()+"\t"+var.getAlternative()+"\t"+refAndAltMatchAlleles);
					}
				}else{
					// error?
				}
			}
		}
		System.err.println("Done.");
	}
	
	
	private static boolean allelesMatchdbSNP(String alleleA, String alleleB, String ref, String alt){
		boolean result = false;
		if(alleleA.equals(ref) || alleleA.equals(alt))
			if(alleleB.equals(ref) || alleleB.equals(alt))
				result = true;
		return(result);
	}
	
	
	public static void main(String[] args) throws IOException {
		//args = new String[]{"/Users/robk/WORK/YALE_offline/BrainSeq/ALLELE_SPECIFIC_EXPRESSION/SNP_files/Chip/humanomni2.5-4v1_h.csv"};
		//args = new String[]{"/Users/robk/WORK/YALE_offline/BrainSeq/ALLELE_SPECIFIC_EXPRESSION/SNP_files/Chip/humanomni2.5-4v1_h.csv","/Users/robk/WORK/YALE_offline/ANNOTATIONS/dbSNP/00-All.vcf"};
		
		if(args.length == 0)
			System.err.println("Usage: java -Xmx2G ParseAllelesFromIlluminaManifest <IlluminaArrayManifest.csv> <dbSNP_annotation.vcf>");
		else if(args.length == 1)
			ParseAllelesFromIlluminaManifest.parseManifest(args[0]);
		else if(args.length == 2)
			ParseAllelesFromIlluminaManifest.parseManifest(readManifest(args[0]), args[1], "SAMPLE");
		
	}

}
