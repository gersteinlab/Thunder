package transcriptome;

import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.Writer;
import java.util.HashMap;

import main.Thunder;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;

public class ModifyHeadersInTophatTranscriptome {

	public ModifyHeadersInTophatTranscriptome(String annotationPath) throws IOException{
		// Read transcript/gene map from annotation
		System.out.print("Reading transcript-gene mapping info from annotation file: "+annotationPath);
		readGTFAnnotation(annotationPath);
		System.out.println("Done");
	}

	
	/**
	 * Read and modify tophat headers
	 * @param tophatTranscriptomePath
	 * @throws IOException 
	 */
	public void modifyFasta(String tophatTranscriptomePath) throws IOException{
		BufferedReader in = new BufferedReader(new FileReader(tophatTranscriptomePath));
		Writer out = new BufferedWriter(new FileWriter(tophatTranscriptomePath+".modified"));
		
		String tmp, transcriptID;
		
		while((tmp=in.readLine()) != null){
			if(tmp.startsWith(">") | tmp.startsWith("@")){
				transcriptID = tmp.split(" ")[1];
				out.write(">"+transcriptID+" GN="+transcript2Gene.get(transcriptID)+"\n");
			}else{
				out.write(tmp+"\n");
			}
		}
		in.close();
		out.flush();
		out.close();
	}

	
	
	private HashMap<String, String> transcript2Gene = new HashMap<String, String>();
	private void readGTFAnnotation(String path) throws IOException{
		BufferedReader buffer = new BufferedReader(new FileReader(path));
		String line = "";
		String transcriptID, geneID;
		String[] bits;
		while((line=buffer.readLine()) != null){
			if(!line.startsWith("#")){
				bits = line.split(" |\t");
				if(bits.length >= 12){

					transcriptID = bits[11].substring(1, bits[11].length()-2);
					geneID = bits[9].substring(1, bits[9].length()-2);

					if(!transcript2Gene.containsKey(transcriptID)){
						transcript2Gene.put(transcriptID, geneID);
					}
				}
			}
		}
		buffer.close();
	}


	/**
	 * Specify the available command-line interface options
	 * @return
	 */
	//@SuppressWarnings("static-access")
	public static Options getCmdLineOptions(){
		Options options = new Options();
		//		options.addOption(OptionBuilder.withArgName("outputPath").hasArg().withDescription("output sequence lengths to a file [if not specified, lengths are printed to stdout]").create(Thunder.OPT_PATH_OUTPUT));
		//		options.addOption(OptionBuilder.withArgName("maxExpectedLength").hasArg().withDescription("print sequences exceeding this maximum expected size (for debugging + QC)").create("m"));
		return options;
	}





	public static void main(String[] args) throws Exception{

//		args = new String[]{"ModifyTophatHeaders", "/Users/robk/WORK/YALE_offline/ANNOTATIONS/gencode.v18.annotation.gtf", "NULL"};

		CommandLine cmdArgs = Thunder.parseArgs(args, getCmdLineOptions());
		int nArgs = cmdArgs.getArgList().size();
		
		if(nArgs >= 3){
			ModifyHeadersInTophatTranscriptome engine = new ModifyHeadersInTophatTranscriptome(cmdArgs.getArgList().get(1).toString());
			
			String tmpFile = "";
			for(int i=2;i<nArgs;i++){
				tmpFile = cmdArgs.getArgList().get(i).toString();
				System.out.print("Modifying file: "+tmpFile+"...");
				engine.modifyFasta(tmpFile);
				System.out.println("Done");
			}

		}else{
			HelpFormatter formatter = new HelpFormatter();
			formatter.printHelp(200, Thunder.THUNDER_EXE_COMMAND+" ModifyTophatHeaders <gtf annotation> <tophat transcriptome fasta> [additional tophat transcriptome fasta files...]", "", getCmdLineOptions(), "");
			System.out.println();
		}
	}

}
