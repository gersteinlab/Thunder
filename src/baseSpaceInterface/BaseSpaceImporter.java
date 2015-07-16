package baseSpaceInterface;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileFilter;
import java.io.FileReader;
import java.io.FilenameFilter;
import java.io.IOException;
import java.util.HashMap;
import java.util.Iterator;

import org.apache.commons.io.filefilter.DirectoryFileFilter;
import org.json.JSONArray;
import org.json.JSONObject;

public class BaseSpaceImporter {

	/**
	 * 
	 * @param genomeVersion
	 * @return
	 */
	public static String getGenomePath(String genomeVersion){
		//return hg19 by default
		String genomePath = "/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa";
		
		if(genomeVersion.equalsIgnoreCase("hg19"))
			genomePath = "/genomes/Homo_sapiens/UCSC/hg19/Sequence/WholeGenomeFasta/genome.fa";
		else if(genomeVersion.equalsIgnoreCase("mm9"))
			genomePath = "/genomes/Mus_musculus/UCSC/mm9/Sequence/WholeGenomeFasta/genome.fa";
		else if(genomeVersion.equalsIgnoreCase("mm10"))
			genomePath = "/genomes/Mus_musculus/UCSC/mm10/Sequence/WholeGenomeFasta/genome.fa";
		
		return genomePath;
	}
	
	
	/**
	 * 
	 * @return
	 * @throws IOException
	 */
	public static String readOutputProjectIDFromJSON() throws IOException{
		return(readOutputProjectIDFromJSON("/data/input/AppSession.json"));
		//return readOutputProjectIDFromJSON((new File("/data")).listFiles()[0].getAbsolutePath() + "/data/input/AppSession.json");
	}
	
	/**
	 * 
	 * @param jsonPath
	 * @return
	 * @throws IOException
	 */
	public static String readOutputProjectIDFromJSON(String jsonPath) throws IOException{
		
		BufferedReader in = new BufferedReader(new FileReader(jsonPath));
		String line = "";
		String str = new String();
		while((line=in.readLine()) != null){
			str += line;
		}
		in.close();

		
	    // build the JSON object
	    JSONObject obj = new JSONObject(str);
	    String projectID = "";
	    
	    JSONArray items = obj.getJSONObject("Properties").getJSONArray("Items");
	    JSONObject tmp;
	    for(int i=0;i<items.length();i++){
	    	tmp = (JSONObject)items.get(i);
	    	
	    	// Sets the project ID for output
	    	if(tmp.getString("Name").equals("Output.Projects")){
		    	projectID = ((JSONObject) tmp.getJSONArray("Items").get(0)).get("Id").toString();
	    	}
	    	//else if(tmp.getString("Name").equals("Input.AppResults")){
	    	//	System.out.println(tmp.getString("Name"));
		    //	System.out.println(tmp.toString());
	    	//}
	    }
	    
	    //System.out.println("\n\nProject ID = "+projectID);
		
		return projectID;
	}
	
	
	public static File getAppResultInputDirectory(){
		//File[] subdirs = (new File("/data")).listFiles((FileFilter) DirectoryFileFilter.DIRECTORY);
		//File tmp = new File(subdirs[0], "/data/input/appresults/");
		File tmp = new File("/data/input/appresults/");
		File[] tmp2 = tmp.listFiles((FileFilter) DirectoryFileFilter.DIRECTORY);
		return tmp2[0];
	}
	
	public static File getAppResultOutputDirectory(String resultID){
		//File[] subdirs = (new File("/data")).listFiles((FileFilter) DirectoryFileFilter.DIRECTORY);
		//File appResultOutputDir = new File(subdirs[0], "/data/output/appresults/");
		File appResultOutputDir = new File("/data/output/appresults/");
		
		File outputDir = new File(appResultOutputDir, "/"+resultID+"/customProteome/");
		outputDir.mkdirs();
		
		return outputDir;
	}
	
	
	/**
	 * 
	 * @return
	 */
	public static File getDefaultCuffmergeDirectory(){
		return new File(getAppResultInputDirectory(), "/differential/cuffmerge");
	}
	
	/**
	 * 
	 * @return
	 */
	public static File getDefaultCuffdiffDirectory(){
		return new File(getAppResultInputDirectory(), "/differential/cuffdiff");
	}
	
	/**
	 * 
	 * @return
	 */
	//public static String getDefaultOutputDirectory(){
	//	return (new File("/data")).listFiles()[0].getAbsolutePath() + "/data/output/appresults/";
	//}
	
	
	/**
	 * 
	 * @param basePath
	 * @return
	 */
	public static HashMap<String, File[]> searchForInputFiles(File basePath){
		return(searchForInputFiles(new String[0], basePath));
	}

	
	/**
	 * 
	 * @param fileExtensions
	 * @param basePath
	 * @return
	 */
	public static HashMap<String, File[]> searchForInputFiles(String[] fileExtensions, File basePath){
		HashMap<String, File[]> inputFiles = new HashMap<String, File[]>();
		if(basePath.isDirectory()){
			if(fileExtensions.length > 0)
				inputFiles.put(basePath.getAbsolutePath(), basePath.listFiles(new MyFilenameFilter(fileExtensions)));
			else
				inputFiles.put(basePath.getAbsolutePath(), basePath.listFiles());
		}
		return inputFiles;
	}
	
	
	
	public static void main(String[] args) throws IOException {
		
		// TEST parsing project ID from JSON
		System.out.println("Project ID = "+readOutputProjectIDFromJSON("/Users/robk/Box Sync/Work for other people/BaseSpaceApp/TEST.json"));
		//System.out.println("Project ID = "+readOutputProjectIDFromJSON());
		
		
		// TEST searching for input files		
		String[] fileTypes = new String[]{""};
		fileTypes = new String[]{".gtf"};
		fileTypes = new String[]{".diff"};
		fileTypes = new String[]{".gtf",".diff"};
		fileTypes = new String[]{".merged.gtf"};
		
		HashMap<String, File[]> mainfiles = searchForInputFiles(fileTypes, new File("/Users/robk/Box Sync/Work for other people/BaseSpaceApp/FakeInputAppResults"));
		//HashMap<String, File[]> mainfiles = searchForInputFiles(new String[0], new File("/Users/robk/Box Sync/Work for other people/BaseSpaceApp/FakeInputAppResults"));
		System.out.println("mainfiles.size() = "+mainfiles.size());
		Iterator<String> it = mainfiles.keySet().iterator();
		while(it.hasNext()){
			File[] files = mainfiles.get(it.next());
			System.out.println("files.length = "+files.length);
			for(int i=0;i<files.length;i++){
				System.out.println(files[i]);
			}
		}
				
	}

}

class MyFilenameFilter implements FilenameFilter{
	private String[] exts;
	MyFilenameFilter(String[] extensions){
		exts = extensions;
	}
    public boolean accept(File dir, String name) {
    	boolean result = false;
    	for(int i=0;i<exts.length;i++)
    		if(name.toLowerCase().endsWith(exts[i]))
    			result = true;
        return result;
    }
}

