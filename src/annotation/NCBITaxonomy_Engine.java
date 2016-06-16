package annotation;

import java.io.BufferedReader;
import java.io.FileReader;
import java.io.IOException;
import java.util.HashMap;

import org.apache.commons.cli.ParseException;

import utils.IO_utils;

public class NCBITaxonomy_Engine {

	private HashMap<String, Integer> _nodeName2nodeIndex = new HashMap<String, Integer>();
	private HashMap<Integer, NCBITaxonomy_Node> _nodeIndex2node = new HashMap<Integer, NCBITaxonomy_Node>();
	private HashMap<Integer, String> _nodeIndex2nodeName = new HashMap<Integer, String>();

	public boolean containsNode(String nodeName){ return _nodeName2nodeIndex.containsKey(nodeName); }
	public NCBITaxonomy_Node getNode(String nodeName){ return _nodeIndex2node.get(_nodeName2nodeIndex.get(nodeName)); }
	public NCBITaxonomy_Node getNode(int nodeIndex){ return _nodeIndex2node.get(nodeIndex); }
	public NCBITaxonomy_Node getRootNode(){ return _nodeIndex2node.get(1); }
	

	public NCBITaxonomy_Engine(String taxonomyBasePath) throws IOException{
		readTaxonomy(taxonomyBasePath);
	}

	public void readNames(String path) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line;
		while((line=br.readLine())!=null){
			//System.out.println(line);
			String[] bits = line.split("\\|");
			if(bits[3].trim().equals("scientific name")){
				Integer nodeIndex = Integer.valueOf(bits[0].trim());
				String nodeName_original = bits[1].trim();
				
				String nodeName = nodeName_original.toLowerCase();
				nodeName = nodeName.replace("sp.", "sp");
				nodeName = nodeName.replace("str.", "str");
				nodeName = nodeName.replace("-", " ");
				nodeName = nodeName.replace("/", " ");
				nodeName = nodeName.replace(":", " ");
				nodeName = nodeName.replace("_", " ");
				nodeName = nodeName.replace(".", " ");
				nodeName = nodeName.replace("(", "");
				nodeName = nodeName.replace(")", "");
				nodeName = nodeName.replace("#", "");
				nodeName = nodeName.trim();
				//System.out.println(nodeName +"\t"+ nodeName.replace(".", ""));
				
				
				_nodeName2nodeIndex.put(nodeName, nodeIndex);
				_nodeIndex2nodeName.put(nodeIndex, nodeName);
				_nodeIndex2node.put(nodeIndex, new NCBITaxonomy_Node(nodeName, nodeName_original, nodeIndex));
			}
		}
		br.close();
	}


	public void readNodes(String path) throws IOException{
		BufferedReader br = new BufferedReader(new FileReader(path));
		String line;
		while((line=br.readLine())!=null){
			String[] bits = line.split("\\|");

			//System.out.println(bits[0].trim()+"\t"+bits[1].trim()+"\t"+bits[2].trim());

			Integer thisIndex = Integer.valueOf(bits[0].trim());
			Integer parentIndex = Integer.valueOf(bits[1].trim());

			// set the level description for this node
			getNode(thisIndex).setNodeLevel(bits[2].trim());
			
			// if it is not root...
			if(thisIndex.intValue() > 1){
				// add parent for the child
				getNode(thisIndex).setParent(getNode(parentIndex));

				// add child for the parent
				getNode(parentIndex).addChild(getNode(thisIndex));
			}
		}
		br.close();
	}


	public void readTaxonomy(String basePath) throws IOException{
		IO_utils.printLineErr("Reading taxonomy: node names");
		//System.out.println("\n\nNAMES:");
		readNames(basePath+"/names.dmp");

		IO_utils.printLineErr("Reading taxonomy: node hierarchy");
		//System.out.println("\n\nNODES:");
		readNodes(basePath+"/nodes.dmp");

		//IO_utils.printLineErr("Done");
	}

	
	
	


	public static void main(String[] args) throws ParseException, IOException {

		String taxonomyPath = "/Users/robk/WORK/YALE_offline/ANNOTATIONS/taxdump";

		new NCBITaxonomy_Engine(taxonomyPath);

	}
}
