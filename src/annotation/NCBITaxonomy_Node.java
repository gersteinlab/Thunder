package annotation;

import java.util.ArrayList;
import java.util.HashSet;
import java.util.Iterator;

public class NCBITaxonomy_Node {


	private String _nodeName; 
	private String _taxonomyName; 
	private int _nodeIndex;
	private String _nodeLevel = null;

	private int _directReadCount = 0;
	private int _inheritedReadCount = 0;

	private boolean _hasParent = false;
	private NCBITaxonomy_Node _parent = null;

	private boolean _hasChildren = false;
	private ArrayList<NCBITaxonomy_Node> _children = new ArrayList<NCBITaxonomy_Node>();

	public NCBITaxonomy_Node(String nodeName){ _nodeName = nodeName; }
	public NCBITaxonomy_Node(String nodeName, String taxonomyName, int nodeIndex){ 
		_nodeName = nodeName; 
		_taxonomyName = taxonomyName; 
		_nodeIndex = nodeIndex;
	}
	

	/**
	 * Describe this node
	 * @return
	 */
	public String getName(){ return _nodeName; }
	public int getNodeIndex(){ return _nodeIndex; }
	public String getTaxonomyName(){ return _taxonomyName; }
	public void setNodeLevel(String level){ _nodeLevel = level; }
	public String getLevel(){ return _nodeLevel; }

	/**
	 * Set parent node
	 * @param parent
	 */
	public void setParent(NCBITaxonomy_Node parent){ _parent = parent; _hasParent = true; }
	public boolean hasParent(){ return _hasParent; }
	public NCBITaxonomy_Node getParent(){ return _parent; }

	/**
	 * Deal with child nodes
	 * @param child
	 */
	public void addChild(NCBITaxonomy_Node child){ _children.add(child); _hasChildren = true; }
	public boolean hasChildren(){ return _hasChildren; }
	public int countChildren(){ return _children.size(); }
	public ArrayList<NCBITaxonomy_Node> getChildren(){ return _children; }

	/**
	 * Get set of all lowest-level 'leaf' nodes below the current node 
	 * @return
	 */
	public HashSet<String> getAllLeafSpeciesForNode(){
		return getAllLeafSpeciesForNode(this);
	}
	public static HashSet<String> getAllLeafSpeciesForNode(NCBITaxonomy_Node node){
		HashSet<String> result = new HashSet<String>();
		//if(node.hasChildren())
		//	System.out.println(node.getName()+"\t"+node.getLevel()+"\t"+node.hasChildren());

		if(node.hasChildren()){
			ArrayList<NCBITaxonomy_Node> children = node.getChildren();
			Iterator<NCBITaxonomy_Node> it = children.iterator();
			while(it.hasNext())
				result.addAll(getAllLeafSpeciesForNode(it.next()));
		}else{
			result.add(node.getName());
		}
		return result;
	}


	/**
	 * Add a read that is directly assigned to this node.  Also add it to all parents of this node 
	 * @param readID
	 */
	public void addRead(String readID){
		_directReadCount ++;
		if(hasParent()){
			_parent.addInheritedRead(readID);
			//_inheritedReadCount --;  // need to remove the read ffrom the current node to avoid double-counting
		}
	}

	/**
	 * Add a read that is assigned to a child node to allow counts to accumulate up the taxonomy
	 * @param readID
	 */
	public void addInheritedRead(String readID){
		_inheritedReadCount ++;
		if(hasParent())
			_parent.addInheritedRead(readID);
	}

	public int getReadCount(){ return _directReadCount; }
	public int getInheritedReadCount(){ return _inheritedReadCount; }



	public void printSummary(String indents, int minReads){

		if(_inheritedReadCount >= minReads  ||  _directReadCount >= minReads){
			System.out.println(indents+"\t"+(indents.length()-1)+"\t"+_nodeLevel+"\t"+_taxonomyName+"\t"+_directReadCount+"\t"+_inheritedReadCount);
			Iterator<NCBITaxonomy_Node> it = getChildren().iterator();
			while(it.hasNext())
				it.next().printSummary(indents+">", minReads);
		}
	}
}
