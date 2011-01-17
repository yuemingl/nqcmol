package nqcmol.hierarchical;

import java.util.ArrayList;
import java.util.List;

/**
 * @author Ahmad
 */
public class Node {

	public Node left = null;
	public Node right = null;
	List<Integer> patternIndexes;
	
	public double distanceBetweenLeftAndRightClusters;
	
	public String name;
	
	private static int counter = 0;
	
	public Node(String name, int patternIndex) {
		this.name = name;
		this.distanceBetweenLeftAndRightClusters = 0;
		patternIndexes = new ArrayList<Integer>();
		patternIndexes.add(patternIndex);
	}
	
	public Node(Node left, Node right, double distanceBetweenLeftAndRightClusters) {
		this.left = left;
		this.right = right;
		// this.name = "[" + left.name + "," + right.name + "]";
		this.name = "" + counter++;
		this.distanceBetweenLeftAndRightClusters = distanceBetweenLeftAndRightClusters;
		
		patternIndexes = new ArrayList<Integer>();
		for (Integer index : left.getPatternIndexes()) {
			patternIndexes.add(index);
		}
		for (Integer index : right.getPatternIndexes()) {
			patternIndexes.add(index);
		}
	}
	
	List<Integer> getPatternIndexes() {
		return patternIndexes;
	}
	
	public void scale(double factor) {
		distanceBetweenLeftAndRightClusters /= factor;
		if(left != null)
			left.scale(factor);
		if(right != null)
			right.scale(factor);
	}
	
	@Override
	public String toString() {
		return this.name = "[" + left.name + "," + right.name + ": " + distanceBetweenLeftAndRightClusters + "]";
	}

    public void setName(String name){
        this.name=name;
    }
}
