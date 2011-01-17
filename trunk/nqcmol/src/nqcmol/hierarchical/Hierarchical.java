package nqcmol.hierarchical;

import java.util.ArrayList;
import java.util.List;
import nqcmol.hierarchical.SparseMatrix.DistanceInfo;

/**
 * Agglomerative hierarchical clustering algorithm. The linkage criteria
 * is set in the constructor of the algorithm.
 *   
 * @author Ahmad Faheem & Amr Ghoneim 
 */
public class Hierarchical  {

	//Vector[] patterns;
	//LinkageCriterion linkageCriterion;
	
	List<Node> nodes;
	
	SparseMatrix sm;
	
	public Hierarchical(Vector[] patterns, LinkageCriterion linkageCriterion) {
		//this.patterns = patterns;
		//this.linkageCriterion = linkageCriterion;

		// create a cluster for each pattern.
        nodes = new ArrayList<Node>();
		for (int i = 0; i < patterns.length; i++) {
			nodes.add(new Node("C" + i, i));
		}

		sm = new SparseMatrix(patterns, linkageCriterion);
	}

    public Hierarchical(double[][] patterns, LinkageCriterion linkageCriterion) {
		//this.patterns = patterns;
		//this.linkageCriterion = linkageCriterion;

		// create a cluster for each pattern.
        nodes = new ArrayList<Node>();
		for (int i = 0; i < patterns.length; i++) {
			nodes.add(new Node("C" + i, i));
		}

		sm = new SparseMatrix(patterns, linkageCriterion);
	}

    public void setNodeName(int i,String name){
        this.nodes.get(i).setName(name);
    }
	
	public List<Integer> partition() {	
		
		int iterationCount = 0;
        String s="";

		for (DistanceInfo distanceInfo = sm.getClosestPair(); 
			distanceInfo != null; 
			distanceInfo = sm.getClosestPair()) {
			
			// we need to merge clusters distanceInfo.row and distanceInfo.column
			// into one cluster
			// this happens in two steps:
			//  1. merging the distances
			//  2. merging the clusters themselves that are stored in clusters list
			sm.merge(distanceInfo.row, distanceInfo.column);
			
			Node newCluster = new Node(nodes.get(distanceInfo.row), nodes.get(distanceInfo.column), distanceInfo.distance);
			nodes.remove(distanceInfo.column);
			nodes.remove(distanceInfo.row);
			nodes.add(distanceInfo.row, newCluster);
			
			iterationCount++;

            s=newCluster.toString();
            System.out.println(newCluster);
		}

        //System.out.println(clusters.size());
		System.out.println(s);
		
		return null;
	}
	
	public Node getRootCluster() {
		return nodes.get(0);
	}
	
	public void adjustRange() {
		double factor;
		if(this.sm.average < 10)
			factor = .001;
		else if(this.sm.average < 1000)
			factor = .01;
		else
			factor = .08;
/*		else
			factor = .01;
		if(this.sm.average < 100000)
			factor = .08;
		else
			factor = .01;
*/		getRootCluster().scale(factor);
	}
//
//	public static void main(String[] args) {
//		Vector v0 = new Vector(0);
//		Vector v1 = new Vector(1);
//		Vector v2 = new Vector(3);
//		Vector v3 = new Vector(5);
//		Vector v4 = new Vector(7);
//		Vector v5 = new Vector(10);
//
//		List<Vector> patternList = new ArrayList<Vector>();
//		patternList.add(v0);
//		patternList.add(v1);
//		patternList.add(v2);
//		patternList.add(v3);
//		patternList.add(v4);
//		patternList.add(v5);
//
//		Hierarchical hierarchical = new Hierarchical(patternList.toArray(new Vector[patternList.size()]), LinkageCriterion.SINGLE);
//		hierarchical.partition();
//
//	}
	
}
