/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

/**
 *
 * @author nqc
 */
public class WaterCluster extends Cluster {
	public int getChargedOxygen(boolean bUpdate){
		if(bUpdate) getPairwiseBond();

		for(int i=0;i<nAtoms;i++)	if(!IsHydrogen(i)){
			int count=0;
			for(int j=0;j<nAtoms;j++)
				if(IsHydrogen(j) && (pairwise[i][j]==PairwiseType.HYD_NEAR) )	count++;

			if(count!=2){
				return i;
			}
		}

		return -1;
	}

}
