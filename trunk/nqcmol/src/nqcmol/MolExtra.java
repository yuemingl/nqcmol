/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import org.openscience.cdk.modeling.forcefield.*;

/**
 *
 * @author nqc
 */
public class MolExtra {
	static public IPotentialFunction SetupPotential(String sPotential){
		IPotentialFunction pot=null;
		if(sPotential.contentEquals("LJ")){
			pot=new LennardJonesFunction();
		}

		if(sPotential.contentEquals("OSS2")){
			pot=new OSS2Function();
		}
		return pot;
	}
	
}
