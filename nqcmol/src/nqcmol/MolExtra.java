/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.potential.*;

/**
 *
 * @author nqc
 */
public class MolExtra {
	static public Potential SetupPotential(String sPotential){
		Potential pot=null;
		if(sPotential.contentEquals("LJ")){
			pot=new LennardJonesPotential();
		}

		if(sPotential.contentEquals("OSS2")){
			pot=new OSS2Potential();
		}
		return pot;
	}
	
}
