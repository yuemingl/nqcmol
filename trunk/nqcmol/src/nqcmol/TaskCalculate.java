/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

//import com.thoughtworks.xstream.XStream;
import nqcmol.cluster.Cluster;
import nqcmol.cluster.MolExtra;
import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.*;


import nqcmol.potential.*;
//import org.xml.sax.SAXException;
//import org.xml.sax.helpers.AttributesImpl;

/**
 *
 * @author nqc
 */
public class TaskCalculate extends Task{

	@Option(name="-ff",usage="Potential model. [LJ}",metaVar="POTENTIAL")
    String sPotential="LJ";

	@Option(name="-a",usage="Parameter file for potential if applicant.",metaVar="FILE")
    String sFileParam="";

	@Option(name="-u",usage="Unit of energy if applicant: Hartree, kcal/mol, kJ/mol, eV, au.",metaVar="STRING")
    String sUnit="";

	@Option(name="-med",usage="Optimization methods if applicant. METHOD must be one of DFPMIN(quasi-newton), CG (conjugate gradient)",metaVar="METHOD")
	String sMethod="DFPMIN";

	

	protected Potential pot;

	Cluster mol = new Cluster();

//	@Override
//	public void ParseArguments(String[] args){
//		super.ParseArguments(args);
//	}

	@Override
	public String getName() {
		return "FitPotential";
	}

	@Override
	protected void Initialize() {
		super.Initialize();
        pot = MolExtra.SetupPotential(sPotential, sFileParam, sUnit,sMethod);
        xmllog.writeNormalText(pot.Info(1));
	}
	
}
