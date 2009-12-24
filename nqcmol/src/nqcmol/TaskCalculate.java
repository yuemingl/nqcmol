/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

//import com.thoughtworks.xstream.XStream;
import java.io.*;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.*;


import nqcmol.potential.*;
import nqcmol.tools.XmlWriter;
//import org.xml.sax.SAXException;
//import org.xml.sax.helpers.AttributesImpl;

/**
 *
 * @author nqc
 */
public class TaskCalculate extends Task{

	@Option(name="-ff",usage="Specify the potential",metaVar="POTENTIAL")
    String sPotential="LJ";

	@Option(name="-a",usage="parameter file for potential if applicant",metaVar="FILE")
    String sFileParam="";

	@Option(name="-u",usage="unit of energy if applicant: Hartree, kcal/mol, kJ/mol, eV, au ",metaVar="STRING")
    String sUnit="";

	@Option(name="-med",usage="optimization methods if applicant. METHOD must be one of DFPMIN(quasi-newton), CG (conjugate gradient)",metaVar="METHOD")
	String sMethod="DFPMIN";

	@Option(name="-nOpts",usage="maximum number of evaluations",metaVar="INTEGER")
	int nOpts=0;

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
		try {			
			pot = MolExtra.SetupPotential(sPotential, sFileParam, sUnit,sMethod);
			xmllog.writeNormalText(pot.Info(1));

		} catch (IOException ex) {
			Logger.getLogger(TaskCalculateEnergy.class.getName()).log(Level.SEVERE, null, ex);
		}
	}
	
}
