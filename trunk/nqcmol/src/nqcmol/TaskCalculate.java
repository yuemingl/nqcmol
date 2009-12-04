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

	@Option(name="-u",usage="unit of energy if applicant",metaVar="FILE")
    String sUnit="";


	protected Potential pot;

	Cluster mol = new Cluster();
	Scanner fileIn=null;
	FileWriter fileOut=null;

	@Override
	protected void Initialize() {
		super.Initialize();
		try {
			fileIn = new Scanner(new File(sFileIn));
			if(!sFileOut.isEmpty()) fileOut=new FileWriter(new File(sFileOut));
			
			pot = MolExtra.SetupPotential(sPotential, sFileParam, sUnit);
			xmllog.writeNormalText(pot.Info(1));

		} catch (IOException ex) {
			Logger.getLogger(TaskCalculateEnergy.class.getName()).log(Level.SEVERE, null, ex);
		}

	}
	
}
