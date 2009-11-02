/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.ArrayList;
import java.util.List;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.kohsuke.args4j.*;

import javax.vecmath.GVector;
import org.openscience.cdk.*;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.io.*;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.modeling.forcefield.*;

// SAX classes.
import org.xml.sax.*;
import org.xml.sax.helpers.*;
//JAXP 1.1
import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.stream.*;
import javax.xml.transform.sax.*;

/**
 *
 * @author nqc
 */
public class MCalculate {

	@Option(name="-ff",usage="Specify the potential",metaVar="POTENTIAL")
    String sPotential="LJ";

	@Option(name="-i",usage="input file name",metaVar="FILE")
    String sFileIn="test/LJlm/lj-pool.xyz";

	@Option(name="-o",usage="output file name",metaVar="FILE")
    String sFileOut="";

	@Option(name="-nRuns",usage=" number of interations",metaVar="INTNUM")
    int nRuns=1;

	@Option(name="-nScale",usage=" number of interations",metaVar="INTNUM")
    int nScale=0;

    @Option(name="-h",usage="Print out the help")
    boolean isHelp= false;
	
	CmdLineParser parser=null;

	public void ParseArguments(String[] args){
		parser = new CmdLineParser(this);
		try {
			parser.parseArgument(args);
		} catch (CmdLineException ex) {
			Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		}		
	}


	public void ener(String[] args)  {		
		ParseArguments(args);
		if(isHelp){
			parser.printUsage(System.out);
			return;
		}

		FileOutputStream fout = null;
		
		try {			
			StreamResult streamResult = null;
			
			if(sFileOut.isEmpty()){				
				streamResult = new StreamResult(System.out);
			}else{
				fout = new FileOutputStream(sFileOut);
				streamResult = new StreamResult(fout);
			}
				
			SAXTransformerFactory tf = (SAXTransformerFactory) SAXTransformerFactory.newInstance();
			// SAX2.0 ContentHandler.
			TransformerHandler hd = tf.newTransformerHandler();
			Transformer serializer = hd.getTransformer();
			serializer.setOutputProperty(OutputKeys.ENCODING, "ISO-8859-1");
			serializer.setOutputProperty(OutputKeys.INDENT, "yes");
			hd.setResult(streamResult);
			
			hd.startDocument();

			AttributesImpl atts = new AttributesImpl();
			atts.clear();


			XYZReader reader=new XYZReader( new FileReader(sFileIn));
			ChemFile chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
			IChemSequence seq = chemFile.getChemSequence(0);



			IPotentialFunction ljfunc = MolExtra.SetupPotential(sPotential);

			atts.addAttribute("", "", "NoOfClusters", "", Integer.toString( seq.getChemModelCount()));
			atts.addAttribute("", "", "Potential", "", ljfunc.toString());
			hd.startElement("", "", "nqc_ener", atts);
			hd.startElement("","","Note",null);
			String sTmp="Duration is measured in miliseconds. Speed is the time for one evaluation.";
			hd.characters(sTmp.toCharArray(),0,sTmp.length());
			hd.endElement("", "", "Note");

			XYZWriter xyzWriter = new XYZWriter(System.out);
			//CMLWriter xyzWriter =new CMLWriter(System.out);

			for (int i = 0; i < seq.getChemModelCount(); i++) {
				


				IMolecule mol = seq.getChemModel(i).getMoleculeSet().getMolecule(0);
				xyzWriter.write(mol);

				GVector molR  = new GVector(mol.getAtomCount()*3);
				molR = ForceFieldTools.getCoordinates3xNVector(mol);

				int myruns=(int)((nScale>0)?nRuns*Math.pow((double)nScale/mol.getAtomCount(),2):nRuns);
				
				double energy=0;

				long duration=System.currentTimeMillis();

				for(int k=0; k <myruns ; k++)	energy=ljfunc.energyFunction(molR);

				duration= System.currentTimeMillis() - duration ;
				
				atts.clear();
				atts.addAttribute("", "", "Tag", "", Integer.toString(i));
				atts.addAttribute("", "", "nRuns", "", Integer.toString(myruns));
				atts.addAttribute("","", "Energy", "", Double.toString(energy));
				atts.addAttribute("","", "Duration", "", Long.toString(duration));
				atts.addAttribute("","", "Speed", "", Double.toString((double)duration/myruns));


//				GeometricMinimizer gm = new GeometricMinimizer();
//
//				gm.setConvergenceParametersForSDM(100, 0.00001);
//		 		gm.steepestDescentsMinimization(molR, ljfunc);
//				molR.set(gm.getSteepestDescentsMinimum());
//
//				energy=ljfunc.energyFunction(molR);
//				atts.addAttribute("","", "Mininum", "", Double.toString(energy));


				hd.startElement("", "", "bench", atts);
				hd.endElement("", "", "bench");
			}

			xyzWriter.close();
			hd.endElement("", "", "nqc_ener");
			hd.endDocument();
		} catch (IOException ex) {
			Logger.getLogger(NQCMol.class.getName()).log(Level.SEVERE, null, ex);
		} catch (TransformerConfigurationException ex) {
				Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		} catch (CDKException ex) {
				Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		} catch (SAXException ex) {
				Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		} finally {
			try {
				fout.close();
			} catch (IOException ex) {
				Logger.getLogger(NQCMol.class.getName()).log(Level.SEVERE, null, ex);
			}
		}

	}


}
