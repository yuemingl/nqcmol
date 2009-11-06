/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.kohsuke.args4j.*;


// SAX classes.
import org.xml.sax.*;
import org.xml.sax.helpers.*;
//JAXP 1.1
import javax.xml.transform.*;
import javax.xml.transform.stream.*;
import javax.xml.transform.sax.*;

import nqcmol.potential.*;

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

			//IPotentialFunction ljfunc = MolExtra.SetupPotential(sPotential);
			//OSS2Function ljfunc=new OSS2Function();

			Potential pot=MolExtra.SetupPotential(sPotential);

			atts.addAttribute("", "", "Potential", "", pot.getEquation());
			hd.startElement("", "", "nqc_ener", atts);
			hd.startElement("","","Note",null);
			String sTmp="Duration is measured in miliseconds. Speed is the time for one evaluation.";
			hd.characters(sTmp.toCharArray(),0,sTmp.length());
			hd.endElement("", "", "Note");


			Cluster mol=new Cluster();


			Scanner scanner= new Scanner(new File(sFileIn));
			int i=0;
			while(mol.Read(scanner, "xyz")){
				//mol.Write(System.out,"xyz");

				int myruns=(int)((nScale>0)?nRuns*Math.pow((double)nScale/mol.getNAtoms(),2):nRuns);
				
				double energy=0;

				long duration=System.currentTimeMillis();

				pot.setCluster(mol);

				for(int k=0; k <myruns ; k++)	energy=pot.getEnergy(true);

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
				i++;
				//break;
			}
			hd.endElement("", "", "nqc_ener");
			hd.endDocument();
		} catch (IOException ex) {
			Logger.getLogger(NQCMol.class.getName()).log(Level.SEVERE, null, ex);
		} catch (TransformerConfigurationException ex) {
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
