/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
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

import org.apache.commons.cli.*;


/**
 *
 * @author nqc
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {
       	try {
            Options opt = new Options();
			CommandLineParser parser = new GnuParser();
			CommandLine cl = parser.parse( opt, args);

            if ( cl.hasOption('h') ) {
                HelpFormatter f = new HelpFormatter();
                f.printHelp("OptionsTip", opt);
            }
            else {
               if(cl.hasOption("Tener")){
				
			   }
                //System.out.println(cl.getOptionValue("dsn"));
					 ener(args);
            }
        }
        catch (ParseException e) {
            e.printStackTrace();
        }


			//BufferedReader buff = new BufferedReader(new FileReader("test/LJlm/lj13.xyz"));
//			FileReader buff = new FileReader("test/LJlm/lj13.xyz");
//			XYZReader reader=new XYZReader(buff);
//			ChemFile chemFile;
//
//			chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
//			IChemSequence seq = chemFile.getChemSequence(0);
//			IChemModel model = seq.getChemModel(0);
//			IMoleculeSet som = model.getMoleculeSet();
//			IMolecule mol= som.getMolecule(0);
//
//			System.out.print(mol);
//
//
//			FileWriter writer = new FileWriter("test/LJlm/new-lj13.xyz");
//			XYZWriter xyzWriter = new XYZWriter(writer);
//			xyzWriter.write(mol);
//			xyzWriter.close();
//			writer.close();
//
//			System.out.println("\n Now evaluate the energy \n");
//			LennardJonesFunction pot = new LennardJonesFunction();
//			double energy = pot.energyFunctionOfAMolecule(mol);
//
//			System.out.println(" Energy = "+ Double.toString(energy));
//
//			GVector molecule3Coordinates  = new GVector(mol.getAtomCount()*3);
//			molecule3Coordinates=ForceFieldTools.getCoordinates3xNVector(mol);
			//GeometricMinimizer gm = new GeometricMinimizer();
			//gm.setConvergenceParametersForCGM(100, 0.00001);
	 		//gm.conjugateGradientMinimization(molecule3Coordinates, pot);

    }

	private static void ener(String[] args) throws TransformerConfigurationException, SAXException, CDKException {
		FileOutputStream writer = null;
		try {			
			writer = new FileOutputStream("test/test.xml");

			StreamResult streamResult = new StreamResult(writer);
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
			atts.addAttribute("", "", "Potential", "", "Lennard-Jones");
			

			XYZReader reader=new XYZReader( new FileReader("test/LJlm/lj-pool.xyz"));
			ChemFile chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
			IChemSequence seq = chemFile.getChemSequence(0);
	

			IPotentialFunction ljfunc = new LennardJonesFunction();

			atts.addAttribute("", "", "NoOfClusters", "", Integer.toString( seq.getChemModelCount()));
			atts.addAttribute("", "", "Potential", "", ljfunc.toString());
			hd.startElement("", "", "nqc_ener", atts);

			XYZWriter xyzWriter = new XYZWriter(System.out);

			
			for (int i = 0; i < seq.getChemModelCount(); i++) {

				atts.clear();
				atts.addAttribute("", "", "id", "", Integer.toString(i));
				

				IMolecule mol = seq.getChemModel(i).getMoleculeSet().getMolecule(0);
				xyzWriter.write(mol);

				GVector molR  = new GVector(mol.getAtomCount()*3);
				molR = ForceFieldTools.getCoordinates3xNVector(mol);
				
				double energy=ljfunc.energyFunction(molR);
				atts.addAttribute("","", "Energy", "", Double.toString(energy));


				GeometricMinimizer gm = new GeometricMinimizer();
				//gm.setConvergenceParametersForCGM(100, 0.00001);
				//gm.conjugateGradientMinimization(molR, ljfunc);
				//molR.set(gm.getConjugateGradientMinimum());

				gm.setConvergenceParametersForSDM(100, 0.00001);
		 		gm.steepestDescentsMinimization(molR, ljfunc);
				molR.set(gm.getSteepestDescentsMinimum());
				
				energy=ljfunc.energyFunction(molR);
				atts.addAttribute("","", "Mininum", "", Double.toString(energy));


				hd.startElement("", "", "bench", atts);
				hd.endElement("", "", "bench");
			}

			xyzWriter.close();
			hd.endElement("", "", "nqc_ener");
			hd.endDocument();
		} catch (IOException ex) {
			Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
		} finally {
			try {
				writer.close();
			} catch (IOException ex) {
				Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
			}
		}

	}
}
