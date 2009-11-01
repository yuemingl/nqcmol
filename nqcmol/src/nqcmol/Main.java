/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.Scanner;
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

				ener();
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

	private static void ener() throws TransformerConfigurationException, SAXException {
		FileOutputStream writer = null;
		try {
			//writer = new FileWriter("test/test.xml");
			writer = new FileOutputStream("test/test.xml");

			StreamResult streamResult = new StreamResult(writer);
			SAXTransformerFactory tf = (SAXTransformerFactory) SAXTransformerFactory.newInstance();
			// SAX2.0 ContentHandler.
			TransformerHandler hd = tf.newTransformerHandler();
			Transformer serializer = hd.getTransformer();
			serializer.setOutputProperty(OutputKeys.ENCODING, "ISO-8859-1");
			serializer.setOutputProperty(OutputKeys.DOCTYPE_SYSTEM, "users.dtd");
			serializer.setOutputProperty(OutputKeys.INDENT, "yes");
			hd.setResult(streamResult);
			hd.startDocument();
			AttributesImpl atts = new AttributesImpl();
// USERS tag.
			hd.startElement("", "", "USERS", atts);
// USER tags.
			String[] id = {"PWD122", "MX787", "A4Q45"};
			String[] type = {"customer", "manager", "employee"};
			String[] desc = {"Tim@Home", "Jack&Moud", "John D'oé"};
			for (int i = 0; i < id.length; i++) {
				atts.clear();
				atts.addAttribute("", "", "ID", "CDATA", id[i]);
				atts.addAttribute("", "", "TYPE", "CDATA", type[i]);
				hd.startElement("", "", "USER", atts);
				hd.characters(desc[i].toCharArray(), 0, desc[i].length());
				hd.endElement("", "", "USER");
			}
			hd.endElement("", "", "USERS");
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
