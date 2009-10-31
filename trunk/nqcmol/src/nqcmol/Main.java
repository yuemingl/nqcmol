/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.openscience.cdk.*;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.io.*;
import org.openscience.cdk.interfaces.*;


/**
 *
 * @author nqc
 */
public class Main {

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) {
        Scanner scan;
        try {
            // TODO code application logic here
            scan = new Scanner(new File("test/LJlm/lj13.xyz"));
	    while(scan.hasNext()){
		System.out.print(scan.nextLine()+"\n");
	    }
        } catch (FileNotFoundException ex) {
            Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
        }
	
	try {
	    BufferedReader buff = new BufferedReader(new FileReader("test/LJlm/lj13.xyz"));
	   
	    XYZReader reader=new XYZReader(buff);
	    ChemFile chemFile;
	    
		chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
		IChemSequence seq = chemFile.getChemSequence(0);
		IChemModel model = seq.getChemModel(0);
		IMoleculeSet som = model.getMoleculeSet();
		IMolecule mol= som.getMolecule(0);

		System.out.print(mol);


		FileWriter writer = new FileWriter("test/LJlm/new-lj13.xyz");
		XYZWriter xyzWriter = new XYZWriter(writer);
		xyzWriter.write(mol);
		xyzWriter.close();
		writer.close();
	   
	} catch (FileNotFoundException ex) {
	    Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
	} catch (CDKException ex) {
		Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
	}
	catch (IOException ex) {
	    Logger.getLogger(Main.class.getName()).log(Level.SEVERE, null, ex);
	}

	
    }

}
