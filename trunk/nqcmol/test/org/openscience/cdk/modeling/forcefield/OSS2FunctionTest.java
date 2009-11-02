/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.openscience.cdk.modeling.forcefield;

import java.io.FileNotFoundException;
import java.io.FileReader;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.vecmath.GMatrix;
import javax.vecmath.GVector;
import org.junit.AfterClass;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.ChemFile;
import org.openscience.cdk.ChemObject;
import static org.junit.Assert.*;
import org.openscience.cdk.exception.CDKException;
import org.openscience.cdk.interfaces.IAtomContainer;
import org.openscience.cdk.interfaces.IChemSequence;
import org.openscience.cdk.io.XYZReader;

/**
 *
 * @author nqc
 */
public class OSS2FunctionTest {

    public OSS2FunctionTest() {
    }

	@BeforeClass
	public static void setUpClass() throws Exception {
	}

	@AfterClass
	public static void tearDownClass() throws Exception {
	}

	
	/**
	 * Test of energyFunctionOfAMolecule method, of class OSS2Function.
	 */
	@Test
	public void testEnergyFunctionOfAMolecule() {
		try {
			System.out.println("energyFunctionOfAMolecule");
			IAtomContainer molecule = null;
			XYZReader reader = new XYZReader(new FileReader("test/OSS2lm/pW8oss2.first"));
			ChemFile chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
			IChemSequence seq = chemFile.getChemSequence(0);
			
			molecule=seq.getChemModel(0).getMoleculeSet().getMolecule(0);
			OSS2Function instance = new OSS2Function(molecule);
			double expResult = 0.0;
			double result = instance.energyFunctionOfAMolecule(molecule);

			//assertEquals(expResult, result);
		} catch (CDKException ex) {
			Logger.getLogger(OSS2FunctionTest.class.getName()).log(Level.SEVERE, null, ex);
		} catch (FileNotFoundException ex) {
			Logger.getLogger(OSS2FunctionTest.class.getName()).log(Level.SEVERE, null, ex);
		}
		
		//assertEquals(expResult, result);
		
	}	

}