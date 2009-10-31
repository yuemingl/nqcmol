/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.openscience.cdk.modeling.forcefield;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.vecmath.GMatrix;
import javax.vecmath.GVector;
import org.junit.After;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
import org.junit.Test;
import org.openscience.cdk.*;
import static org.junit.Assert.*;
import org.openscience.cdk.interfaces.*;
import org.openscience.cdk.io.XYZReader;

/**
 *
 * @author nqc
 */
public class LennardJonesFunctionTest {

	IMolecule mol= null;

	LennardJonesFunction ljfunc = new LennardJonesFunction();

	GVector molecule3Coordinates  = null;

    public LennardJonesFunctionTest() {
    }
	

    @Before
    public void setUp() throws Exception {		
			FileReader buff = new FileReader("test/LJlm/lj13.xyz");
			XYZReader reader=new XYZReader(buff);
			ChemFile chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
			IChemSequence seq = chemFile.getChemSequence(0);
			IChemModel model = seq.getChemModel(0);
			IMoleculeSet som = model.getMoleculeSet();
			mol= som.getMolecule(0);

			System.out.println("\n Test function with the following structures\n");
			System.out.println(mol);

			molecule3Coordinates  = new GVector(mol.getAtomCount()*3);
			molecule3Coordinates = ForceFieldTools.getCoordinates3xNVector(mol);
    }

    @After
    public void tearDown() {
		System.out.println("\n Finish the test! \n");
    }

	/**
	 * Test of energyFunctionOfAMolecule method, of class LennardJonesFunction.
	 */
	@Test
	public void testEnergyFunctionOfAMolecule() {
		System.out.println("energyFunctionOfAMolecule");
		double expResult = -44.326801;
		double result = ljfunc.energyFunctionOfAMolecule(mol);

		System.out.printf(" Calculated energy = %f . Expected Results = %f \n",result,expResult);
		Assert.assertEquals(expResult, result,0.001);


		System.out.println("setEnergyGradient");
		ljfunc.setEnergyGradient(molecule3Coordinates);

		System.out.println("set2ndOrderErrorApproximateGradient");
		ljfunc.set2ndOrderErrorApproximateGradient(molecule3Coordinates);
	}	

}