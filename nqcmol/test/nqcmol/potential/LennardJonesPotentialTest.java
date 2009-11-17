/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.potential;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.Cluster;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;

/**
 *
 * @author nqc
 */
public class LennardJonesPotentialTest {
	private Cluster cluster;

    public LennardJonesPotentialTest() {
    }

	
    @Before
    public void setUp() {
		Scanner is = null;
		try {
			System.out.print(" Read Cluster \n");
			cluster = new Cluster();
			is = new Scanner(new File("test/LJlm/lj13.xyz"));
			Assert.assertTrue(cluster.Read(is,"xyz"));

			System.out.print(" Testing Write XYZ method in Cluster \n");
			cluster.Write(System.out,"xyz");
		} catch (FileNotFoundException ex) {
			Logger.getLogger(LennardJonesPotentialTest.class.getName()).log(Level.SEVERE, null, ex);
		}  finally {
			is.close();
		}
    }

 

	/**
	 * Test of Energy_ method, of class LennardJonesPotential.
	 */
	@Test
	public void testEnergy_() {
		System.out.println("Test Energy_");
		LennardJonesPotential instance = new LennardJonesPotential();
		double expResult =  -44.326801;
		double result = instance.Energy_(cluster.getCoords());
		assertEquals(expResult, result, 0.0001);
	}

	/**
	 * Test of Gradient_ method, of class LennardJonesPotential.
	 */
	@Test
	public void testGradient_() {
		System.out.println("Validate Gradient_");
		LennardJonesPotential instance = new LennardJonesPotential();
		instance.setCluster(cluster);
		System.out.println(instance.ValidateGradient(1e-4));
		//instance.Gradient_(cluster.getCoords(),cluster.getGradient());
		// TODO review the generated test code and remove the default call to fail.
		//fail("The test case is a prototype.");
	}

}