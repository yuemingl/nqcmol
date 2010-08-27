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
import junit.framework.Assert;
import nqcmol.Cluster;
import org.junit.Before;
import org.junit.Test;
import static org.junit.Assert.*;


/**
 *
 * @author chinh
 */
public class KKYPotentialTest {

    private Cluster cluster;

    private KKYPotential pot = new KKYPotential();

    public KKYPotentialTest() {
    }

    @Before
    public void setUp(){
        Scanner is = null;
		try {
			System.out.print(" Read Cluster \n");
			cluster = new Cluster();
			is = new Scanner(new File("test/nqcmol/potential/KKY-W5-test.xyz"));
			Assert.assertTrue(cluster.Read(is,"xyz"));

			//System.out.print(" Testing Write XYZ method in Cluster \n");
			//cluster.Write(System.out,"xyz");

            System.out.print(" Testing setParam method in Potential \n");
            pot.setParam("test/nqcmol/potential/KKY-water.param");


		} catch (FileNotFoundException ex) {
			Logger.getLogger(LennardJonesPotentialTest.class.getName()).log(Level.SEVERE, null, ex);
		}  finally {
			is.close();
		}
    }
    

    /**
     * Test of Energy_ method, of class KKYPotential.
     */
    @Test
    public void testEnergy_() {
        System.out.println("Testing Energy_");
        
        double expResult = cluster.getEnergy();
        pot.setCluster(cluster);
		double result = pot.Energy_(cluster.getCoords());
        System.out.println("Expected result = "+expResult + " Calculated result = "+result);
		assertEquals(expResult, result, 0.005);
    }

    /**
     * Test of Gradient_ method, of class KKYPotential.
     */
    @Test
    public void testGradient_() {
        System.out.println("Testing Gradient_");
        double[] _p = null;
        double[] gradient = null;
        pot.setCluster(cluster);
		System.out.println(pot.ValidateGradient(1e-4));
    }    
}