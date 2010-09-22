/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
import java.io.File;
import java.io.IOException;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.Assert;
import org.junit.Before;
import org.junit.Test;

/**
 *
 * @author nqc
 */
public class ClusterTest {

    public ClusterTest() {
    }

	Cluster mol;

    @Before
    public void setUp() {
		
    }

	
	/**
	 * Test of Read method, of class Cluster.
	 */
	@Test
	public void testRead() {
		//FileReader is = null;
		Scanner is = null;
		try {
			System.out.print(" Testing Read method in Cluster \n");
			mol = new Cluster();
			//is = new FileReader("test/OSS2lm/Wn-test.xyz");
			is= new Scanner(new File("test/OSS2lm/Wn-test.xyz"));

			while(is.hasNext()){
				Assert.assertTrue(mol.Read(is,"xyz"));
			
				System.out.print(" Testing Write XYZ method in Cluster \n");
				mol.Write(System.out,"xyz");
			}
		} catch (IOException ex) {
			Logger.getLogger(ClusterTest.class.getName()).log(Level.SEVERE, null, ex);
		}  finally {
			is.close();
		}
	}

	
}