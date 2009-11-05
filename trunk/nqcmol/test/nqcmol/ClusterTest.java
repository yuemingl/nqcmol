/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.FileInputStream;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.io.InputStream;
import java.io.OutputStream;
import java.io.Writer;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.junit.AfterClass;
import org.junit.Assert;
import org.junit.Before;
import org.junit.BeforeClass;
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
		FileInputStream is = null;
		try {
			mol = new Cluster();
			is = new FileInputStream("test/OSS2lm/Wn-test.xyz");
			Assert.assertTrue(mol.Read(is,"xyz"));
			mol.Write(System.out,"xyz");
		} catch (IOException ex) {
			Logger.getLogger(ClusterTest.class.getName()).log(Level.SEVERE, null, ex);
		}  finally {
			try {
				is.close();
			} catch (IOException ex) {
				Logger.getLogger(ClusterTest.class.getName()).log(Level.SEVERE, null, ex);
			}
		}
	}

	
}