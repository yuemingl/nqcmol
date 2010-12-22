/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.cluster.Cluster;


/**
 *
 * @author nqc
 */
public class TaskValidateGradient extends TaskCalculate {
	@Override
	public String getName(){
		return "ValidateGradient";
	}

    static final public String Option="grad";

    static final public String Descriptions="\t "+Option+" \t - "+ "validate gradients. Perform both analytical (if applicant) and numerical gradients.\n";

    public static void main(String[] args) throws IOException  {
        new TaskValidateGradient().Execute(args);
	}

	@Override
	protected void Process() {		
        try {
            int i = 0;
            Scanner fileIn = new Scanner(new File(sFileIn));
            Cluster mol = new Cluster();
            while (mol.Read(fileIn, sFormatIn)) {
                //mol.Write(System.out,"xyz");
                pot.setCluster(mol);
                System.out.print(pot.ValidateGradient(1.0E-4));
                i++;
            }
            fileIn.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(TaskValidateGradient.class.getName()).log(Level.SEVERE, null, ex);
        }
	}	
}
