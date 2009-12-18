/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.*;


import nqcmol.potential.*;
import nqcmol.symmetry.PointGroup;
import nqcmol.symmetry.Symmetry;
import nqcmol.tools.XmlWriter;

/**
 *
 * @author nqc
 */
public class TaskTest extends TaskCalculate {

	@Override
	public String getName(){
		return "Test";
	}

	@Override
	protected void Process() {
		try {			
			//System.out.println("What the fuck!"+sFileIn);
			int i = 0;
			
			while (mol.Read(fileIn, sFormatIn)) {
				TaskFitPotential fit=new TaskFitPotential();
				pot.setCluster(mol);
				
				fit.func.setF(pot);
				fit.optimize();
				xmllog.writeEntity("Cluster");
				xmllog.endEntity();
				xmllog.flush();
				i++;
				break;
			}
			fileIn.close();
		} catch (IOException ex) {
			Logger.getLogger(TaskTest.class.getName()).log(Level.SEVERE, null, ex);
		} 
	}

	
}
