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
import nqcmol.symmetry.Symmetry;
import nqcmol.tools.XmlWriter;

/**
 *
 * @author nqc
 */
public class TaskTest extends TaskCalculate {
	@Option(name="-nRuns",usage=" number of interations",metaVar="INTNUM")
    int nRuns=1;

	@Option(name="-nScale",usage=" number of interations",metaVar="INTNUM")
    int nScale=0;

	@Option(name="-grad",usage="Benchmark gradients or not.")
    boolean isGrad= false;

	@Override
	public String getName(){
		return "Test";
	}

	@Override
	protected void Process() {
		try {
			
			//System.out.println("What the fuck!"+sFileIn);
			int i = 0;
			//JSONArray jsChild2=new JSONArray();
			while (mol.Read(fileIn, sFormatIn)) {
				Symmetry sym=new Symmetry();
				Symmetry.CalculateSymmetry(sym, mol);
				
				xmllog.writeEntity("Bench");
				xmllog.writeAttribute("id", Integer.toString(i));
				xmllog.writeAttribute("nAtoms", Integer.toString(pot.getCluster().getNAtoms()));
				
				xmllog.endEntity();
				xmllog.flush();
				i++;
			}
			fileIn.close();
		} catch (IOException ex) {
			Logger.getLogger(TaskTest.class.getName()).log(Level.SEVERE, null, ex);
		} 
	}

	
}
