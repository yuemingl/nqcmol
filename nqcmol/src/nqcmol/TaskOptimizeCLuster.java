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
import nqcmol.tools.XmlWriter;


/**
 *
 * @author nqc
 */
public class TaskOptimizeCLuster extends TaskCalculate {

	@Override
	public String getName(){
		return "Optimize";
	}

	@Override
	protected void Process() {
		try {
//			xmllog.writeEntity("Note");
//			xmllog.writeText(" Time is measured in seconds");
//			xmllog.endEntity();

			int i = 0;
			while (mol.Read(fileIn, sFormatIn)) {
				//mol.Write(System.out,"xyz");
				double oldEnergy = mol.getEnergy();
				pot.Optimize(mol);
				double newEnergy = mol.getEnergy();
				//				MultivariateRealOptimizer opt=new NelderMead();
				//				PotentialFitting fit=new PotentialFitting();
				//				fit.setF(pot);
				//				double newEnergy=0;
				//				RealPointValuePair r = null;//=new RealPointValuePair();
				//				try {
				//					r = opt.optimize(fit, GoalType.MINIMIZE, mol.getCoords());
				//					newEnergy =r.getValue();
				//				} catch (FunctionEvaluationException ex) {
				//					Logger.getLogger(TaskSingleCluster.class.getName()).xmllog(Level.SEVERE, null, ex);
				//				} catch (OptimizationException ex) {
				//					Logger.getLogger(TaskSingleCluster.class.getName()).xmllog(Level.SEVERE, null, ex);
				//				}
				if (fileOut!=null) {
					mol.Write(fileOut, sFormatOut);
				}
				xmllog.writeEntity("Cluster").writeAttribute("id", Integer.toString(i));
				xmllog.writeAttribute("nAtoms", Integer.toString(mol.getNAtoms()));
				xmllog.writeAttribute("OldEnergy", Double.toString(oldEnergy));
				xmllog.writeAttribute("NewEnergy", Double.toString(newEnergy));
				xmllog.writeAttribute("NEvals", Integer.toString(pot.getNEvals()));
				xmllog.writeAttribute("RMSGrad", Double.toString(pot.getRMSGrad()));
				xmllog.writeAttribute("MaxRMSGrad", Double.toString(pot.getMaxGrad()));
				xmllog.endEntity();
				xmllog.flush();
				i++;
			}
			fileIn.close();
		} catch (IOException ex) {
			Logger.getLogger(TaskOptimizeCLuster.class.getName()).log(Level.SEVERE, null, ex);
		}
		
	}

	
}
