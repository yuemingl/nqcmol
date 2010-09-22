/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.*;



/**
 *
 * @author nqc
 */
public class TaskOptimizeCluster extends TaskCalculate {
    @Option(name="-ftol",usage="Energy tolerance. [1e-9]",metaVar="DOUBLE")
    double ftol=1e-9;

    @Option(name="-gtol",usage="RMS of gradient tolerance.[1e-6]",metaVar="DOUBLE")
    double gtol=1e-6;

    @Option(name="-mstep",usage="Maximum step size of each iteration. [0.1]",metaVar="DOUBLE")
    double mstep=1e-1;

    @Option(name="-nOpts",usage="Maximum number of iterations in local optimization. [0]",metaVar="INTEGER")
	int nOpts=0;

	@Override
	public String getName(){
		return "Optimize";
	}

    static final public String Option="opt";

    static final public String Descriptions="\t "+Option+" \t - "+ "optimize clusters\n";

	@Override
	protected void Process() {				
        pot.setEnergyTol(ftol);
        pot.setGradientTol(gtol);
        pot.setMaxStepSize(mstep);
        pot.setMaxEvals(nOpts);

        int i = 0;
        while (mol.Read(fileIn, sFormatIn)) {
            //mol.Write(System.out,"xyz");

            Cluster oldMol=(Cluster) mol.clone();
           // for(int k=0;k<12;k++) oldMol.getUSRsig()[k]=0;
            //MTools.PrintArray(oldMol.getUSRsig());
            oldMol.CalcUSRsignature();


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
            if ((fileOut!=null)&& !mol.isNAN()) {
                mol.Write(fileOut, sFormatOut);
                //MTools.PrintArray(mol.getUSRsig());
            }
            mol.CalcUSRsignature();

            xmllog.writeEntity("Cluster").writeAttribute("id", Integer.toString(i));
            xmllog.writeAttribute("nAtoms", Integer.toString(mol.getNAtoms()));
            xmllog.writeAttribute("OldEnergy", Double.toString(oldEnergy));
            xmllog.writeAttribute("NewEnergy", Double.toString(newEnergy));
            xmllog.writeAttribute("NEvals", Integer.toString(pot.getNEvals()));
            xmllog.writeAttribute("RMSGrad", Double.toString(pot.getRMSGrad()));
            xmllog.writeAttribute("MaxRMSGrad", Double.toString(pot.getMaxGrad()));
            xmllog.writeAttribute("Similarity", String.format("%1.2f",mol.CalcUSRSimilarity(oldMol)));
            xmllog.writeAttribute("Error", Boolean.toString(mol.isNAN()));
            xmllog.endEntity();
            xmllog.flush();
            i++;
        }
        fileIn.close();
		
	}

	
}
