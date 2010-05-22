/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.*;





/**
 *
 * @author nqc
 */
public class TaskOptimizeCluster extends TaskCalculate {
    @Option(name="-ftol",usage=" energy tolerance: default is 1e-9",metaVar="DOUBLE")
    double ftol=1e-9;

    @Option(name="-gtol",usage=" RMS of gradient tolerance: default is 1e-6",metaVar="DOUBLE")
    double gtol=1e-6;

    @Option(name="-mstep",usage=" Maximum step size of each iteration: default is 0.1",metaVar="DOUBLE")
    double mstep=1e-1;

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

            pot.setEnergyTol(ftol);
            pot.setGradientTol(gtol);
            pot.setMaxStepSize(mstep);
            
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
		} catch (IOException ex) {
			Logger.getLogger(TaskOptimizeCluster.class.getName()).log(Level.SEVERE, null, ex);
		}
		
	}

	
}
