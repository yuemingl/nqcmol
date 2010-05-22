/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.IOException;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.potential.HarmonicVibration;


/**
 *
 * @author nqc
 */
public class TaskHarmonicVibrationAnalysis extends TaskCalculate {
	@Override
	public String getName(){
		return "HarmonicVibrationAnalysis";
	}


	@Override
	protected void Process() {
		try {
			int i = 0;
			while (mol.Read(fileIn, sFormatIn))
            if(mol.getNAtoms()>2){
				//mol.Write(System.err, "xyz");
				//calculate Hessian
//				pot.getGradient(mol);
//				System.err.printf("Gradient \n");
//				MTools.PrintArray(mol.getGradient());
				pot.getNumericalHessian(mol, 1e-5);
//				System.err.printf("Hessian \n");
//				MTools.PrintArray(mol.getHessian());
				HarmonicVibration vib = new HarmonicVibration(mol);
				vib.CalcFreq(); //calculate frequencies of molecule
				//mol.Write(System.err, "xyz");
				//write a report
				xmllog.writeEntity("Cluster").writeAttribute("id", Integer.toString(i));
				xmllog.writeAttribute("nAtom", Integer.toString(mol.getNAtoms()));
				xmllog.writeAttribute("SmallestFreq", Double.toString(mol.getFreqs(0)));
				xmllog.writeAttribute("LargestFreq", Double.toString(mol.getFreqs(mol.getFreqs().length - 1)));
				xmllog.endEntity();
				xmllog.flush();

                if ((fileOut!=null)&& !mol.isNAN()) {
					mol.Write(fileOut, sFormatOut);
				}
				i++;

			}
			fileIn.close();
		} catch (IOException ex) {
			Logger.getLogger(TaskHarmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	
}
