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



/**
 *
 * @author nqc
 */
public class TaskCalculateBindingEnergy extends TaskCalculate {
	@Option(name="-ref",usage="to calculate binding energies. It MUST be in F_XYZ format and in the following order: neutral, protonated, deprotonated. ",metaVar="FILE")
    String sFileRef="";

	@Override
	public String getName(){
		return "CalculateBindingEnergy";
	}

	@Override
	protected void Process() {
			try {
			xmllog.writeAttribute("RefFile",sFileRef);
			xmllog.writeEntity("Note");
			xmllog.writeText("Calculate binding energies of water clusters");
			xmllog.endEntity();
			xmllog.writeNormalText(pot.Info(1));

			Scanner scanRef=new Scanner(new File(sFileRef));
			Cluster[] molRef= new Cluster[3];
			double[] oldEnergyRef=new double[3];
			double[] newEnergyRef=new double[3];

			xmllog.writeEntity("ReadReference");

			for(int i=0;i<3;i++){
				molRef[i]=new Cluster();
				xmllog.flush();
				if(!molRef[i].Read(scanRef, sFormatIn)){
					Logger.getLogger(TaskCalculate.class.getName()).log(Level.SEVERE, "Error in reading ref file");
					return;
				}

				oldEnergyRef[i]=molRef[i].getEnergy();
				pot.Optimize(molRef[i]);
				newEnergyRef[i]=pot.getEnergy(true);
				xmllog.writeEntity("Reference");
				xmllog.writeAttribute("id", Integer.toString(i));
				xmllog.writeAttribute("nAtoms", Integer.toString(molRef[i].getNAtoms()));
				xmllog.writeAttribute("OldEnergy", Double.toString(oldEnergyRef[i]));
				xmllog.writeAttribute("NewEnergy", Double.toString(newEnergyRef[i]));
				xmllog.endEntity();
			}
			xmllog.endEntity();
			xmllog.flush();

			scanRef.close();

			int i = 0;
			//JSONArray jsChild2=new JSONArray();
			while (mol.Read(fileIn, sFormatIn)) {
				//mol.Write(System.out,"xyz");

				double oldEnergy=mol.getEnergy();
				pot.setCluster(mol);
				double newEnergy=pot.getEnergy(true);

				/*
				int iCharged=mol.getNAtoms()%3;
				int index=0;
				switch(iCharged){
					case 0: index=2;break;
					case 1: index=1;break;
					case 2: index=0;break;
				}

				double oldBE= oldEnergy - (mol.getNonHydrogenNum()-1)*oldEnergyRef[2] - oldEnergyRef[index];
				double newBE= newEnergy - (mol.getNonHydrogenNum()-1)*newEnergyRef[2] - newEnergyRef[index];
				*/
				int iCharged=mol.getNAtoms()%2;
				int index=iCharged;

				double oldBE= oldEnergy - (mol.getNonHydrogenNum()-1)*oldEnergyRef[0] - oldEnergyRef[index];
				double newBE= newEnergy - (mol.getNonHydrogenNum()-1)*newEnergyRef[0] - newEnergyRef[index];

				xmllog.writeEntity("Cluster");
				xmllog.writeAttribute("id", Integer.toString(i));
				xmllog.writeAttribute("nAtoms", Integer.toString(mol.getNAtoms()));
				xmllog.writeAttribute("OldEnergy", Double.toString(oldEnergy));
				xmllog.writeAttribute("OldBindingEnergy", Double.toString(oldBE));
				xmllog.writeAttribute("NewEnergy", Double.toString(newEnergy));
				xmllog.writeAttribute("NewBindingEnergy", Double.toString(newBE));
				xmllog.writeAttribute("DeltaBindingEnergy", Double.toString(newBE-oldBE));

				xmllog.endEntity();
				i++;

				//break;
			}
		} catch (IOException ex) {
			Logger.getLogger(TaskCalculate.class.getName()).log(Level.SEVERE, null, ex);
		}
	}	
}
