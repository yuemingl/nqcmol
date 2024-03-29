/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
import nqcmol.cluster.MolExtra;
import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.Vector;
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

    static final public String Option="be";

    static final public String Descriptions="\t "+Option+" \t - "+ "Calculate binding energies\n";

    public static void main(String[] args) throws IOException  {
        new TaskCalculateBindingEnergy().Execute(args);
	}

	@Override
	protected void Process() {
			try {
			xmllog.writeAttribute("RefFile",sFileRef);
			xmllog.writeEntity("Note");
			xmllog.writeText("Calculate binding energies of water clusters");
			xmllog.endEntity();			

			Scanner scanRef=new Scanner(new File(sFileRef));

			ArrayList<Cluster> molRef= new ArrayList<Cluster>();
			while (scanRef.hasNext()) {
					Cluster tmpMol = new Cluster();
					if(!tmpMol.Read(scanRef, sFormatIn)){
						Logger.getLogger(TaskCalculate.class.getName()).log(Level.SEVERE, "Error in reading ref file");
						break;
					}

					tmpMol.CorrectOrder();
					if(tmpMol.getNAtoms()==0) break;
					molRef.add(tmpMol);
			}

			double[] oldEnergyRef=new double[molRef.size()];
			double[] newEnergyRef=new double[molRef.size()];

			xmllog.writeEntity("ReadReference");


			for(int i=0;i<molRef.size();i++){
				Cluster m=molRef.get(i);
				
				oldEnergyRef[i]=m.getEnergy();
				pot.Optimize(m);
				newEnergyRef[i]=pot.getEnergy(true);
				xmllog.writeEntity("Cluster");
				xmllog.writeAttribute("id", Integer.toString(i));
				xmllog.writeAttribute("Tag", m.getTag());
				xmllog.writeAttribute("nAtoms", Integer.toString(m.getNAtoms()));
				xmllog.writeAttribute("OldEnergy", Double.toString(oldEnergyRef[i]));
				xmllog.writeAttribute("NewEnergy", Double.toString(newEnergyRef[i]));
				xmllog.endEntity();
			}
			xmllog.endEntity().flush();

			scanRef.close();

			int i = 0;
			//JSONArray jsChild2=new JSONArray();

            Scanner fileIn = new Scanner(new File(sFileIn));
            Cluster mol=new Cluster();
			while (mol.Read(fileIn, sFormatIn)) {
				//mol.Write(System.out,"xyz");

				double oldEnergy=mol.getEnergy();
				pot.setCluster(mol);
				double newEnergy=pot.getEnergy(true);

				int index=MolExtra.indexForBE(mol,molRef);

				double oldBE= oldEnergy - (mol.getNonHydrogenNum()-1)*oldEnergyRef[0] - oldEnergyRef[index];
				double newBE= newEnergy - (mol.getNonHydrogenNum()-1)*newEnergyRef[0] - newEnergyRef[index];

				xmllog.writeEntity("Cluster");
				xmllog.writeAttribute("id", Integer.toString(i));
				xmllog.writeAttribute("nAtoms", Integer.toString(mol.getNAtoms()));
				xmllog.writeAttribute("OldE", Double.toString(oldEnergy));
				xmllog.writeAttribute("OldBE", Double.toString(oldBE));
				xmllog.writeAttribute("NewE", Double.toString(newEnergy));
				xmllog.writeAttribute("NewBE", Double.toString(newBE));
				xmllog.writeAttribute("DeltaBE", Double.toString(newBE-oldBE));

				xmllog.endEntity();
				i++;

				//break;
			}
            fileIn.close();
		} catch (IOException ex) {
			logger.severe(ex.getMessage());
		}
	}	
}
