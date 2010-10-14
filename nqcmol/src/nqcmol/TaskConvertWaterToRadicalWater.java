/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
import java.io.*;
import java.util.Scanner;

/**
 *
 * @author nqc
 */
public class TaskConvertWaterToRadicalWater extends Task{	
	@Override
	public String getName(){
		return "ConvertWaterToRadicalWater";
	}

    static final public String Option="W2radW";

    static final public String Descriptions="\t "+Option+" \t - "+ "Convert neutral water clusters to radical water clustesr by chopping off a dangling hydrogen bond\n";

    public static void main(String[] args) throws IOException  {
        new TaskConvertWaterToRadicalWater().Execute(args);
	}

	@Override
	protected void Initialize() {
		super.Initialize();		
		//xmllog.writeAttribute("USRthreshold", Double.toString(lim));
		//xmllog.writeAttribute("deltaE", Double.toString(deltaE));
	}

	@Override
	protected void Process() {
		try {
			xmllog.writeEntity("Note");
			xmllog.writeText("Convert neutral water clusters to radical water clustesr by chopping off a dangling hydrogen bond\n");
			xmllog.endEntity();

			Scanner scanner = new Scanner(new File(sFileIn));
            if(!sFileOut.isEmpty()){
				fileOut=new FileWriter(new File(sFileOut));
			}

			int i = 0;
            Cluster mol = new Cluster();
			while (scanner.hasNext()) {
				
				mol.Read(scanner, sFormatIn);
				
				if(mol.getNAtoms()<=0) break;

                xmllog.writeEntity("Cluster");
				xmllog.writeAttribute("id", Integer.toString(i));
				xmllog.writeAttribute("Tag", mol.getTag());
				xmllog.writeAttribute("nAtoms", Integer.toString(mol.getNAtoms()));

                int count=0;
                for(int h1=0;h1<mol.getNAtoms();h1++)   if(mol.IsHydrogen(h1)){
                        boolean IsHydFar=false;
                        boolean IsHydNear=false;

                        for(int o1=0;o1<mol.getNAtoms();o1++) if(mol.IsNonHydrogen(o1)){
                            Cluster.PairwiseType bondType=mol.getPairwiseBond(h1, o1);
                            if(bondType==Cluster.PairwiseType.HYD_FAR) IsHydFar=true;
                            if(bondType==Cluster.PairwiseType.HYD_NEAR) IsHydNear=true;
                        }

                        if((IsHydFar==false)&&(IsHydNear==true)){//h1 is free OH
                            xmllog.writeEntity("FreeOH");
                                xmllog.writeAttribute("id", Integer.toString(h1));

                                if(fileOut!=null){
                                        //if(nC!=1){ //filter out free OH attached to 1-coordinate oxygen
                                        Cluster copyMol=new Cluster();
                                        copyMol.set(mol);
                                        String newTag=mol.getTag();
                                        newTag=newTag+String.format("-%d", count);
                                        copyMol.setTag(newTag);
                                        copyMol.eraseAtom(h1);
                                        copyMol.Write(fileOut, sFormatIn);
                                        count++;
                                }
                              xmllog.endEntity().flush();
                        }
                }

                xmllog.endEntity().flush();

				
				//if(i>2)pop.get(i-1).Write(System.err,MTools.getExtension(sFileOut));
				i++;				
			}			

		} catch (IOException ex) {
			//Logger.getLogger(TaskSingleCluster.class.getName()).log(Level.SEVERE, null, ex);
		}
	}


}
