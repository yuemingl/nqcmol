/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;
import java.util.Vector;
import java.util.logging.*;
import nqcmol.tools.MTools;
import nqcmol.tools.XmlWriter;
import org.kohsuke.args4j.*;

/**
 *
 * @author nqc
 */
public class TaskRemoveDuplicateCluster extends Task{

	@Option(name="-lim",usage="threshold of USR similarity",metaVar="DOUBLE")
    double lim=2;

	@Option(name="-deltaE",usage="threshold of energy difference",metaVar="DOUBLE")
    double deltaE=0.1;

	@Override
	public String getName(){
		return "RemoveDuplicateCluster";
	}

	@Override
	protected void Initialize() {
		super.Initialize();		
		xmllog.writeAttribute("USRthreshold", Double.toString(lim));
		xmllog.writeAttribute("deltaE", Double.toString(deltaE));		
	}

	@Override
	protected void Process() {
		try {
			xmllog.writeEntity("Note");
			xmllog.writeText("Remove duplicate clusters based on USR\n");
			xmllog.endEntity();

			Scanner scanner = new Scanner(new File(sFileIn));

			int i = 0;
			Vector<Cluster> pop=new Vector<Cluster>();
			while (scanner.hasNext()) {
				Cluster mol = new WaterCluster();
				mol.Read(scanner, sFormatIn);
				
				if(mol.nAtoms<=0) break;
				mol.CalcUSRsignature();
				pop.add(mol);
				//if(i>2)pop.get(i-1).Write(System.err,MTools.getExtension(sFileOut));
				i++;
				
			}			
			xmllog.writeEntity("NumOfInputClusters").writeText(Integer.toString(pop.size())).endEntity();

			boolean[] bb=new boolean[pop.size()];
			for(i=0;i<pop.size();i++) bb[i]=true;
			for(i=0;i<pop.size();i++) if(bb[i]){
				int c1=1;
				for(int j=i+1;j<pop.size();j++){
					double sim=pop.get(i).CalcUSRSimilarity(pop.get(j));
					//cout<<fixed<<setprecision(2)<<sim<<" ";
					double d=Math.abs(pop.get(i).getEnergy()-pop.get(j).getEnergy());
					if(verbose>=1){
							xmllog.writeEntity("Compare").writeAttribute("i", Integer.toString(i));
							xmllog.writeAttribute("j", Integer.toString(j));
							xmllog.writeAttribute("Similarity", Double.toString(sim));
							xmllog.writeAttribute("deltaE", Double.toString(d));
							xmllog.endEntity();
						}
					
					if((sim>=lim)&&(d<deltaE)){						
						//t++;
						bb[j]=false;
						//if(pop[i].energy<pop[j].energy) bb[j]=false;
						//else bb[i]=false;
						c1++;
					}
				}
				//cout<<t<<". "<<i<<" "<<c1<<endl;
			}


			FileWriter fileOut=null;

			if(!sFileOut.isEmpty()){
				fileOut=new FileWriter(new File(sFileOut));
			}

			//dump non duplicate to file if needed
			int	count=0;
			for(i=0;i<pop.size();i++)
					if(bb[i]) {
						if(!sFileOut.isEmpty()){
							pop.get(i).setTag(count);
							//System.err.print(pop.get(i).getNAtoms());
							pop.get(i).Write(fileOut,MTools.getExtension(sFileOut));
						}
						count++;
					}else{
//						if(fnDup!=""){
//							pop[i].Write(odup,outmode);
//						}
					}
			xmllog.writeEntity("NumOfDistincClusters").writeText(Integer.toString(count)).endEntity();

		} catch (IOException ex) {
			//Logger.getLogger(TaskSingleCluster.class.getName()).log(Level.SEVERE, null, ex);
		}
	}


}
