/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.WaterCluster;
import nqcmol.cluster.Cluster;
import java.io.*;
import java.util.ArrayList;
import java.util.Scanner;
import java.util.Vector;
import nqcmol.tools.MTools;
import org.kohsuke.args4j.*;

/**
 *
 * @author nqc
 */
public class TaskRemoveDuplicateCluster extends Task{

	@Option(name="-lim",usage="Threshold of USR similarity",metaVar="DOUBLE")
    double lim=2;

	@Option(name="-deltaE",usage="Threshold of energy difference",metaVar="DOUBLE")
    double deltaE=0.1;

    @Option(name="-spanE",usage="Energy span. Zero value is the lowest energy of isomers. All isomers out of span will be removed. Negative means not used. [-1] ",metaVar="DOUBLE")
    double spanE=-1;

    @Option(name="-retag",usage="Change tags of structures to 0,1,2,etc")
    boolean isReTag=false;

    @Option(name="-sortatom",usage="Sort atoms in descending order respect to atomic number.")
    boolean isSortAtom=false;

    @Option(name="-prefix",usage="Sort atoms in descending order respect to atomic number.",metaVar="STRING")
    String prefix="";

	@Override
	public String getName(){
		return "RemoveDuplicateCluster";
	}

    static final public String Option="sim";

    static final public String Descriptions="\t "+Option+" \t - "+ "Remove duplicate clusters based on USR\n";

    public static void main(String[] args) throws IOException  {
        new TaskRemoveDuplicateCluster().Execute(args);
	}

	@Override
	protected void Initialize() {
		super.Initialize();		
		xmllog.writeAttribute("USRthreshold", Double.toString(lim));
		xmllog.writeAttribute("deltaE", Double.toString(deltaE));
        xmllog.writeAttribute("spanE", Double.toString(spanE));
        xmllog.writeAttribute("reTag", Boolean.toString(isReTag));
        xmllog.writeAttribute("SortAtom", Boolean.toString(isSortAtom));        
	}

	@Override
	protected void Process() {
		try {
			xmllog.writeEntity("Note");
			xmllog.writeText(getName());
			xmllog.endEntity();

			Scanner fileIn = new Scanner(new File(sFileIn));

			int i = 0;
			ArrayList<Cluster> pop=new ArrayList<Cluster>();
            double minE=1e100;
			while (fileIn.hasNext()) {
				Cluster mol = new Cluster();
				mol.Read(fileIn, sFormatIn);
				
				if(mol.getNAtoms()<=0) break;
                minE=Math.min(mol.getEnergy(),minE);
				mol.CalcUSRsignature();
				pop.add(mol);
				//if(i>2)pop.get(i-1).Write(System.err,MTools.getExtension(sFileOut));
				i++;
				
			}
            fileIn.close();

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
            int countOutOfRange=0;
			for(i=0;i<pop.size();i++)
					if(bb[i]) {
						if(!sFileOut.isEmpty()){
							//if(Integer.parseInt(pop.get(i).getTag())==0)
							if(isReTag)	pop.get(i).setTag(Integer.toString(count));
                            if(isSortAtom) pop.get(i).CorrectOrder();
                            if(!prefix.isEmpty()) pop.get(i).setTag(prefix+"-"+i);
							//System.err.print(pop.get(i).getNAtoms());
                            if( ((spanE>0)&&(pop.get(i).getEnergy()<minE+spanE)) || (spanE<0)) {
                                pop.get(i).Write(fileOut,MTools.getExtension(sFileOut));                                
                            }else{
                                countOutOfRange++;
                            }
						}
                        count++;
					}else{
//						if(fnDup!=""){
//							pop[i].Write(odup,outmode);
//						}
					}
            
            fileOut.close();
			xmllog.writeEntity("Summary");
                xmllog.writeAttribute("NumOfDistincClusters",Integer.toString(count));
                xmllog.writeAttribute("NumOfOutOfRange",Integer.toString(countOutOfRange));
            xmllog.endEntity();

		} catch (IOException ex) {
			logger.severe(ex.getMessage());
		}
	}


}
