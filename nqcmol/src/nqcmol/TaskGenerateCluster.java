/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import nqcmol.tools.MTools;
import org.kohsuke.args4j.*;


/**
 *
 * @author nqc
 */
public class TaskGenerateCluster extends Task {	
	@Option(name = "-num", usage = "Number of sampling configurations along each direction", metaVar = "INTEGER")
	int num=3;

	@Option(name = "-deltaX", usage = "Maximum displacement along one direction (in Angstrom)", metaVar = "DOUBLE")
	double deltaX=0.2;

//    @Option(name = "-header", usage = "Header file for Gaussian input", metaVar = "DOUBLE")
//	double deltaX=0.2;


	@Override
	public String getName(){
		return "GenerateCluster";
	}

    static final public String Option="alongVib";

    static final public String Descriptions="\t "+Option+" \t - "+ "Generate clusters by translating along vibrational normal modes\n";
	
	@Override
	protected void Process() {			
		try {			
			xmllog.writeAttribute("deltaX", Double.toString(deltaX));
			xmllog.writeAttribute("num", Integer.toString(num));

            Cluster cluster=new Cluster();

            while(cluster.Read(fileIn, "g03")){
				xmllog.writeEntity("Cluster").writeAttribute("Tag", cluster.getTag());
				GenFromNormalModes(cluster);
				xmllog.endEntity().flush();
			}
		} catch (IOException ex) {
			Logger.getLogger(TaskGenerateCluster.class.getName()).log(Level.SEVERE, null, ex);
		}
			
	}	
	
	void GenFromNormalModes(Cluster cluster){
        //creating direction vectors
        double[][] lCART=cluster.getNormalModeVectors();
        double[] x0=cluster.getCoords();
        
        Vector<double[]> dirList=new Vector<double[]>();

        for(int m1=0;m1<lCART.length;m1++){
            //if(level==1){
            double[] dir=lCART[m1].clone();
            dirList.add(dir);
    //		}else	for(int m2=m1+1;m2<k.nFreq;m2++)
    //					if(level==2){
    //						VEC_PLUS_VEC(vec,k.lCART[m1],k.lCART[m2],k.nAtom*DGR,1.0,1.0);
    //						vecList.push_back(vec);
    //					}else	for(int m3=m2+1;m3<k.nFreq;m3++){
    //								VEC_PLUS_VEC(vec,k.lCART[m1],k.lCART[m2],k.nAtom*DGR,1.0,1.0);
    //								VEC_PLUS_VEC(vec,vec,k.lCART[m3],k.nAtom*DGR,1.0,1.0);
    //								vecList.push_back(vec);
    //							}
        }

       Cluster newCluster=(Cluster) cluster.clone();
       newCluster.setFreqs(null);
       double[] newCoords=newCluster.getCoords();

        for(int m=0;m<dirList.size();m++){
            MTools.NORMALIZE(dirList.elementAt(m));
            for(int i=0;i<num;i++){
                //temporary scaling factor
                double beta=(2.0*i/(num-1)-1)*deltaX;

                MTools.VEC_PLUS_VEC(newCoords,x0,dirList.elementAt(m),1.0,beta);

                String tag=cluster.getTag()+String.format("-M=%d-beta=%1.5f",m,beta);

                newCluster.setTag(tag);


                if (fileOut != null) {
                    newCluster.Write(fileOut, sFormatOut);
                }
            }
        }
        }
}
