/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
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
	@Option(name = "-num", usage = "Number of sampling configurations along each direction. [3]", metaVar = "INTEGER")
	int num=3;

	@Option(name = "-deltaX", usage = "Maximum displacement along one direction (in Angstrom). [0.2]", metaVar = "DOUBLE")
	double deltaX=0.2;


    @Option(name = "-level", usage = "Level of combination. [1]", metaVar = "INTEGER")
	int level=1;

    @Option(name = "-m", usage = "Specify a particular normal mode for translating. Count from 0. Disable if < 0. [-1]", metaVar = "INTEGER")
	int mode=-1;

//    @Option(name = "-scaleF", usage = "Scale displacement according to k constant. if scaleF > 0, d(f)=deltaX*(k/k(scaleF))^0.5 .[-1]", metaVar = "DOUBLE")
//	double scaleF=-1;

//   @Option(name = "-header", usage = "Header file for Gaussian input", metaVar = "DOUBLE")
//	double deltaX=0.2;


	@Override
	public String getName(){
		return "GenerateCluster";
	}

    static final public String Option="aVib";

    static final public String Descriptions="\t "+Option+" \t - "+ "Generate clusters by translating along vibrational normal modes of Gaussian output files\n";

    public static void main(String[] args) throws IOException  {
        new TaskGenerateCluster().Execute(args);
	}

	@Override
	protected void Process() {				
        xmllog.writeAttribute("deltaX", Double.toString(deltaX));
        xmllog.writeAttribute("num", Integer.toString(num));
        xmllog.writeAttribute("level", Integer.toString(level));

        Cluster cluster=new Cluster();

        while(cluster.Read(fileIn, "g03")){
            xmllog.writeEntity("Cluster").writeAttribute("Tag", cluster.getTag());
            int no=GenFromNormalModes(cluster);
            xmllog.writeAttribute("NoOfStruct", Integer.toString(no));
            xmllog.endEntity().flush();
        }
	}	
	
	int GenFromNormalModes(Cluster cluster){
        int no=0;
        //creating direction vectors
        double[][] lCART=cluster.getVibrationData().getNormalModeVectors();
        double[] x0=cluster.getCoords();
        
        Vector<double[]> dirList=new Vector<double[]>();
        double[] dir;

        if(mode<0)
        for(int m1=0;m1<lCART.length;m1++){
            if(level==1){
                dir=lCART[m1].clone();
                dirList.add(dir);
    		}else	for(int m2=m1+1;m2<lCART.length;m2++)
    					if(level==2){
                            dir=new double[lCART[m1].length];
    						MTools.VEC_PLUS_VEC(dir,lCART[m1],lCART[m2],1.0,1.0);
    						dirList.add(dir);
    					}else	for(int m3=m2+1;m3<lCART.length;m3++){
                                    dir=new double[lCART[m1].length];
    								MTools.VEC_PLUS_VEC(dir,lCART[m1],lCART[m2],1.0,1.0);
    								MTools.VEC_PLUS_VEC(dir,dir,lCART[m3],1.0,1.0);
    								dirList.add(dir);
    							}
        }
        else{
            dir=lCART[mode].clone();
            dirList.add(dir);
        }

       Cluster newCluster=(Cluster) cluster.clone();
       newCluster.setVibrationData(null);
       double[] newCoords=newCluster.getCoords();
       //double[] kConst=cluster.getVibration().getForceConst();

       for(int m=0;m<dirList.size();m++){
            MTools.NORMALIZE(dirList.elementAt(m));
            for(int i=0;i<num;i++){
                //temporary scaling factor
                double beta=(2.0*i/(num-1)-1)*deltaX;
//                double maxDisp=deltaX*kConst[m]
//                double beta=(2.0*i/(num-1)-1)*maxDisp;

                MTools.VEC_PLUS_VEC(newCoords,x0,dirList.elementAt(m),1.0,beta);

                String tag=cluster.getTag()+String.format("-M=%d-beta=%1.5f",m,beta);

                newCluster.setTag(tag);


                if (fileOut != null) {
                    newCluster.Write(fileOut, sFormatOut);
                    no++;
                }
            }
        }
       return no;
      }
}
