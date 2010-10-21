/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.ga;

import java.io.*;
import java.util.Random;
import java.util.Scanner;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.cluster.Cluster;
import nqcmol.tools.MTools;
import org.kohsuke.args4j.Option;



/**
 *
 * @author nqc
 */
public class TaskMPIParallelTempering extends TaskMPIExecutor {
    
    @Option(name = "-seed", usage = "Seed for random number generation.[current time]", metaVar = "INTEGER")
	long seed = System.currentTimeMillis();

    @Option(name = "-T", usage = "Temperature file. If not specified, temperature of each trajectories is calculated automatically within [Tmin Tmax].", metaVar = "STRING")
	String sFileTemp="";

    @Option(name = "-Tmin", usage = "Lower bound of temperatures. [100]", metaVar = "DOUBLE")
	double Tmin=100; //!< range of temperatures

    @Option(name = "-Tmax", usage = "upper bound of temperatures. [100]", metaVar = "DOUBLE")
	double Tmax=100; //!< range of temperatures
    

    @Override
	public String getName() {
		return "MPIParalelTempering";
	}

    static final public String Option="mpipt";

    static final public String Descriptions="\t "+Option+" \t - "+ "Running MPI Parallel Tempering\n";
    
    public static void main(String[] args){
        new TaskMPIParallelTempering().Execute(args);
    }

    void copyClusterToBuffer(Cluster m,double[] dest,int destpos){
         int pos=destpos;
         dest[pos]=m.getNAtoms(); pos++;
         //dest[pos]=Double.parseDouble(m.getTag()); pos++;
         dest[pos]=m.getEnergy(); pos++;
         int[] Nz=m.getAtomicNumber();
         for(int i=0;i<Nz.length;i++){
             dest[pos]=Nz[i];  pos++;
         }

         System.arraycopy(m.getCoords(), 0, dest, pos,m.getCoords().length);       pos+=m.getCoords().length;
    }

    void copyBufferToCluster(double[] src,int srcpos,Cluster m){
        int pos=srcpos;
        m.setNAtoms((int) src[pos]); pos++;
        //m.setTag(String.format("%d",src[pos])); pos++;
        m.setEnergy(src[pos]);  pos++;
        int[] Nz=m.getAtomicNumber();
        for(int i=0;i<Nz.length;i++){
           Nz[i]=(int) src[pos];  pos++;
        }

        System.arraycopy(src, pos, m.getCoords(),0,m.getCoords().length);       pos+=m.getCoords().length;
    }

    @Override
    int iConverged(){
        throw new UnsupportedOperationException("Not yet implemented");
    }


    @Override
    public double[][] MasterInitialize(){
        //initialize randome generator for swapping
        ranSwap=new Random(seed);

        LoadTemperature();

         logger.fine(name+": LoadStructure: start reading");
         LoadInitialStructure();

         for(int i=0;i<ncpu;i++)
            logger.fine(name+": LoadStructure: process "+i+"th will get structure with tag="+pool.get(i).getTag());

         //arrange the initialized data into package
         double[][] answer=new double[ncpu][];
         for(int i=0;i<ncpu;i++){
             Cluster m=pool.get(i);
             answer[i]=new double[1+1+1+1+m.getAtomicNumber().length+m.getCoords().length];

             int pos=0;
             answer[i][pos]=T[i];
             copyClusterToBuffer(m,answer[i],1);
         }

         SlaveInitialize(answer[MASTER]);

         return answer;
    }

    @Override
    public double[] SlaveInitialize(double[] data){
        MTools.PrintArray(data);
        mc=new MonteCarloSimulator();
        mc.setName(String.format("c%2d",myrank));
        mc.Initialize(args, seed+2*(myrank+1),seed+3*(myrank+1));

        //read the initialized data from package
        int pos=0;
        mc.setTemperature(data[pos]); pos++;

        Cluster m=new Cluster();

        copyBufferToCluster(data,pos,m);

        mc.setCurrentEnsemble(m);

        mc.getCurrentEnsemble().Write(System.out, sFormatOut);

        //print out MonteCarlo info

       xmllog.writeNormalText(mc.XMLInfo(3));

        double[] answer=null;
        return answer;
    }

    @Override
    public double[][] MasterWork() {
        throw new UnsupportedOperationException("Not yet implemented");
    }

    private double[] SlaveWork() {
        throw new UnsupportedOperationException("Not yet implemented");
    }    
    
    MonteCarloSimulator mc;

    Vector<Cluster> initpop,pool;

	private double[] T; //!< array of temperatures
    
    Random ranSwap; //!< random number generator for swapping

    String getMyRank(){
        return String.format("c%2d",myrank);
    }
    
    void LoadTemperature(){
         xmllog.writeEntity("LoadTemperature");
            T=new double[ncpu];
            if(!sFileTemp.isEmpty()){ //if file for temperature exists, read it
                xmllog.writeAttribute("Input",sFileTemp);
                try {
                    Scanner fileTemp = new Scanner(new File(sFileTemp));
                    for (int i = 0; i < ncpu; i++) {
                        T[i] = fileTemp.nextDouble();
                    }
                    fileTemp.close();
                } catch (FileNotFoundException ex) {
                    Logger.getLogger(TaskMPIParallelTempering.class.getName()).log(Level.SEVERE, null, ex);
                }
            }else{
                xmllog.writeAttribute("Input","FromFormula");
                double alpha=Math.pow(Tmax/Tmin,1.0/(ncpu-1));
                for(int i=0;i<ncpu;i++) T[i]=Tmin*Math.pow(alpha,i);
            }



            for(int i=0;i<ncpu;i++){
                xmllog.writeEntity("CPU");
                    xmllog.writeAttribute("Rank",Integer.toString(i)).writeAttribute("T", Double.toString(T[i]));
                xmllog.endEntity().flush();
            }
        xmllog.endEntity().flush();
	}

	void LoadInitialStructure(){        		
        try {
            //read initial structures from file
            xmllog.writeEntity("LoadStructure");
                xmllog.writeAttribute("InputStructure", sFileIn).flush();
                fileIn = new Scanner(new File(sFileIn));

                initpop = new Vector<Cluster>();

                Cluster tmpCluster  = new Cluster();
                while (tmpCluster.Read(fileIn, sFormatIn)) {
                    //tmpCluster.Write(System.out, sFormatOut);
                    if (tmpCluster.getNAtoms() != 0) {
                        initpop.add((Cluster) tmpCluster.clone());
                    }else break;

                }
                fileIn.close();
                xmllog.writeAttribute("NumOfStruct", Integer.toString(initpop.size())).flush();

                //prepare random structures for slave
                pool=new Vector<Cluster>();
                for (int i = 0; i < ncpu; i++) {
                    int iSelect= i%initpop.size();
                    pool.add((Cluster) initpop.elementAt(iSelect));
                }
            xmllog.endEntity().flush();
        } catch (FileNotFoundException ex) {
           logger.severe(ex.getMessage());
        }
	}
}

