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
import mpi.Datatype;
import mpi.Intracomm;
import nqcmol.*;
import mpi.MPI;
import nqcmol.cluster.Cluster;
import nqcmol.tools.XmlWriter;
import org.kohsuke.args4j.Option;



/**
 *
 * @author nqc
 */
public class TaskMPIParallelTempering extends Task {
    
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

    

    @Override
    public void Execute(String[] args) {
        try {
            this.args = MPI.Init(args);
            comm=MPI.COMM_WORLD;
            myrank =comm.Rank();
            ncpu = comm.Size();
            name=String.format("c%2d",myrank);
            //logger.fine(name+": Test MPI Scatter");
            double[] a = new double[]{1225};
            double[] b = new double[]{0};
            if (myrank == MASTER) {
                logger.fine(name+": Test MPI Scatter:send data " + a[0]);
            }
            comm.Scatter(a, 0, 1, MPI.DOUBLE, b, 0, 1, MPI.DOUBLE, MASTER);
            logger.fine(name+": Test MPI Scatter: receive " + b[0]);
            stdwriter = new BufferedWriter(new FileWriter("main-" + name));

            ParseArguments(this.args);
            xmllog = new XmlWriter(stdwriter);
                xmllog.writeEntity(getName());
                Initialize();
                Process();
                Finalize();
                xmllog.endEntity().close();
            stdwriter.close();
        } catch (IOException ex) {
            Logger.getLogger(TaskMPIParallelTempering.class.getName()).log(Level.SEVERE, null, ex);
        }
    }



    @Override
    protected void Process() {
        super.Process();
    }

    @Override
    protected void Finalize() {
        super.Finalize();
        MPI.Finalize();
    }

    @Override
    protected void Initialize() {
        xmllog.writeEntity("Initialize");
            mc=new MonteCarloSimulation();
            mc.setName(String.format("c%2d",myrank));
            mc.Initialize(args, seed+2*(myrank+1),seed+3*(myrank+1));
		
//			mainLog.Get("[mess=sucessfully initilized OF, trying to evaluate initial structure ...]");
//			double ener=mc->CalcEnergy(mc->ens);
//			mainLog.Get(str(format("[mess=successfully evaluate initial structure E=%lf..]")%ener));
//            mainLog.Get(mc->Info(7),"InitInfo");
            LoadTemperature();
            LoadInitialStructure();
        xmllog.endEntity().flush();
    }
    
    //================
    static final int MASTER=0;
    int myrank=MASTER;
    String name="";
    int ncpu=0;
    Intracomm comm;
    
    MonteCarloSimulation mc;

    Vector<Cluster> initpop,pool;

	double[] T; //!< array of temperatures
    
    Random ranSwap; //!< random number generator for swapping

    Datatype newtype;
	/*!
	 *  \brief Create a custom type and register to MPI
	 *  \param bCommit true for Master, false for Slave
	 */
	int CreateMyType(boolean bCommit){
//		int n1=DGR+DGR+SIZE_MAX*DGR+SIZE_MAX*DGR+2+12;
//		int n2=1+SIZE_MAX+1+1+TYPE_MAX;

//		Datatype type1=Datatype.Contiguous(n1,MPI.DOUBLE);
//		.Datatype type2=INT.Create_contiguous(n2);
//		MPI.Datatype T[]={type1,type2};
//		int B[]={1,1};
//		Aint D[]={0,(n1)*sizeof(double)};
//		newtype=Datatype::Create_struct(2,B,D,T);
//		if(bCommit)
//			newtype.Commit();
		return 0;
	}

    void BufferToCluster(double[] data,Cluster mol){
//        mol.setNAtoms((int) data[0]);
//        mol.setEnergy(data[1]);
//        for(int i=0;i<mol.getNAtoms()*3;i++){
//				//cout<<i<<" ++ "<<xmol[i]<<"\n";
//				mol.getCoords()[i]=data[1+i];
//			}
//			mol.SetEnergy(energy);
    }

    String getMyRank(){
        return String.format("c%2d",myrank);
    }
    
    void LoadTemperature(){
         xmllog.writeEntity("LoadTemperature");
            if(myrank==MASTER){ //Master reading the temperature
                T=new double[ncpu];
                if(!sFileTemp.isEmpty()){ //if file for temperature exists, read it
                    xmllog.writeAttribute("InputTemperature",sFileTemp);
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
                    xmllog.writeAttribute("InputTemperature","FromFormula");
                    double alpha=Math.pow(Tmax/Tmin,1.0/(ncpu-1));
                    for(int i=0;i<ncpu;i++) T[i]=Tmin*Math.pow(alpha,i);
                }

                    for(int i=0;i<ncpu;i++)
                        logger.fine(name+": LoadTemperature: process "+i+"th will get T= "+T[i]);

                   logger.fine(name+": Scattering temperatures");
            }else  xmllog.writeAttribute("InputTemperature","FromMASTER");

            double[] tmpT=new double[]{0};
            comm.Scatter(T,0,1,MPI.DOUBLE,tmpT,0,1,MPI.DOUBLE,MASTER);//and scatter to everyone
            mc.setTemperature(tmpT[0]);

            xmllog.writeAttribute("Temperature",""+tmpT[0] );
        xmllog.endEntity().flush();
	}

	void LoadInitialStructure(){
        xmllog.writeEntity("LoadStructure");
		if(myrank==MASTER){
            try {
                //initialize randome generator for swapping
                ranSwap=new Random(seed);

                //read initial structures from file
                xmllog.writeAttribute("InputStructure", sFileIn);
                fileIn = new Scanner(new File(sFileIn));
                
                initpop = new Vector<Cluster>();
                while (!fileIn.hasNext()) {
                    Cluster tmpCluster = new Cluster();
                    tmpCluster.Read(fileIn, sFormatIn);
                    if (tmpCluster.getNAtoms() != 0) {
                        initpop.add(tmpCluster);
                    }else break;
                }
                fileIn.close();
                xmllog.writeAttribute("NumOfStruct", "" + initpop.size());


                //prepare random structures for slave
                pool=new Vector<Cluster>();
                for (int i = 0; i < ncpu; i++) {
                    //if(i!=0) j=ranSwap.genInt(0,initpop.size());	else j=0;
//                    if (initpop.size() == ncpu) {
//                        j = i;
//                    } else {
//                        j = 0; //if number of input structures are not equal to the number of replica, get only the first one
//                    } //if number of input structures are not equal to the number of replica, get only the first one

                    int iSelect= i%initpop.size();
                    Cluster tmpCluster = (Cluster) initpop.elementAt(iSelect).clone();
                    pool.add(tmpCluster);

                    logger.fine(name+": will send "+iSelect+" th struct to "+i);
                }
            } catch (FileNotFoundException ex) {
               logger.severe(ex.getMessage());
            }

        }else{
            xmllog.writeAttribute("InputStructure", "FromMASTER");
        }


        //scatter structures to slaves
        logger.fine(name+": scattering initial structure");

		//comm.Scatter(pool,0,1,newtype,&buff[0],1,newtype,MASTER);
		// slaves get it and correct the order
        
        logger.fine(name+": receive initial structure");
        mc.setCurrentEnsemble(null);
		//mc->ens.Write(mainLog.Stream_(),F_XYZ);
        
        xmllog.endEntity().flush();
	}



}

