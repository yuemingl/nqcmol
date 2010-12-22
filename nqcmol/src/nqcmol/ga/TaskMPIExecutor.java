/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.ga;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import mpi.Intracomm;
import nqcmol.*;
import mpi.MPI;
import mpi.Status;
import nqcmol.tools.MTools;
import nqcmol.tools.XmlWriter;



/**
 *
 * @author nqc
 */
public class TaskMPIExecutor extends Task {    
    
    @Override
	public String getName() {
		return "MPIExecutor";
	}
//
//    static final public String Option="mpipt";
//
//    static final public String Descriptions="\t "+Option+" \t - "+ "Running MPI Executor\n";
//

    @Override
    public void Execute(String[] args) {
        try {
            this.args = MPI.Init(args);
            comm=MPI.COMM_WORLD;
            myrank =comm.Rank();
            ncpu = comm.Size();
            name=String.format("c%02d",myrank);
                logger.fine(name+": Test MPI Scatter");
                double[] a = new double[]{1225,234};
                double[] b = new double[]{0};
                if (myrank == MASTER) {
                    logger.fine(name+": Test MPI Scatter:send data " + a[0]);
                }
                comm.Scatter(a, 0, 1, MPI.DOUBLE, b, 0, 1, MPI.DOUBLE, MASTER);
                logger.fine(name+": Test MPI Scatter: receive " + b[0]);
                
            stdwriter = new BufferedWriter(new FileWriter("main-" + name+".xml"));

            ParseArguments(this.args);
            
            xmllog = new XmlWriter(stdwriter);
                xmllog.writeEntity(getName());
                Initialize();
//                Process();
                Finalize();
                xmllog.endEntity().close();
            stdwriter.close();
            MPI.Finalize();
        } catch (IOException ex) {
            Logger.getLogger(TaskMPIExecutor.class.getName()).log(Level.SEVERE, null, ex);
        }
    }



    @Override
    protected void Process() {  
        //send data and process
        while(true){
            double[] data=SlaveWork();
            double[][] poolData=null;

            if( IsMaster()) poolData=new double[ncpu][data.length];
            comm.Gather(data, 0, data.length,MPI.DOUBLE, poolData, 0, data.length, MPI.DOUBLE, MASTER);//master gather the data from slaves

            if( IsMaster()){
                poolData=MasterWork();
            }

            comm.Scatter(poolData,0, data.length,MPI.DOUBLE, data, 0, data.length, MPI.DOUBLE,MASTER);//then scatter back to slaves
            break;
        }
    }

    @Override
    protected void Initialize() {
        xmllog.writeEntity("Initialize");
            //initialize
            if( IsMaster()){
                double[][] poolData=MasterInitialize();
                //MTools.PrintArray(poolData);

                logger.fine(name+": scattering data for iniatialization");
                for(int i=0;i<ncpu;i++)
                    if(i!=myrank)
                        comm.Send(poolData[i], 0, poolData[i].length, MPI.DOUBLE, i, 0);
            }else{
                Status prob=comm.Probe(MASTER,0);
                double[] data=new double[prob.Get_count(MPI.DOUBLE)];
                comm.Recv(data, 0, data.length,  MPI.DOUBLE, MASTER, 0);
                SlaveInitialize(data);
            }
            //MTools.PrintArray(poolData);
//            double[] data=new double[poolData[0].length];
//            comm.Scatter(poolData,0, data.length,MPI.DOUBLE, data, 0, data.length, MPI.DOUBLE,MASTER);

            //SlaveInitialize(data);
        xmllog.endEntity().flush();
    }
    
    //================
    static final int MASTER=0;

    int myrank=MASTER;

    String name="";

    int ncpu=0;

    Intracomm comm;

    boolean IsMaster(){
        return myrank==MASTER;
    }

    int iConverged(){
        throw new UnsupportedOperationException("Not yet implemented");
    }
    

    public double[][] MasterInitialize(){
        throw new UnsupportedOperationException("Not yet implemented");
    }

    public double[] SlaveInitialize(double[] data){
        throw new UnsupportedOperationException("Not yet implemented");
    }

    public double[][] MasterWork() {
        throw new UnsupportedOperationException("Not yet implemented");
    }

    private double[] SlaveWork() {
        throw new UnsupportedOperationException("Not yet implemented");
    }
}

