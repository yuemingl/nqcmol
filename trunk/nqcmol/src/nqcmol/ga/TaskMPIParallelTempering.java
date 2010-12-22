/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nqcmol.ga;

import java.io.*;
import java.util.ArrayList;
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
    String sFileTemp = "";

    @Option(name = "-Tmin", usage = "Lower bound of temperatures in Kelvin. [100]", metaVar = "DOUBLE")
    double Tmin = 100; //!< range of temperatures

    @Option(name = "-Tmax", usage = "Upper bound of temperatures in Kelvin. [100]", metaVar = "DOUBLE")
    double Tmax = 300; //!< range of temperatures

    @Option(name = "-nRuns", usage = "Length of equibrium MC (with exchange replica),  counted by MC pass. [100]", metaVar = "INTEGER")
    int nRuns = 100;

    @Option(name = "-nWarmUps", usage = "Length of warmup MC (without exchange replica),  counted by MC pass. [0]", metaVar = "INTEGER")
    int nWarmUps = 0;

    @Option(name = "-nSwaps", usage = "Determine how frequent for swapping, counted by MCpass. [10]", metaVar = "DOUBLE")
    int nSwaps = 10;//!< determine how frequent for swapping, counted by MCpass

    @Option(name = "-nMaxData", usage = "Maximum number of data sampled, counted by MC pass. when reached, some data will be discarded (determined by throw)[100]", metaVar = "INTEGER")
    int nMaxData; //!< when the sampled data reaches nMaxData, Cv and Histogram will be reported and then discard some data (determined by nDataThrow)

    @Option(name = "-throw", usage = "Amount of data will be discarded for every cycle, if throw=0, discard all. if throw=<1, throw*TotalPoint will be discarded. [1]", metaVar = "DOUBLE")
    double fDataThrow = 1;

    @Option(name = "-structSampleRate", usage = "Amount of data will be discarded for every cycle, if throw=0, discard all. if throw=<1, throw*TotalPoint will be discarded. [1]", metaVar = "DOUBLE")
    double structSampleRate = 0.01;

     @Option(name = "-dataSampleRate", usage = "Amount of data will be discarded for every cycle, if throw=0, discard all. if throw=<1, throw*TotalPoint will be discarded. [1]", metaVar = "DOUBLE")
    double dataSampleRate = 0.01;

    @Override
    public String getName() {
        return "MPIParalelTempering";
    }

    static final public String Option = "mpipt";
    static final public String Descriptions = "\t " + Option + " \t - " + "Running MPI Parallel Tempering\n";

    public static void main(String[] args) {
        new TaskMPIParallelTempering().Execute(args);
    }

    void copyClusterToBuffer(Cluster m, double[] dest, int destpos) {
        int pos = destpos;
        dest[pos] = m.getNAtoms();
        pos++;
        //dest[pos]=Double.parseDouble(m.getTag()); pos++;
        dest[pos] = m.getEnergy();
        pos++;
        int[] Nz = m.getAtomicNumber();
        for (int i = 0; i < Nz.length; i++) {
            dest[pos] = Nz[i];
            pos++;
        }

        System.arraycopy(m.getCoords(), 0, dest, pos, m.getCoords().length);
        pos += m.getCoords().length;
    }

    void copyBufferToCluster(double[] src, int srcpos, Cluster m) {
        int pos = srcpos;
        m.setNAtoms((int) src[pos]);
        pos++;
        //m.setTag(String.format("%d",src[pos])); pos++;
        m.setEnergy(src[pos]);
        pos++;
        int[] Nz = m.getAtomicNumber();
        for (int i = 0; i < Nz.length; i++) {
            Nz[i] = (int) src[pos];
            pos++;
        }

        System.arraycopy(src, pos, m.getCoords(), 0, m.getCoords().length);
        pos += m.getCoords().length;
    }

    @Override
    int iConverged() {
        throw new UnsupportedOperationException("Not yet implemented");
    }

    @Override
    public double[][] MasterInitialize() {
        //initialize randome generator for swapping
        ranSwap = new Random(seed);

        LoadTemperature();

        logger.fine(name + ": LoadStructure: start reading");
        LoadInitialStructure();

        for (int i = 0; i < ncpu; i++) {
            logger.fine(name + ": LoadStructure: process " + i + "th will get structure with tag=" + pool.get(i).getTag());
        }

        //arrange the initialized data into package
        double[][] answer = new double[ncpu][];
        for (int i = 0; i < ncpu; i++) {
            Cluster m = pool.get(i);
            answer[i] = new double[1 + 1 + 1 + 1 + m.getAtomicNumber().length + m.getCoords().length];

            int pos = 0;
            answer[i][pos] = T[i];
            copyClusterToBuffer(m, answer[i], 1);
        }

        SlaveInitialize(answer[MASTER]);

        return answer;
    }

    @Override
    public double[] SlaveInitialize(double[] data) {
        double[] answer = null;
        try {
            MTools.PrintArray(data);
            mc = new MonteCarloSimulator();
            mc.setName(getMyRank());
            mc.initialize(args, xmllog, new FileWriter(new File("struct-" + name)), seed + 2 * (myrank + 1), seed + 3 * (myrank + 1));
            //read the initialized data from package
            int pos = 0;
            mc.setTemperature(data[pos]);
            pos++;
            Cluster m = new Cluster();
            copyBufferToCluster(data, pos, m);
            mc.setCurrentEnsemble(m);
            mc.getCurrentEnsemble().Write(System.out, sFormatOut);
            //print out MonteCarlo info
            xmllog.writeNormalText(mc.XMLInfo(3));

        } catch (IOException ex) {
            Logger.getLogger(TaskMPIParallelTempering.class.getName()).log(Level.SEVERE, null, ex);
        }
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
    ArrayList<Cluster> initpop, pool;
    private double[] T; //!< array of temperatures
    Random ranSwap; //!< random number generator for swapping

    String getMyRank() {
        return String.format("c%2d", myrank);
    }

    void LoadTemperature() {
        xmllog.writeEntity("LoadTemperature");
        T = new double[ncpu];
        if (!sFileTemp.isEmpty()) { //if file for temperature exists, read it
            xmllog.writeAttribute("Input", sFileTemp);
            try {
                Scanner fileTemp = new Scanner(new File(sFileTemp));
                for (int i = 0; i < ncpu; i++) {
                    T[i] = fileTemp.nextDouble();
                }
                fileTemp.close();
            } catch (FileNotFoundException ex) {
                Logger.getLogger(TaskMPIParallelTempering.class.getName()).log(Level.SEVERE, null, ex);
            }
        } else {
            xmllog.writeAttribute("Input", "FromFormula");
            double alpha = Math.pow(Tmax / Tmin, 1.0 / (ncpu - 1));
            for (int i = 0; i < ncpu; i++) {
                T[i] = Tmin * Math.pow(alpha, i);
            }
        }



        for (int i = 0; i < ncpu; i++) {
            xmllog.writeEntity("CPU");
            xmllog.writeAttribute("Rank", Integer.toString(i)).writeAttribute("T", Double.toString(T[i]));
            xmllog.endEntity().flush();
        }
        xmllog.endEntity().flush();
    }

    void LoadInitialStructure() {
        try {
            //read initial structures from file
            xmllog.writeEntity("LoadStructure");
            xmllog.writeAttribute("InputStructure", sFileIn).flush();
            Scanner fileIn = new Scanner(new File(sFileIn));

            initpop = new ArrayList<Cluster>();

            Cluster tmpCluster = new Cluster();
            while (tmpCluster.Read(fileIn, sFormatIn)) {
                //tmpCluster.Write(System.out, sFormatOut);
                if (tmpCluster.getNAtoms() != 0) {
                    initpop.add((Cluster) tmpCluster.clone());
                } else {
                    break;
                }

            }
            fileIn.close();
            xmllog.writeAttribute("NumOfStruct", Integer.toString(initpop.size())).flush();

            //prepare random structures for slave
            pool = new ArrayList<Cluster>();
            for (int i = 0; i < ncpu; i++) {
                int iSelect = i % initpop.size();
                pool.add((Cluster) initpop.get(iSelect));
            }
            xmllog.endEntity().flush();
        } catch (FileNotFoundException ex) {
            logger.severe(ex.getMessage());
        }
    }

    void RunParalelTempering(boolean bExchange) {

        int nMaxPasses;

        if (bExchange) {
            nMaxPasses = nRuns;
            mc.setName("CYC1");
        } else {
            nMaxPasses = nWarmUps;
            mc.setName("CYC0");
        }


        int nMaxMoves = nMaxPasses * mc.nMCperPass();


        int nStructSampleRate =  (structSampleRate < 1)? (int) (structSampleRate * nMaxMoves): (int) structSampleRate ;

        int nDataSampleRate =  (dataSampleRate < 1)? (int) (dataSampleRate * nMaxMoves): (int) dataSampleRate ;
        
//        int realSwapLogRate = (int) pr.swapLogRate;
//        if (swapLogRate < 1) {
//            realSwapLogRate = (int) (swapLogRate * MaxPasses / nSwap);
//        }
//        if (realSwapLogRate == 0) {
//            realSwapLogRate = 1e2;
//        }

        int nTotalSampledData = nMaxMoves / nStructSampleRate;


        long duration = System.currentTimeMillis();

        mc.reset();


        xmllog.writeEntity(mc.getName());

        xmllog.writeAttribute("TotalMCmoves", ""+nMaxMoves);
        xmllog.writeAttribute("DataSampleRate", ""+nDataSampleRate);
        xmllog.writeAttribute("StructSampleRate", ""+nStructSampleRate);
        //xmllog.writeAttribute("SwapRate", ""+nStructSampleRate);
        xmllog.writeAttribute("TotalSampledData", ""+nTotalSampledData);

        xmllog.writeNormalText(mc.XMLInfo(5));

        //Data_str wcTemp,wcTemp2;

//        if ((bExchange)) { //master report to main streams
//            for (int i = 0; i < ncpu; i++) {
//                SwapR[i] = 0;
//            }
//            sTmp = "Step ErgodicConv";
//            for (int i = 0; i < ncpu; i++) {
//                sTmp += str(format(" %lf") % T[i]);
//            }
//            masterLog.SetFrom_(mc.Name_());
//            masterLog.Get(sTmp);
//        }

        int iMove = 0; //number of MC moves

        while (iMove < nMaxMoves) {
            for (int i = 0; i < (nSwaps); i++) { //running 1 interval
                iMove+=mc.performMCMoves(iMove, nSwaps*mc.nMCperPass(),nDataSampleRate,nStructSampleRate);
//                if (bExchange) {
//                    if ((mc.hist.size() >= TotalSampleData) || nPass >= MaxPasses - 1) {
//                        mc.CalcTherProp();
//                        mainLog.Get(str(format("[Step=%d] [nCvlog=%d] ") % nMove % nCvlog) + mc.Info(1), "Cvlog");
//                        mc.hist.WriteHist(mainLog.Stream_(), nCvlog);
//                        mc.hist.ThrowData(pr.nDataThrow);
//                        //mainLog.Get(str(format("[Step=%d] ")%nMove)+mc.Info(1),"CvLog");
//                        nCvlog++;
//                    }
//                }
           }

//                
            }

//            if (bExchange) {//start swapping
//                //calculate ergodicity and swapping rate
//                if ((nSwap % realSwapLogRate == 0) || (nPass == MaxPasses - 1)) { //if log event occurs, output swapping rate and ErgodicConv to main stream
//                    double meanE[
//                    MAX_ncpu]
//                    ;
//                     mc  .hist.CalcAll
//                    ();
//                    double tmpMean=mc.hist.mean;//cout<<myrank<<" +++ "<<mc.hist.mean<<endl;
//                    COMM_WORLD.Gather( & tmpMean, 1, DOUBLE, meanE, 1, DOUBLE, MASTER);
//
//                    if (myrank == MASTER) {
//                        #
//                        ifdef PTMC_DEBUG
//                        sTmp = " MeanE = ";
//                        for (int i = 0; i < ncpu; i++) {
//                            sTmp += str(format(" %lf") % meanE[i]);
//                        }
//                        masterLog.Get(sTmp);
//                        #
//                        endif
//
//
//                        double ErgodicConv = Ergodicity(meanE);
//
//                        sTmp = str(format("%d %lf") % nMove % ErgodicConv);
//                        for (int i = 0; i < ncpu; i++) {
//                            if (nSwap) {
//                                sTmp += str(format(" %1.4lf") % (SwapR[i] / (2.0 * nSwap)));
//                            } else {
//                                sTmp += " 0";
//                            }
//                        }
//                        masterLog.Get(sTmp);
//                    }
//
//                    if (access(fnDead.c_str(), F_OK) == 0) { //if dead file occurs, quitting
//                        mc.CalcTherProp();
//                        mainLog.Get(str(format("[Step=%d] [nCvlog=%d] ") % nMove % nCvlog) + mc.Info(1), "Cvlog");
//                        mc.hist.WriteHist(mainLog.Stream_(), nCvlog);
//                        mc.hist.ThrowData(pr.nDataThrow);
//                        //mainLog.Get(str(format("[Step=%d] ")%nMove)+mc.Info(1),"CvLog");
//                        nCvlog++;
//                        mainLog.Get("  quitting because of DEAD file ");
//                        break;
//                    }
//                }
//
//                //after running, passing ensemble to MASTER
//                buff.Get(mc.ens);
//                #
//                ifdef PTMC_DEBUG
//                mainLog.Get("passing ensemble to MASTER");
//                #
//                endif COMM_WORLD
//                .Gather( & buff, 1, newtype, pool, 1, newtype, MASTER);
//
//                //MASTER swapping ensembles for the next run
//                #
//                ifdef PTMC_DEBUG
//                mainLog.Get("got ensembles from slaves, now swapping them");
//                #
//                endif
//
//                switch (pr.iExchangeMode) {
//                    case EM_ADJACENT:
//                        if (myrank == MASTER) {
//                            SwapReplica_adjacent(nSwap);
//                        }
//                        break;
//                    case EM_ALLEXCHANGE:
//                        if (myrank == MASTER) {
//                            SwapReplica_allexchange(nSwap);
//                        }
//                        break;
//                    case EM_FEEDBACK:
//                        //if(myrank==MASTER)
//                        SwapReplica_feedback(nSwap);
//                        break;
//                    case EM_RANDOMWALK:
//                        if (myrank == MASTER) {
//                            SwapReplica_randomwalk(nSwap);
//                        }
//                        break;
//                }
//                nSwap++;
//
//                //MASTER pass ensemble back to slaves
//                #
//                ifdef PTMC_DEBUG
//                mainLog.Get("MASTER pass ensemble back to slaves");
//                #
//                endif COMM_WORLD
//                 .Scatter(pool,1, newtype, &  buff, 1, newtype, MASTER);
//                buff.Migrate(mc.ens);//get it and correct the order of O and H
//
//            }
            //if(access(fnDead.c_str(), F_OK)==0){ //if dead file occurs, quitting
            //	flog<<" +++ "<<myrank<<": quitting because of DEAD file "<<fnDead<<endl;
            //	break;
            //}
        //}


        duration=System.currentTimeMillis()-duration;

            xmllog.writeNormalText(mc.XMLInfo(4));

            xmllog.writeEntity("ConvergenceInfo");
                xmllog.writeAttribute("ElapsedTime", ""+duration/1000);
                xmllog.writeAttribute("TotalMCmoves", ""+iMove);
                xmllog.writeAttribute("Speed", ""+duration/1000/iMove);
            xmllog.endEntity();
        xmllog.endEntity();

    }
}

