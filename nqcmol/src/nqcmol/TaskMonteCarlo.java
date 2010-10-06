/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.ClusterOperation;
import nqcmol.cluster.Cluster;
import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.tools.MTools;
import org.apache.commons.math.stat.descriptive.SummaryStatistics;
import org.kohsuke.args4j.*;
import java.util.Random;
import java.util.Vector;



/**
 *
 * @author nqc
 */
public class TaskMonteCarlo extends TaskCalculate {
    @Option(name="-ftol",usage="Energy tolerance. [1e-9]",metaVar="DOUBLE")
    double ftol=1e-9;

    @Option(name="-gtol",usage="RMS of gradient tolerance.[1e-6]",metaVar="DOUBLE")
    double gtol=1e-6;

    @Option(name="-mstep",usage="Maximum step size of each iteration. [0.1]",metaVar="DOUBLE")
    double mstep=1e-1;

    @Option(name="-nOpts",usage="Maximum number of iterations in local optimization. [0]",metaVar="INTEGER")
	int nOpts=0;

	@Override
	public String getName(){
		return "MonteCarlo";
	}

    static final public String Option="mc";

    static final public String Descriptions="\t "+Option+" \t - "+ "Monte Carlo simulation\n";

    public static void main(String[] args) throws IOException  {
        new TaskMonteCarlo().Execute(args);
	}

	@Override
	protected void Process() {		
            LoadEnvParam();
            Work();	
	}


	//String fnCtl="cnt-pgaWater.ctl";//!< control file of environment parameters
	String fnDead="dead";
	String fnDump="dump";
   

    File fDump;
    File fDead;

    @Option(name="-n",usage="Index of input structure. [0]",metaVar="INTEGER")
    int iN=0;

    @Option(name="-nRuns",usage="Maximum step size of each iteration. [100]",metaVar="INTEGER")
    int nRuns=10;

    @Option(name="-nLogRate",usage="Logging rate. [0.01]",metaVar="DOUBLE")
    double nLogRate=0.01;

    @Option(name="-WallTime",usage="Walltime in seconds. 0 means not be used. [0]",metaVar="INTEGER")
    long WallTime=0;

    @Option(name="-T",usage="Temperature in Kevin. [300]",metaVar="DOUBLE")
    double T=300;

    @Option(name="-fnArchive",usage="File name for archival. [archive.xyz]",metaVar="STRING")
    String fnArchive="archive.xyz";

    @Option(name="-arvSize",usage="File name for archival. [10000]",metaVar="DOUBLE")
    int arvSize=10000;

    @Option(name="-arvSim",usage="Similarity limit used in duplicate checking in archival. [0.96]",metaVar="DOUBLE")
    double arvSim=0.96;

    @Option(name="-arvEnergyGap",usage="Energy gap used in duplicate checking in archival. [0.002]",metaVar="DOUBLE")
    double arvEnergyGap=0.002;

    @Option(name="-arvRmsGrad",usage="Largest RMS gradient to be accepted in archival. [0.0001]",metaVar="DOUBLE")
    double arvRmsGrad=0.0001;

    //int arvCheckFreq;


	int LoadEnvParam(){                 
        //seed=inp.getLong("seed");	if(pr.seed==0) pr.seed=time(0);
        //ranBoltz.SetSeed_(pr.seed);
        pot.setEnergyTol(ftol);
        pot.setGradientTol(gtol);
        pot.setMaxStepSize(mstep);
        pot.setMaxEvals(nOpts);

        int i=0;
        while(mol.Read(fileIn, sFormatIn)){
            if(i==iN) break;
            i++;
        }
        fileIn.close();
        if(pot.getUnit().contentEquals("Hartree")) kB=3.16682968e-6;//Hartree*K^-1
        if(pot.getUnit().contentEquals("kcal/mol")) kB=0.1986e-3;//kcal/mol*K^-1
        if(pot.getUnit().contentEquals("kJ/mol")) kB=8.309424e-4;//kJ/mol*K^-1

        kBT=kB*T;

        xmllog.writeEntity("Setting");
        xmllog.writeAttribute("nRuns",Integer.toString(nRuns));
        xmllog.writeAttribute("nLogRate",Double.toString(nLogRate));
        xmllog.writeAttribute("WallTime",Long.toString(WallTime));
        xmllog.writeAttribute("Temperature",Double.toString(T));
        xmllog.writeAttribute("kBT",Double.toString(kBT));
        xmllog.endEntity().flush();
		//string sTmp=str(format("[nRuns=%d] [WallTime=%d] [ETarget=%lf] [nLogRate=%1.3lf]") %pr.nRuns %pr.WallTime %pr.ETarget %pr.nLogRate);
		//sTmp+=str(format("[seed=%d]") %pr.seed);

		//mainLog.Get(sTmp);
		return 0;
	}


    SummaryStatistics rateAccept=new SummaryStatistics ();
    SummaryStatistics nAccept_arv=new SummaryStatistics ();
    SummaryStatistics nDuplicate_arv=new SummaryStatistics ();
    SummaryStatistics nIncorrect_arv=new SummaryStatistics ();
    SummaryStatistics nLargeRMS_arv=new SummaryStatistics ();
    SummaryStatistics nNegativeFreq_arv=new SummaryStatistics ();

    Random ranBoltz=new Random();

    ClusterOperation oper=new ClusterOperation();

    Vector<Cluster> arv=new Vector<Cluster>();

    double kB=1.0;
    double kBT=1.0;

     Cluster trialMol=new Cluster();
//    int iConverged(){ //!< check whether system converges
//		if (iStep>=pr.nRuns)  return 1;
//		if(access(fnDead.c_str(), F_OK)==0) return 3;
//		if((time(0)-timer)> pr.WallTime)	return 4;
//		return 0;
//	}

    int iConverged(int iStep,long initialTime){ //!< check whether system converges
		if (iStep>=nRuns)  return 1;
		if (fDead.exists()) return 3;
		if(WallTime>0)
            if((System.currentTimeMillis()-initialTime)/1000000> WallTime)	return 4;
		return 0;
	}

	void Logging(){       
            xmllog.writeEntity("MCStep");
            xmllog.writeAttribute("id", Long.toString(rateAccept.getN()));
            xmllog.writeAttribute("CurrentE", Double.toString(mol.getEnergy()));
            xmllog.writeAttribute("TrialE", Double.toString(trialMol.getEnergy()));
            xmllog.writeAttribute("AcceptR", String.format("%1.2f",rateAccept.getMean()));
            //xmllog.writeAttribute("RefuseR", Double.toString(1-rateAccept.getMean()));
            xmllog.writeAttribute("nArchive", Integer.toString(arv.size()));
            xmllog.writeAttribute("AcceptR_arv",String.format("%1.2f",nAccept_arv.getMean()));
            xmllog.writeAttribute("DupR_arv", String.format("%1.2f",nDuplicate_arv.getMean()));
            xmllog.writeAttribute("IncorrectR", String.format("%1.2f",nIncorrect_arv.getMean()));
            xmllog.writeAttribute("LargeR", String.format("%1.2f",nLargeRMS_arv.getMean()));
            //xmllog.writeAttribute("NegFreq", String.format("%1.2f",nNegativeFreq_arv.getMean()));
            xmllog.endEntity().flush();        
	}

    void IncreaseCounter(int d){        
            nDuplicate_arv.addValue(MTools.booleanToInt(d==-1));
            nIncorrect_arv.addValue(MTools.booleanToInt(d==-2));
            nLargeRMS_arv.addValue(MTools.booleanToInt(d==-3));
            nNegativeFreq_arv.addValue(MTools.booleanToInt(d==-4));
            nAccept_arv.addValue(MTools.booleanToInt(d==0));
    }
    
    int UpdateArchive(Cluster m){
        Cluster p=(Cluster)m.clone();
        p.CalcUSRsignature();

        if(p.isNAN())  return -2;//incorrect structure
		if(p.getRmsGrad()>arvRmsGrad){ return -3;}

        //	if(arvCheckFreq){
        //		double f=bCheckFreq(tmp);	//cout<<"frequency = "<<f<<endl;
        //		if(f<0) return -4;
        //	}

        for(int i=0;i<arv.size();i++){
            double sim=arv.get(i).CalcUSRSimilarity(p);
            if((sim>=arvSim)&&(Math.abs(arv.get(i).getEnergy() - p.getEnergy()) < arvEnergyGap))
                return -1; //duplicate
        }

        if(arv.size()<arvSize){
            for(int i=0;i<arv.size();i++)
                if(arv.get(i).getEnergy() >p.getEnergy()){
                    arv.insertElementAt(p, i);
                    return 0;
                }
            arv.add(p);
        }

	return 0;
}

	void Work(){       
            //Individual tmp;
            fDump = new File(fnDump);   if (fDump.exists()) {       fDump.delete(); }
            fDead = new File(fnDead);   if (fDead.exists()) {       fDead.delete(); }

            int realLogRate = (int) nLogRate;
            if (nLogRate < 1) realLogRate = (int) (nLogRate*nRuns);
            if (realLogRate == 0) realLogRate = 1;

            xmllog.writeEntity("MC");
            xmllog.writeAttribute("ExpectedMCmoves", Integer.toString(nRuns));
            xmllog.writeAttribute("realLogRate", Integer.toString(realLogRate));
            //Logging();

            int iStep=0;
            int iConv = 0; //iStep=1;
            long duration=System.currentTimeMillis();

            mol.setFreqs(null);
           

            while (true) {                               
                trialMol.Get(mol);
                if (iStep!=0) {
                    oper.RandomMolMove(trialMol);
                }

                //trialMol.Write(System.out, "xyz");

                //pot.setCluster(trialMol);
				if(nOpts!=0) pot.Optimize(trialMol);
                else pot.getEnergy(trialMol);

				double deltaE = trialMol.getEnergy()-mol.getEnergy();

                //System.out.printf("Step = %d  newE=%f  oldE=%f detalE=%f rate=%f\n",iStep, trialMol.getEnergy(),mol.getEnergy(),deltaE,Math.exp(deltaE/kBT));
                //cout<<" dsdf 12 "<<ga.bf<<" "<<ga.bf->Dim_()<<endl;
                if (ranBoltz.nextDouble() < Math.exp(-deltaE/kBT)) {
                    mol.Get(trialMol);
                    mol.Center();
                    //mol.Write(System.out, "xyz");
                    rateAccept.addValue(1);
                } else {
                    rateAccept.addValue(0);
                }
                
                int d = UpdateArchive(trialMol);
                IncreaseCounter(d);           

                 if (iStep % realLogRate == 0) {
                    Logging();
                    if(fDump.exists()) {
                        WriteArchive();
                        fDump.delete();
                    }
                    iConv = iConverged(iStep,duration);
                    if (iConv != 0) {
                        break;
                    }
                }
                
                iStep++;
            }

            duration=(System.currentTimeMillis()-duration);
            String reason="Undetermined";
            switch (iConv) {
                case 1:
                    reason="Exceeding max runs";
                break;
                case 2:
                    reason="Reach the target";
                    break;
                case 3:
                    reason="Found dead file";
                    break;
                case 4:
                    reason="Reach time limitation";
                break;
            }

            xmllog.writeEntity("ConvergenceInfo");
            xmllog.writeAttribute("TotalMCmoves", Integer.toString(iStep));
            xmllog.writeAttribute("Duration", Double.toString(duration/1000.0));
            xmllog.writeAttribute("Speed", String.format("%1.3f",duration/1000.0/iStep));
            xmllog.writeAttribute("Reason", reason);
            xmllog.endEntity();
            xmllog.endEntity().flush();

            WriteArchive();        
	}

    private void WriteArchive() {
        if(arv.size()>0){
            try {
                FileWriter fArchive = new FileWriter(new File(fnArchive));
                for(int i=0;i<arv.size();i++){
                    arv.get(i).setTag(Integer.toString(i));
                    arv.get(i).Write(fArchive, "xyz");
                }
                fArchive.close();
            } catch (IOException ex) {
                Logger.getLogger(TaskMonteCarlo.class.getName()).log(Level.SEVERE, null, ex);
            }
        }
    }
}
