/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.ga;

import java.io.IOException;
import java.io.Writer;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.cluster.Cluster;
import nqcmol.cluster.MolExtra;
import nqcmol.potential.Potential;
import nqcmol.tools.XmlWriter;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.kohsuke.args4j.*;

/**
 *
 * @author chinh
 */
public class MonteCarloSimulator {

    @Option(name="-ff",usage="Potential model. [LJ}",metaVar="POTENTIAL")
    String sPotential="LJ";

    @Option(name="-a",usage="Parameter file for potential if applicant.",metaVar="FILE")
    String sFileParam="";

    @Option(name="-u",usage="Unit of energy if applicant: Hartree, kcal/mol, kJ/mol, eV, au.",metaVar="STRING")
    String sUnit="";

    @Option(name="-AtomdR",usage="AtomdR",metaVar="DOUBLE")
    double AtomdR=0.2; //!< maximum displacement of atomic moves

    @Option(name="-MoldR",usage="AtomdR",metaVar="DOUBLE")
    double MoldR=0.1;	//!< maximum displacement of molecular moves

    @Option(name="-MoldPhi",usage="MoldPhi (in degree)",metaVar="DOUBLE")
    double MoldPhi=10; //!< maximum rotation of molecular moves


    double ScalingBase;
    double ScalingIndex;


    protected String[] args=null;

    protected CmdLineParser parser = null;

    double ConfinedRadius; //!< for confining

    Random ranMov=new Random(); //!< random number generator for moving
    Random ranBoltz=new Random(); //!< random number generator for Boltzmann factor


    double meanE=0,Cv=0;
    double AcceptR,RejectR,DisR; //overall accept rate

    double MolAcceptR,MolRejectR,MolDisR; //!< molecular accepted, rejected, disassociated rate
    double AtomAcceptR,AtomRejectR,AtomDisR; //!< atomic accepted, rejected, disassociated rate

   // protected static Logger logger=Logger.getLogger(MonteCarloSimulation.class.getName());

    double T=0;

    public double getTemperature() {
        return T;
    }

    public void setTemperature(double T) {
        beta=1.0/(T*getkB());
        this.T = T;
    }

    double beta=0;

    double kB=0;//!< Boltzmann constants

    double getkB(){//!< get Boltzmann constant
        return kB;
    }

    void Initialize(String[] args,long seedMov,long seedBoltz){//!< initilize based on control file
        pot=MolExtra.SetupPotential(sPotential, sFileParam, sUnit, "DFPMIN");
        
        String unit=pot.getUnit();
        if(unit.contentEquals("Hartree")) kB=3.16681536e-6;
        else if(unit.contentEquals("kcal/mol")) kB=1.98587767e-4;
        else if(unit.contentEquals("kJ/mol")) kB=3.16681536e-6;
        else if(unit.contentEquals("eV")) kB=3.16681536e-6;
        else kB=1.0;


        ScaleParam(ScalingBase,ScalingIndex);//tune the parameter according to temperature
        ranMov.setSeed(seedMov);
        ranBoltz.setSeed(seedBoltz);
    }

    /*!
     *  \brief Perform Monte Carlo simulations
    *  \return Number of MC moves (NOT pass) = MC pass * number of MC moves/pass
    */
//    int Run(int iMove, //!< counter for MC moves
//              int MaxMove, //!< maximum number of MC moves
//              Logger log, //!< stream for logging                                                                                                                            int nLogRate=0 //!< logging rate
//    ){
//        double oldE,newE;
//        for(int i=0;i<MaxMoves;i++){
//            TryMCMove(oldE,newE,iMove);	hist.Push(ens.energy);
//
//            if(nLogRate) if((iMove%nLogRate==0)||(iMove==MaxMoves-1)){
//                CalcTherProp();
//                log.Get(str(format("[Step=%d] ")%iMove)+Info(1),"log");
//                ens.tag=iMove;
//                ens.Write(log.Stream_(),F_XYZ);
//            }
//            iMove++;
//        }
//        return iMove;
//    }


    /**
    * \brief Scaling atomic and molecular moves for temperature independent jumping rate
    *  new parameter = \f$ (old parameter)\exp{\alpha(T/T_0)} \f$ */
    void ScaleParam(double To,double alpha){
        double hs=T/To;
        AtomdR*=Math.pow(hs,alpha);
        MoldR*=Math.pow(hs,alpha);
        MoldPhi*=Math.pow(hs,alpha);
    }

    String name="MC";   //!< name of this MC instance

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    void CalcTherProp(){
        meanE=hist.getMean();
        Cv=beta*beta*hist.getVariance()/ens.getNAtoms();
    }

    double CalcEnergy(Cluster p){ //calculate energy of individual
        pot.getEnergy(p);
        return p.getEnergy();
    }


	protected Cluster ens=new Cluster(); //!< the current ensemble

    public Cluster getCurrentEnsemble() {
        return ens;
    }

    public void setCurrentEnsemble(Cluster ens) {
        this.ens = ens;
    }

    //=======================
    /**
	 * Return info of this class in xml format
	 * @param verbose =1: equation name and unit, 2= setting info, 3= status
	 */
    public String XMLInfo(int verbose){
        String info = "";
        try {
            //!< print setting parameters
            Writer writer = new java.io.StringWriter();
            XmlWriter xmlwriter = new XmlWriter(writer);
            xmlwriter.writeEntity("MonteCarloSimulator");
            if ((verbose & 1) != 0) {
                // current thermodynamic properties
            }
            
            if ((verbose & 2) != 0) {
                //setting parameter
                xmlwriter.writeAttribute("Temperature", Double.toString(T));
                xmlwriter.writeAttribute("kB", Double.toString(getkB()));
                xmlwriter.writeAttribute("beta", Double.toString(beta));
                xmlwriter.writeAttribute("ScalingBase", Double.toString(ScalingBase));
                xmlwriter.writeAttribute("ScalingIndex", Double.toString(ScalingIndex));
                xmlwriter.writeAttribute("AtomdR", Double.toString(AtomdR));
                xmlwriter.writeAttribute("MoldR", Double.toString(MoldR));
                xmlwriter.writeAttribute("MoldPhi", Double.toString(MoldPhi));
                xmlwriter.writeNormalText(pot.XMLInfo(1));
            }

//             if ((verbose & 4) != 0) {
//                //setting parameter
//                xmlwriter.writeNormalText(pot.XMLInfo(1));
//                xmlwriter.writeAttribute("Temperature", Double.toString(T));
//                xmlwriter.writeAttribute("kB", Double.toString(kB));
//                xmlwriter.writeAttribute("ScalingBase", Double.toString(kBT));
//                xmlwriter.writeAttribute("ScalingIndex_A", Double.toString(kBT));
//                xmlwriter.writeAttribute("ScalingIndex_M", Double.toString(kBT));
//                xmlwriter.writeAttribute("beta", Double.toString(kBT));
//            }


            xmlwriter.endEntity().close();
            info = writer.toString();
            
            //	info+="[name="+name+"] ";
            //	if(level&1){ // current thermodynamic properties
            //		info+=str(format("[E=%1.6f] [Cv=%1.6lf] " ) % ens.GetEnergy()  %Cv);
            //		info+=hist.Info(1);
            //	}
            //
            //	if(level&2){		//setting parameter
            //		info+=str(format("[T=%1.2lf] [kB=%1.15lf] [beta=%1.15lf] ")% T %kB %beta);
            //		info+=str(format("[ScalingBase=%lf] [ScalingIndex_A=%1.2lf] [ScalingIndex_M=%1.2lf] ") %ScalingBase %ScalingIndex_A%ScalingIndex_M);
            //		info+=str(format("[ConfinedRadius=%1.2lf] ") %ConfinedRadius);
            //		info+=str(format("[bMolMove=%d] [bAtomMov=%d] ")% bMolMove %bAtomMove);
            //		info+=hist.Info(2);
            //		//info+=str(format(" [SampleRate=%lf]") %SampleRate);
            //		info+=oper.Info();
            //
            //	}
            //
            //	if(level&4){ // statistical properties
            //		info+=str(format("[meanDeltaE=%1.6lf] [meanRate=%1.6lf] ")%meanDeltaE%exp(-meanDeltaE*beta));
            //		info+=str(format("[allMovR=%4.2lf %4.2lf %4.2lf] ")%(AcceptR)%(RejectR)%(DisR));
            //		if(bAtomMove)
            //		info+=str(format("[AtomMovR=%4.2lf %4.2lf %4.2lf] ")%(AtomAcceptR)%(AtomRejectR)%(AtomDisR));
            //		if(bMolMove)
            //		info+=str(format("[MolMovR=%4.2lf %4.2lf %4.2lf] ")%(MolAcceptR)%(MolRejectR)%(MolDisR));
            //	}
        } catch (IOException ex) {
            Logger.getLogger(MonteCarloSimulator.class.getName()).log(Level.SEVERE, null, ex);
        }
        return info;
}


    int nMCperPass(){ return ens.getNAtoms();} //!< number MC moves/ pass = atomic moves + molecular moves

    void Reset(){
        meanE=Cv=0;
        AcceptR=RejectR=DisR=0; //overall accept rate
        MolAcceptR=MolRejectR=MolDisR=0; // molecular accept rate
        AtomAcceptR=AtomRejectR=AtomDisR=0; // molecular accept rate
        hist.clear();
    }

    DescriptiveStatistics hist=new DescriptiveStatistics();

    Potential pot; //!< pointer to objective function

    /**
     * \brief Function to perform 1 MC moves
     * \return 0 if successfully jumping, 1 if unsucessully, 2 if disassociated structure
     */
//    int TryMCMove(double oldE, //!< energy before moving
//                          double newE, //!< energy after moving                                                                                                                                                                                                                              int iTag //!< index to indicate which kind of moves
//                        ){
//        int i;
//
//        oldE=ens.energy;
//        Cluster trialEns=new Cluster();
//        trialEns.set(ens);
//
//        if(iTag%ens.getNAtoms()==0){
//            double[] vec=new double[3];
//            MTools.generateRandomVector(vec, ranMov, AtomdR);
//        }
//
////
////        string error;
////        if(!inv.bCorrect(error)){
////            DisR++;
////            AtomDisR++;
////            return 2;
////        }
//        newE=CalcEnergy(trialEns);
//        //inv.Write(cout,F_XYZ);
//
//        if(ranBoltz.nextDouble() < Math.exp(-(newE-oldE)/(kBT))){
//            ens.set(trialEns);
//            AcceptR++;
//            AtomAcceptR++;
//            return 0;
//        }else{
//            RejectR++;
//            AtomRejectR++;
//            return 1;
//        }
//    }
//
//
//    void MC_AtomMove(Cluster p, int i){
//
//    }
}