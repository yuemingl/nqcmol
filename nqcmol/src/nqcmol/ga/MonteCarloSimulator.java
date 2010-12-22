/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nqcmol.ga;

import java.io.IOException;
import java.io.OutputStream;
import java.io.Writer;
import java.util.Random;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.cluster.Cluster;
import nqcmol.cluster.MolExtra;
import nqcmol.potential.Potential;
import nqcmol.tools.Histogram;
import nqcmol.tools.MTools;
import nqcmol.tools.XmlWriter;
import org.apache.commons.math.stat.descriptive.DescriptiveStatistics;
import org.kohsuke.args4j.*;

/**
 *
 * @author chinh
 */
public class MonteCarloSimulator {

    @Option(name = "-ff", usage = "Potential model. [LJ}", metaVar = "POTENTIAL")
    String sPotential = "LJ";
    @Option(name = "-a", usage = "Parameter file for potential if applicant.", metaVar = "FILE")
    String sFileParam = "";
    @Option(name = "-u", usage = "Unit of energy if applicant: Hartree, kcal/mol, kJ/mol, eV, au.", metaVar = "STRING")
    String sUnit = "";
    @Option(name = "-AtomdR", usage = "AtomdR", metaVar = "DOUBLE")
    double AtomdR = 0.2; //!< maximum displacement of atomic moves
    @Option(name = "-MoldR", usage = "AtomdR", metaVar = "DOUBLE")
    double MoldR = 0.1;	//!< maximum displacement of molecular moves
    @Option(name = "-MoldPhi", usage = "MoldPhi (in degree)", metaVar = "DOUBLE")
    double MoldPhi = 10; //!< maximum rotation of molecular moves
    double ScalingBase;
    double ScalingIndex;
    protected String[] args = null;
    protected CmdLineParser parser = null;
    double ConfinedRadius; //!< for confining
    Random ranMov = new Random(); //!< random number generator for moving
    Random ranBoltz = new Random(); //!< random number generator for Boltzmann factor
    double meanE = 0, Cv = 0;
    double AcceptR, RejectR, DisR; //overall accept rate
    double MolAcceptR, MolRejectR, MolDisR; //!< molecular accepted, rejected, disassociated rate
    double AtomAcceptR, AtomRejectR, AtomDisR; //!< atomic accepted, rejected, disassociated rate
    // protected static Logger logger=Logger.getLogger(MonteCarloSimulation.class.getName());
    double T = 0;
    XmlWriter xmlLog;
    Writer structLog;

    public double getTemperature() {
        return T;
    }

    public void setTemperature(double T) {
        beta = 1.0 / (T * getkB());
        this.T = T;
    }
    double beta = 0;
    double kB = 0;//!< Boltzmann constants

    double getkB() {//!< get Boltzmann constant
        return kB;
    }

    void initialize(String[] args, XmlWriter xmlLog, Writer structLog, long seedMov, long seedBoltz) {
        //initialize potential
        pot = MolExtra.SetupPotential(sPotential, sFileParam, sUnit, "DFPMIN");

        //calculate Boltzmann's constants
        String unit = pot.getUnit();
        if (unit.contentEquals("Hartree")) {
            kB = 3.16681536e-6;
        } else if (unit.contentEquals("kcal/mol")) {
            kB = 1.98587767e-4;
        } else if (unit.contentEquals("kJ/mol")) {
            kB = 3.16681536e-6;
        } else if (unit.contentEquals("eV")) {
            kB = 3.16681536e-6;
        } else {
            kB = 1.0;
        }


        //scale the parameter according to temperature
        scaleParam(ScalingBase, ScalingIndex);

        //initialize random generator
        ranMov.setSeed(seedMov);
        ranBoltz.setSeed(seedBoltz);

        //initialize logging stream
        this.xmlLog = xmlLog;
        this.structLog = structLog;
    }

    /**
     * \brief Scaling atomic and molecular moves for temperature independent jumping rate
     *  new parameter = \f$ (old parameter)\exp{\alpha(T/T_0)} \f$ */
    void scaleParam(double To, double alpha) {
        double hs = T / To;
        AtomdR *= Math.pow(hs, alpha);
        MoldR *= Math.pow(hs, alpha);
        MoldPhi *= Math.pow(hs, alpha);
    }
    String name = "MC";   //!< name of this MC instance

    public String getName() {
        return name;
    }

    public void setName(String name) {
        this.name = name;
    }

    void calcThermalProperties() {
        meanE = stats.getMean();
        Cv = beta * beta * stats.getVariance() / ens.getNAtoms();
    }

    double calcEnergy(Cluster p) { //calculate energy of individual
        pot.getEnergy(p);
        return p.getEnergy();
    }
    protected Cluster ens = new Cluster(); //!< the current ensemble

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
    public String XMLInfo(int verbose) {
        String info = "";
        try {
            //!< print setting parameters
            Writer writer = new java.io.StringWriter();
            XmlWriter xmlwriter = new XmlWriter(writer);
            xmlwriter.writeEntity("MCSimulatorInfo");
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

    int nMCperPass() { //!< number MC moves/ pass = atomic moves + molecular moves
        return ens.getNAtoms() + ens.getNonHydrogenNum();
    }

    void reset() {
        meanE = Cv = 0;
        AcceptR = RejectR = DisR = 0; //overall accept rate
        MolAcceptR = MolRejectR = MolDisR = 0; // molecular accept rate
        AtomAcceptR = AtomRejectR = AtomDisR = 0; // molecular accept rate
        stats.clear();
    }
    
    Histogram stats = new Histogram();
    Potential pot; //!< pointer to objective function

    /**
     * Perform a series of MC moves
     *
     */
    int performMCMoves(int iCounter, int nMoves, int nDataSampleRate, int nStructSampleRate) {
        int c = iCounter;
        for (int i = 0; i < nMoves; i++) {
            double oldE = ens.getEnergy();
            int iJump = performOneMCMove(c);
            switch(iJump){
//
//                AcceptR++;		MolAcceptR+=(iMovType==1);		AtomAcceptR+=(iMovType==0);		return 0;
//		RejectR++;		MolRejectR+=(iMovType==1);		AtomRejectR+=(iMovType==0);		return 1;
//                case 0: acceptRate++;
//                case 1: rejectRate++;
            }

            double newE = ens.getEnergy();

            if (nDataSampleRate > 0) {
                if (c % nDataSampleRate == 0) {
                    stats.addValue(newE);
                }
            }

            if (nStructSampleRate > 0) {
                if ((c % nStructSampleRate == 0)) {
                    calcThermalProperties();
                    xmlLog.writeEntity("Step").writeAttribute("id", "" + c).writeAttribute("OldE", "" + oldE).writeAttribute("NewE", "" + newE);
                    xmlLog.endEntity().flush();
                    ens.setTag("" + c);
                    ens.Write(structLog, "XYZ");
                }
            }

            c++;
        }
        return c;
    }

    /**
     * Perform 1 MC moves
     * @return 0 if successfully jumping, 1 if unsucessully, 2 if disassociated structure
     */
    int performOneMCMove(int iMove) {
        int answer = 0;
        double oldE = ens.getEnergy();
        Cluster trialEns = new Cluster();
        trialEns.set(ens);

        int iTag = iMove % nMCperPass();
        double[] ranVec = new double[3];
        if (iTag < ens.getNAtoms()) {

            MTools.generateRandomVector(ranVec, ranMov, AtomdR);
            ens.Translate_atom(iTag, ranVec, 1.0);
        } else {
            int[] mol = ens.getMolecule(iTag - iTag);

            double[] ranAxis = new double[3];
            MTools.generateRandomVector(ranAxis, ranMov, 1.0);
            ens.RotateAxis_mol(mol, Math.toRadians(MoldPhi * (2.0 * ranMov.nextDouble() - 1.0)), ranAxis);

            MTools.generateRandomVector(ranVec, ranMov, MoldR);
            ens.Translate_mol(mol, ranVec, 1.0);
        }

//        string error;
//        if(!inv.bCorrect(error)){
//            DisR++;
//            AtomDisR++;
//            return 2;
//        }

        double newE = calcEnergy(trialEns);
        //inv.Write(cout,F_XYZ);

        if (ranBoltz.nextDouble() < Math.exp(-beta * (newE - oldE))) {
            ens.set(trialEns);
        } else {
            answer = 1;
        }
        return answer;
    }
}
