/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.cluster.Cluster;

import org.kohsuke.args4j.*;



/**
 *
 * @author nqc
 */
public class TaskCalculateEnergy extends TaskCalculate {
	@Option(name="-nRuns",usage=" number of interations",metaVar="INTEGER")
    int nRuns=1;

	@Option(name="-nScale",usage=" number of interations",metaVar="INTEGER")
    int nScale=0;

	@Option(name="-grad",usage="Benchmark gradients or not.")
    boolean isGrad= false;

	@Override
	public String getName(){
		return "CalculateEnergy";
	}

    static final public String Option="ener";

    static final public String Descriptions="\t "+Option+" \t - "+ "calculate energy\n";

    public static void main(String[] args) throws IOException  {
        new TaskCalculateEnergy().Execute(args);
	}

	@Override
	protected void Process() {
        try {
            xmllog.writeEntity("Note");
            xmllog.writeText(" Time is measured in seconds");
            xmllog.endEntity();
            int i = 0;

            Cluster mol=new Cluster();
            Scanner fileIn = new Scanner(new File(sFileIn));
            while (mol.Read(fileIn, sFormatIn)) {
                int myruns = (int) ((nScale>0)?nRuns*Math.pow((double)nScale/mol.getNAtoms(),2):nRuns);
                myruns = Math.max(myruns, 1);
                double energy = 0;
                long startTime = System.nanoTime();
                pot.setCluster(mol);
                for (int k = 0; k < myruns; k++) {
                    energy = pot.getEnergy(true);
                }
                double duration = (System.nanoTime() - startTime) * 1e-9;
                double EnergySpeed = 0;
                if (myruns != 0) {
                    EnergySpeed = duration / myruns;
                }
                xmllog.writeEntity("Cluster");
                xmllog.writeAttribute("id", Integer.toString(i));
                xmllog.writeAttribute("Tag", mol.getTag());
                xmllog.writeAttribute("nAtoms", Integer.toString(pot.getCluster().getNAtoms()));
                xmllog.writeAttribute("nRuns", Integer.toString(myruns));
                xmllog.writeAttribute("Energy", Double.toString(energy));
                xmllog.writeAttribute("EnergyDuration", Double.toString(duration));
                xmllog.writeAttribute("EnergySpeed", Double.toString(EnergySpeed));
                if (isGrad) {
                    startTime = System.currentTimeMillis();
                    for (int k = 0; k < myruns; k++) {
                        double[] grad = pot.getGradient(true);
                    }
                    duration = (System.nanoTime() - startTime) * 1e-9;
                    double GradientSpeed = 0;
                    if (myruns != 0) {
                        GradientSpeed = duration / myruns;
                    }
                    xmllog.writeAttribute("GradientDuration", Double.toString(duration));
                    xmllog.writeAttribute("GradientSpeed", Double.toString(GradientSpeed));
                }
                xmllog.endEntity();
                xmllog.flush();
                i++;
            }
            fileIn.close();
        } catch (FileNotFoundException ex) {
            Logger.getLogger(TaskCalculateEnergy.class.getName()).log(Level.SEVERE, null, ex);
        }
	}

	
}
