/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

//import com.thoughtworks.xstream.XStream;
import java.io.*;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.*;


import nqcmol.potential.*;
import nqcmol.tools.MTools;
import nqcmol.tools.XmlWriter;
//import org.xml.sax.SAXException;
//import org.xml.sax.helpers.AttributesImpl;

/**
 *
 * @author nqc
 */
public class MCalculate {

	@Option(name="-ff",usage="Specify the potential",metaVar="POTENTIAL")
    String sPotential="LJ";

	@Option(name="-i",usage="input file name",metaVar="FILE")
    String sFileIn="test/LJlm/lj-pool.xyz";

	@Option(name="-o",usage="output file name",metaVar="FILE")
    String sFileOut="";

	@Option(name="-nRuns",usage=" number of interations",metaVar="INTNUM")
    int nRuns=1;

	@Option(name="-nScale",usage=" number of interations",metaVar="INTNUM")
    int nScale=0;

	@Option(name="-grad",usage="Benchmark gradients or not.")
    boolean isGrad= false;

    @Option(name="-h",usage="Print out the help")
    boolean isHelp= false;
	
	CmdLineParser parser=null;

	public void ParseArguments(String[] args){
		parser = new CmdLineParser(this);
		try {
			parser.parseArgument(args);
		} catch (CmdLineException ex) {
			Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		}		
	}


	public void CalculateEnergy(String[] args)  {
		try {
			ParseArguments(args);
			if (isHelp) {
				parser.printUsage(System.out);
				return;
			}

			Potential pot = MolExtra.SetupPotential(sPotential);

			//System.out.println(" Come here");
			BufferedWriter writer=new BufferedWriter(new OutputStreamWriter(System.out));
			XmlWriter xmlwriter=new XmlWriter(writer);
			xmlwriter.writeEntity("BenchmarkEnergy");
			xmlwriter.writeEntity("Note");
			xmlwriter.writeText(" Time is measured in milliseconds");
			xmlwriter.endEntity();
			xmlwriter.writeNormalText(pot.Info(1));

			
			Cluster mol = new Cluster();
			Scanner scanner = new Scanner(new File(sFileIn));

			//System.out.println("What the fuck!"+sFileIn);
			int i = 0;
			//JSONArray jsChild2=new JSONArray();
			while (mol.Read(scanner, "xyz")) {
				//mol.Write(System.out,"xyz");

				int myruns = (int) ((nScale>0)?nRuns*Math.pow((double)nScale/mol.getNAtoms(),2):nRuns);
				myruns=Math.max(myruns, 1);
				double energy = 0;
				long duration = System.currentTimeMillis();
				pot.setCluster(mol);
				for (int k = 0; k < myruns; k++) {
					energy = pot.getEnergy(true);
				}
				duration = System.currentTimeMillis() - duration;
				double EnergySpeed=0;
				if(myruns!=0) EnergySpeed=duration/myruns;


				
				xmlwriter.writeEntity("Bench");
				xmlwriter.writeAttribute("Tag", Integer.toString(i));
				xmlwriter.writeAttribute("nAtoms", Integer.toString(pot.getCluster().getNAtoms()));
				xmlwriter.writeAttribute("nRuns", Integer.toString(myruns));
				xmlwriter.writeAttribute("Energy", Double.toString(energy));
				xmlwriter.writeAttribute("EnergyDuration", Long.toString(duration));
				xmlwriter.writeAttribute("EnergySpeed", Double.toString(EnergySpeed));


				if(isGrad){
					duration = System.currentTimeMillis();
					for (int k = 0; k < myruns; k++) {
						double[] grad = pot.getGradient(true);
					}
					duration = System.currentTimeMillis() - duration;

					double GradientSpeed=0;
					if(myruns!=0) GradientSpeed=duration/myruns;
					xmlwriter.writeAttribute("GradientDuration", Long.toString(duration));
					xmlwriter.writeAttribute("GradientSpeed", Double.toString(GradientSpeed));
				}


				xmlwriter.endEntity();
//				JSONObject jsBench = new JSONObject();
//				jsBench.put("Tag", i);
//				jsBench.put("nAtoms", pot.getCluster().getNAtoms());
//				jsBench.put("nRuns", myruns);
//				jsBench.put("Energy", energy);
//				jsBench.put("Duration", duration);
//				jsBench.put("Speed", (double) duration / myruns);
//				jsChild1.put("bench", jsBench);
//		

//				atts.clear();
//				atts.addAttribute("", "", "Tag", "", Integer.toString(i));
//				atts.addAttribute("", "", "nRuns", "", Integer.toString(myruns));
//				atts.addAttribute("","", "Energy", "", Double.toString(energy));
//				atts.addAttribute("","", "Duration", "", Long.toString(duration));
//				atts.addAttribute("","", "Speed", "", Double.toString((double)duration/myruns));
//
//
//
//				hd.startElement("", "", "bench", atts);
//				hd.endElement("", "", "bench");
//
				i++;
				
				//break;
			}
//			hd.endElement("", "", "nqc_ener");
//			hd.endDocument();

			xmlwriter.endEntity();
			xmlwriter.close();
			writer.close();
			
			//System.out.print(jsParent);
		} catch (IOException ex) {
			Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		} 
	}

	public void Optimize(String[] args){
		try {
			ParseArguments(args);
			if (isHelp) {
				parser.printUsage(System.out);
				return;
			}

			BufferedWriter writer=new BufferedWriter(new OutputStreamWriter(System.out));
			XmlWriter xmlwriter=new XmlWriter(writer);
			xmlwriter.writeEntity("Optimization");

			Potential pot = MolExtra.SetupPotential(sPotential);

			xmlwriter.writeNormalText(pot.Info(3));

			Cluster mol = new Cluster();
			Scanner scanner = new Scanner(new File(sFileIn));
			//System.out.println("What the fuck!"+sFileIn);
			int i = 0;
			while (mol.Read(scanner, "xyz")) {
				//mol.Write(System.out,"xyz");


				xmlwriter.writeEntity("Opt").writeAttribute("Tag",Integer.toString(i));
				
//				xmlwriter.writeEntity("XYZ");

//				xmlwriter.endEntity();
				xmlwriter.writeAttribute("nAtom",Integer.toString(12));
				pot.Optimize(mol);
				mol.Write(System.err, "xyz");
				xmlwriter.endEntity();
				
				i++;
			}
			xmlwriter.endEntity();
			writer.close();
		} catch (IOException ex) {
			Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		} 
	}

	public void ValidateGradients(String[] args)  {
		try {
			ParseArguments(args);
			if (isHelp) {
				parser.printUsage(System.out);
				return;
			}

			Potential pot = MolExtra.SetupPotential(sPotential);

			//reading cluster  from file
			Cluster mol = new Cluster();
			Scanner scanner = new Scanner(new File(sFileIn));
			int i = 0;
			while (mol.Read(scanner, "xyz")) {
				//mol.Write(System.out,"xyz");
				pot.setCluster(mol);
				System.out.print(pot.ValidateGradient(1.0E-4));
				
				i++;
			}
		} catch (FileNotFoundException ex) {
			Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	public void HarmonicVibrationAnalysis(String[] args){
		try {
			ParseArguments(args);
			if (isHelp) {
				parser.printUsage(System.out);
				return;
			}

			BufferedWriter writer=new BufferedWriter(new OutputStreamWriter(System.out));
			XmlWriter xmlwriter=new XmlWriter(writer);
			xmlwriter.writeEntity("VibrationAnalysis");


			Potential pot = MolExtra.SetupPotential(sPotential);

			xmlwriter.writeNormalText(pot.Info(1));

			Cluster mol = new Cluster();
			Scanner scanner = new Scanner(new File(sFileIn));
			int i = 0;
			while (mol.Read(scanner, "xyz")) {
				//mol.Write(System.out,"xyz");

				//calculate Hessian
//				pot.getGradient(mol);
//				System.err.printf("Gradient \n");
//				MTools.PrintArray(mol.getGradient());
				
				pot.getNumericalHessian(mol,1e-5);
//				System.err.printf("Hessian \n");
//				MTools.PrintArray(mol.getHessian());
				
				HarmonicVibration vib=new HarmonicVibration(mol);
				vib.CalcFreq();//calculate frequencies
				mol.Write(System.err, "xyz");
			


				//write a report
				xmlwriter.writeEntity("Vib").writeAttribute("Tag",Integer.toString(i));
//				xmlwriter.writeEntity("XYZ");
//				xmlwriter.endEntity();
				xmlwriter.writeAttribute("nAtom",Integer.toString(12));
				xmlwriter.endEntity();
				i++;
			}
			xmlwriter.endEntity();
			writer.close();
		} catch (IOException ex) {
			Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

}
