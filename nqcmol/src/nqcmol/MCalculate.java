/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import com.generationjava.io.WritingException;
import com.generationjava.io.xml.XmlWriter;
import com.thoughtworks.xstream.XStream;
import java.io.*;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import javax.xml.transform.OutputKeys;
import javax.xml.transform.Transformer;
import javax.xml.transform.TransformerConfigurationException;
import javax.xml.transform.sax.SAXTransformerFactory;
import javax.xml.transform.sax.TransformerHandler;
import javax.xml.transform.stream.StreamResult;
import org.kohsuke.args4j.*;


import nqcmol.potential.*;
import org.xml.sax.SAXException;
import org.xml.sax.helpers.AttributesImpl;

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


	public void ener(String[] args)  {		
		try {
			ParseArguments(args);
			if (isHelp) {
				parser.printUsage(System.out);
				return;
			}
//			StreamResult streamResult=new StreamResult(System.out);
//			SAXTransformerFactory tf = (SAXTransformerFactory) SAXTransformerFactory.newInstance();
//			// SAX2.0 ContentHandler.
//			TransformerHandler hd = tf.newTransformerHandler();
//			Transformer serializer = hd.getTransformer();
//			serializer.setOutputProperty(OutputKeys.ENCODING, "ISO-8859-1");
//			serializer.setOutputProperty(OutputKeys.INDENT, "yes");
//
//			hd.setResult(streamResult);
//			hd.startDocument();
//			AttributesImpl atts = new AttributesImpl();
//			atts.clear();

			Potential pot = MolExtra.SetupPotential(sPotential);

			//System.out.println(" Come here");
			BufferedWriter writer=new BufferedWriter(new OutputStreamWriter(System.out));
			XmlWriter xmlwriter=new XmlWriter(writer);
			xmlwriter.writeEntity("BenchmarkEnergy");
			xmlwriter.writeAttribute("Potential",pot.getEquation());

			xmlwriter.writeText(pot.Info(1));

//			atts.addAttribute("", "", "Potential", "", pot.getEquation());
//			hd.startElement("", "", "nqc_ener", atts);
//			hd.startElement("","","Note",null);
//			String sTmp="Duration is measured in miliseconds. Speed is the time for one evaluation.";
//			hd.characters(sTmp.toCharArray(),0,sTmp.length());
//			hd.endElement("", "", "Note");

			//XStream xs=new XStream();
			//xs.omitField(this, parser);
			//xs.toXML(this,System.out);

//			JSONObject jsChild1 = new JSONObject();
//			jsChild1.put("Potential", pot.getEquation());
//			jsChild1.put("Note", "Duration is measured in miliseconds. Speed is the time for one evaluation.");
//
//			JSONObject jsParent = new JSONObject();
//			jsParent.put("nqc_ener", jsChild1);
			
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
				double speed=0;
				if(myruns!=0) speed=duration/myruns;

				xmlwriter.writeEntity("Bench");
				xmlwriter.writeAttribute("Tag", Integer.toString(i));
				xmlwriter.writeAttribute("nAtoms", Integer.toString(pot.getCluster().getNAtoms()));
				xmlwriter.writeAttribute("nRuns", Integer.toString(myruns));
				xmlwriter.writeAttribute("Energy", Double.toString(energy));
				xmlwriter.writeAttribute("Duration", Long.toString(duration));
				xmlwriter.writeAttribute("Speed", Double.toString(speed));
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
		} catch (WritingException ex) {
			Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		} 
	}

	public void opt(String[] args){
		
	}

	public void validgrad(String[] args)  {
		try {
			ParseArguments(args);
			if (isHelp) {
				parser.printUsage(System.out);
				return;
			}

			Potential pot = MolExtra.SetupPotential(sPotential);

			Cluster mol = new Cluster();
			Scanner scanner = new Scanner(new File(sFileIn));
			//System.out.println("What the fuck!"+sFileIn);
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

}
