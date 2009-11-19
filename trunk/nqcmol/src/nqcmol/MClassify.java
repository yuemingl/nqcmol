/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.HashMap;
import java.util.Scanner;
import java.util.logging.*;
import nqcmol.tools.MTools;
import nqcmol.tools.XmlWriter;
import org.kohsuke.args4j.*;

/**
 *
 * @author nqc
 */
public class MClassify {
	@Option(name="-i",usage="input file name",metaVar="FILE")
    String sFileIn="test/LJlm/lj-pool.xyz";

	@Option(name="-o",usage="output file name",metaVar="FILE")
    String sFileOut="";

	@Option(name="-pattern",usage="pattern of morphology to filter",metaVar="String")
    String sPattern="";
	
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

	public void ClassifyClusters(String[] args){
		try {
			ParseArguments(args);
			if (isHelp) {
				parser.printUsage(System.out);
				return;
			}

			FileWriter fileOut=null;

			if(!sFileOut.isEmpty())
				fileOut=new FileWriter(new File(sFileOut));

			BufferedWriter writer=new BufferedWriter(new OutputStreamWriter(System.out));
			XmlWriter xml=new XmlWriter(writer);
			xml.writeEntity("Classification");
			xml.writeAttribute("Pattern",sPattern);
			xml.writeEntity("Note");
			xml.writeText("Classify input structures based on morphology. It will output the classification results to screen and output the structures mathching to the pattern if applicant");
			xml.endEntity();

			Cluster mol = new Cluster();
			Scanner scanner = new Scanner(new File(sFileIn));

			HashMap<String, Integer> map=new HashMap<String, Integer>();
			map.put("Total",0);
			int i = 0;
			while (mol.Read(scanner, MTools.getExtension(sFileIn))) {
				//System.out.print(MTools.getExtension(sFileOut));

				String morph=mol.getMorphology(true);
				if(morph.contentEquals(sPattern)){					
					if(!sFileOut.isEmpty())
						mol.Write(fileOut,MTools.getExtension(sFileOut));
				}

				if(map.containsKey(morph)){
					int  newVal= ((Integer)map.get(morph)).intValue();
					map.put(morph, map.get(morph) + 1);
				}else map.put(morph, 1);

				map.put("Total",  map.get("Total") + 1);
				
				xml.writeEntity("Class");
				xml.writeAttribute("id", Integer.toString(i));
				xml.writeAttribute("Morphology", morph);
				xml.writeAttribute("Nulity",Integer.toString(mol.getNumberOfSmallestRingByCauchyFormula(false)));
				xml.writeAttribute("CoordNum",mol.getCoordNum());
				xml.endEntity();
				i++;
				//break;
			}
			xml.writeEntity("Statistics");

			for (String s : map.keySet()) {
				xml.writeAttribute(s,map.get(s).toString());
            }
			xml.endEntity();

			xml.endEntity();
			writer.close();
		} catch (IOException ex) {
			Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		}
		
	}
	

}
