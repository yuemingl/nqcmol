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

	@Option(name="-bWater",usage="classify according to coordination number (for protonated/deprotoanted water cluster")
    boolean bWaterCoordNum=false;
	
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
			xml.writeEntity("Classify");
			xml.writeAttribute("Pattern",sPattern);
			xml.writeEntity("Note");
			xml.writeText("Classify input clusters based on morphology. It will output the classification results to screen and output the structures mathching to the pattern if applicant.\n");
			if(bWaterCoordNum)
				xml.writeText("Classify input water clusters based on coordination number of protonated or deprotonated oxygen atom.");
			xml.endEntity();

			WaterCluster mol = new WaterCluster();
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
					map.put(morph, map.get(morph) + 1);
				}else map.put(morph, 1);


				int iChargedOxygen=mol.getChargedOxygen(false);
				int nCoord=-1;
				if(bWaterCoordNum){
					if(iChargedOxygen!=-1){						
						nCoord=mol.getCoordNum(iChargedOxygen);
						String key="Coord_"+Integer.toString(nCoord);

						if(map.containsKey(key)){
							map.put(key, map.get(key) + 1);
						}else map.put(key, 1);
					}
				}


				map.put("Total",  map.get("Total") + 1);
				
				xml.writeEntity("Cluster");
				xml.writeAttribute("id", Integer.toString(i));
				xml.writeAttribute("Morphology", morph);
				xml.writeAttribute("Nulity",Integer.toString(mol.getNumberOfSmallestRingByCauchyFormula(false)));
				xml.writeAttribute("CoordNumCount",mol.getCoordNumCount());

				if(bWaterCoordNum){
					xml.writeAttribute("ChargedOxygen", Integer.toString(iChargedOxygen));
					xml.writeAttribute("CoordNum",Integer.toString(nCoord));
				}
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
