/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.WaterCluster;
import nqcmol.cluster.Cluster;
import java.io.*;
import java.util.HashMap;
import java.util.Scanner;
import nqcmol.tools.MTools;
import org.kohsuke.args4j.*;

/**
 *
 * @author nqc
 */
public class TaskClassifyCluster extends Task{
	@Option(name="-pattern",usage="pattern of morphology to filter: Cage, MultiRing, DoubleRing, SingleRing, TreeLike, Liner",metaVar="String")
    String sPattern="";

	@Option(name="-bWater",usage="classify according to coordination number (for protonated/deprotonated water cluster")
    boolean bWaterCoordNum=false;

    @Option(name="-clustering",usage="generate hierachical clustering",metaVar="STRING")
    String sFileXMLClustering="";

	@Override
	public String getName(){
		return "ClassifyCluster";
	}

    static final public String Option="class";

    static final public String Descriptions="\t "+Option+" \t - "+ "Classify clusters according to topologies\n";

    public static void main(String[] args) throws IOException  {
        new TaskClassifyCluster().Execute(args);
	}

	@Override
	protected void Initialize() {
		super.Initialize();
        xmllog.writeAttribute("Pattern", sPattern);
        xmllog.writeEntity("Note");
        xmllog.writeText("Classify input clusters based on morphology. It will output the classification results to screen and output the structures mathching to the pattern if applicant.\n");
        xmllog.endEntity();
        if (bWaterCoordNum) {
            xmllog.writeText("Classify input water clusters based on coordination number of protonated or deprotonated oxygen atom.");
        }
	}

	@Override
	protected void Process() {					
        try {
           

            

            HashMap<String, Integer> statCount = new HashMap<String, Integer>();
            statCount.put("Total", 0);

            double totalCoordNum = 0;
            double[] avgCoordNum = new double[Cluster.cNumCountMax];
            MTools.VEC_EQU_NUM(avgCoordNum, 0);
            double avgCompactness = 0;
            int i = 0;

            FileWriter fileOut=null;
            if (!sFileOut.isEmpty()) {
                fileOut = new FileWriter(new File(sFileOut));
            }

            WaterCluster mol = new WaterCluster();
            
            Scanner fileIn = new Scanner(new File(sFileIn));
            while (mol.Read(fileIn, sFormatIn)) {
                mol.getPairwiseBond();
                String morph = mol.getMorphology(false);
                if (morph.contentEquals(sPattern)) {
                    if (fileOut!=null) {
                        mol.Write(fileOut, sFormatOut);
                    }
                }
                
                if (statCount.containsKey(morph)) {
                    statCount.put(morph, statCount.get(morph) + 1);
                } else {
                    statCount.put(morph, 1);
                }
                int iChargedOxygen = mol.getChargedOxygen(false);
                int nCoord = -1;
                if (bWaterCoordNum) {
                    if (iChargedOxygen != -1) {
                        nCoord = mol.getCoordNum(iChargedOxygen);
                        String key = "Coord_" + Integer.toString(nCoord);
                        if (statCount.containsKey(key)) {
                            statCount.put(key, statCount.get(key) + 1);
                        } else {
                            statCount.put(key, 1);
                        }
                    }
                }
                double compactness = mol.getCompactness();
                avgCompactness += compactness;
                //average coordination number
                int[] molCoordNumCount = mol.getCoordNumCount();
                for (int j = 0; j < avgCoordNum.length; j++) {
                    totalCoordNum += molCoordNumCount[j];
                    avgCoordNum[j] += molCoordNumCount[j];
                }
                statCount.put("Total", statCount.get("Total") + 1);
                xmllog.writeEntity("Cluster");
                xmllog.writeAttribute("id", Integer.toString(i));
                xmllog.writeAttribute("Tag", mol.getTag());
                xmllog.writeAttribute("Compactness", Double.toString(compactness));
                xmllog.writeAttribute("Morphology", morph);
                xmllog.writeAttribute("Nulity", Integer.toString(mol.getNumberOfSmallestRingByCauchyFormula(false)));
                xmllog.writeAttribute("CoordNumCount", MTools.VEC_TO_STRING(mol.getCoordNumCount()));
                if (bWaterCoordNum) {
                    xmllog.writeAttribute("ChargedOxygen", Integer.toString(iChargedOxygen));
                    xmllog.writeAttribute("CoordNum", Integer.toString(nCoord));
                }
                xmllog.endEntity().flush();
                i++;
                //break;
            }

            fileIn.close();
            fileOut.close();

            xmllog.writeEntity("Statistics");
            for (String s : statCount.keySet()) {
                xmllog.writeAttribute(s, statCount.get(s).toString());
            }
            xmllog.writeAttribute("AvgCompactness", Double.toString(avgCompactness / i));
            MTools.VEC_MUL_NUM(avgCoordNum, avgCoordNum, 1 / totalCoordNum);
            for (int j = 0; j < avgCoordNum.length; j++) {
                xmllog.writeAttribute("AvgCoordNum_" + Integer.toString(j), Double.toString(avgCoordNum[j]));
            }
            xmllog.endEntity();

            //hierrachical clustering is performed
            if(!sFileXMLClustering.isEmpty()){

            }

        } catch (FileNotFoundException ex) {
            logger.severe(ex.getMessage());
        }  catch (IOException ex) {
            logger.severe(ex.getMessage());
        }
	}

}
