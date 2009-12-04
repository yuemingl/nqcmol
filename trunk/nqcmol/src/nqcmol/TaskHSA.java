/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.HashMap;
import java.util.List;
import java.util.Scanner;
import java.util.Vector;
import java.util.logging.*;
import nqcmol.tools.MTools;
import nqcmol.tools.XmlWriter;
import org.kohsuke.args4j.*;



/*
class PartFunc  extends Cluster {
	//enum REGIME_FUNC{M_CLASSICAL=0, M_QUANTUM=1, M_CLASS_NOFREQ=2};
	//enum UNIT{U_HARTREE=0,U_KCAL=1,U_KJ=2};



	@Override
	public void Clear(){ super.Clear(); group=0; mi=1;lnz=0;};

	double calclnZ(double T){
		int i,j;

		lnz=-1;

		double Na=Math.log(2.0)-Math.log(mi);
		for(i=0;i<cTypeMax;i++)	for(j=1;j<nType[i];j++)	Na+=Math.log(j);

		double E_part=0,Freq_part=0;
		calcZeroE();

		switch(regime){
			case M_CLASSICAL:
				Freq_part+=-Math.log(cBetaPlanck/T)*nFreq;
				for(i=0;i<nFreq;i++){
					if(freq[i]>0){
						Freq_part+=-Math.log(freq[i]); //cout<<inv->freq[i]/T<<endl;
					}else{
						cerr<<" PartFuncPop "<<i<<"th has a negative frequency ="<<freq[i]<<". I will stop here."<<endl;
						return -1;
					}
				}
			break;

			case M_QUANTUM:
				for(i=0;i<nFreq;i++)
					if(freq[i]>0){
						Freq_part+= -Math.log(1.0-exp(-cBetaPlanck*freq[i]/T));
						//cout<<inv->freq[i]/T<<endl;
					}else{
						cerr<<" PartFuncPop "<<i<<"th has a negative frequency ="<<freq[i]<<". I will stop here."<<endl;
						return -1;
					}

			break;

			case M_CLASS_NOFREQ:
			break;
		}
		Ea=energy -ref + ZeroE;

		E_part=-beta*Ea/T;

		lnz=Na + E_part + Freq_part;
	#if DEBUG
			cout<<"Energy("<<Eunit<<") = "<<energy<<" ref = "<<ref<<" at "<<T<<"K ";
			cout<<"Ea= "<<Ea<<" Na= "<<Na<<" ZP= "<<ZeroE<<" Epart= "<<E_part<<" Freq_part= "<<Freq_part<<" lnZ(T)= "<<lnz<<endl;
	#endif

		return lnz;

	}

	double calcZeroE();

	CString calcPointGroup(){
		Genes_str p;
		string SymCode,GroupName;
		Migrate(p);
		SymCalc( p,mi, GroupName,SymCode);
		return GroupName;
	}
	
	//static double InitConst(double ref_=0,int Eunit_=U_HARTREE,int regime_=M_CLASSICAL);

	int group;
	int mi;
	double lnz;

	double ZeroE;
	double Ea; //energy with zero point energy correction

	static double ref=0;
	static int Eunit;
	static int regime;
	static double beta=315774.659661423;
	static double planck=4.55633500389233e-06;

	static final double cBeta[]={315774.659661423,503.219018217726,120.272231887602};//for hartree,kcal/mol, kJ/mol
	static final double cPlanck[]={4.55633500389233e-06,0.00285914300348446,0.0119626543265790};//for hartree,kcal/mol, kJ/mol
	static final double cBetaPlanck=1.43877513515753; //cm*K
};


/**
 *
 * @author nqc
 */
public class TaskHSA extends Task{
	/*
	@Override
	public String getName(){
		return "ClassifyCluster";
	}

	@Override
	protected void Initialize() {
		super.Initialize();
		try {

			xmllog.writeAttribute("Pattern", sPattern);
			xmllog.writeEntity("Note");
			xmllog.writeText("Classify input clusters based on morphology. It will output the classification results to screen and output the structures mathching to the pattern if applicant.\n");
			xmllog.endEntity();
			if (bWaterCoordNum) {
				xmllog.writeText("Classify input water clusters based on coordination number of protonated or deprotonated oxygen atom.");
			}
		} catch (IOException ex) {
			Logger.getLogger(TaskRemoveDuplicateCluster.class.getName()).Math.log(Level.SEVERE, null, ex);
		}
	}

	@Override
	protected void Process() {
		try {
			FileWriter fileOut = null;
			if (!sFileOut.isEmpty()) {
				fileOut = new FileWriter(new File(sFileOut));
			}
						
			WaterCluster mol = new WaterCluster();
			Scanner scanner = new Scanner(new File(sFileIn));
			HashMap<String, Integer> statCount = new HashMap<String, Integer>();
			statCount.put("Total", 0);
			double totalCoordNum = 0;
			double[] avgCoordNum = new double[Cluster.cNumCountMax];
			MTools.VEC_EQU_NUM(avgCoordNum, 0);
			double avgCompactness = 0;
			int i = 0;
			while (mol.Read(scanner, sFormatIn)) {
				//System.out.print(MTools.getExtension(sFileOut));
				mol.getPairwiseBond();
				String morph = mol.getMorphology(false);
				if (morph.contentEquals(sPattern)) {
					if (!sFileOut.isEmpty()) {
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
				//average compactness
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
		} catch (IOException ex) {
			Logger.getLogger(TaskHSA.class.getName()).Math.log(Level.SEVERE, null, ex);
		}		
		
		
	}

*/
}
