/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.potential.GaussianInterfacePotential;
import nqcmol.tools.MTools;
import org.apache.commons.math.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math.optimization.OptimizationException;
import org.apache.commons.math.optimization.fitting.PolynomialFitter;
import org.apache.commons.math.optimization.general.LevenbergMarquardtOptimizer;
import org.kohsuke.args4j.Option;


/**
 *
 * @author nqc
 */
public class TaskAnharmonicVibrationAnalysis extends TaskCalculate {
	private String basisSet="";
	private String chargeAndMultiplicity="0 1";
	private String sFileCheckPoint="checkpoint.chk";

	GaussianInterfacePotential gau=null;



	@Override
	public String getName(){
		return "AnharmonicVibrationAnalysis";
	}

	@Override
	protected void Initialize() {
		super.Initialize();
		
		pot = MolExtra.SetupPotential(sPotential, sFileParam, sUnit,sMethod);
		assert(pot.getEquation().contains("Gaussian"));
		gau=(GaussianInterfacePotential)pot;

		basisSet=gau.getBasisSet();
		chargeAndMultiplicity=gau.getChargeAndMultiplicity();

		String filename=(new File(sFileIn)).getName();
		String ext=(filename.lastIndexOf(".")==-1)?"":filename.substring(filename.lastIndexOf(".")+1,filename.length());
		sFileCheckPoint=filename.replace("."+ext,"");
	}


	@Override
	protected void Process() {
		try {
//			xmllog.writeEntity("Note");
//			xmllog.writeText(" Time is measured in seconds");
//			xmllog.endEntity();
			if (cluster.Read(fileIn, "g03")) {
				//Reading normal modes
				xmllog.writeEntity("Cluster").writeAttribute("Tag", cluster.getTag());
				xmllog.writeAttribute("nAtom", Integer.toString(cluster.getNAtoms()));
				xmllog.writeAttribute("Energy", Double.toString(cluster.getEnergy()));
				if(cluster.getNormalModeVectors()!=null){
					
					if(iMode==-1)
						for(int i=0;i<cluster.getFreqs().length;i++){
							CalculateAnharmonicityFreq(i);
						}
					else
						CalculateAnharmonicityFreq(iMode);
				}

				xmllog.endEntity().flush();
			}
			fileIn.close();
		} catch (IOException ex) {
			Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
		}
	}


	///===================


	@Option(name = "-np", usage = "number of points using in approximation", metaVar = "INTEGER")
	int nPoints=5;
	double[] delta=null;
	double[] energies=null;

	@Option(name = "-mode", usage = "specify mode to be analyzed (count from 0). if -1, all normal modes will be approximated.", metaVar = "INTEGER")
	int iMode=0;

	@Option(name = "-degree", usage = "Maximum degree in approximation. Default is 2.", metaVar = "INTEGER")
	int degree=2;


	@Option(name = "-maxX", usage = "maximum displacement in approximation", metaVar = "INTEGER")
	double deltaX=0.2;

	Cluster cluster=new Cluster();

	void CalculateAnharmonicityFreq(int m) throws IOException{
		xmllog.writeEntity("NormalMode").writeAttribute("id", Integer.toString(m));
			if((m>=cluster.getFreqs().length)||(m<0))
				xmllog.writeText(String.format("Error. Mode %d is not valid because it is not in between [0,%d)",m,cluster.getFreqs().length));
			else{
				xmllog.writeAttribute("Freq", Double.toString(cluster.getFreqs(m)));
				xmllog.writeAttribute("ReducedMass", Double.toString(cluster.getReducedMass(m)));
				xmllog.writeEntity("Vector");
				for (int k = 0; k < cluster.getNAtoms(); k++) {
					xmllog.writeEntity("Atom").writeAttribute("id", Integer.toString(k));
					xmllog.writeAttribute("x", Double.toString(cluster.getNormalModeVectors(m)[k * 3 + 0]));
					xmllog.writeAttribute("y", Double.toString(cluster.getNormalModeVectors(m)[k * 3 + 1]));
					xmllog.writeAttribute("z", Double.toString(cluster.getNormalModeVectors(m)[k * 3 + 2]));
					xmllog.endEntity().flush();
				}
				xmllog.endEntity().flush();
				//generate gaussian input
				xmllog.writeEntity("GaussianCalling");
				String sGaussianInput = GenGaussianInputFromNormalModes(m);
				gau.setCluster(cluster);
				energies = gau.getEnergies(sGaussianInput);
				for (int j = 0; j < energies.length; j++) {
					xmllog.writeEntity("Step").writeAttribute("id", Integer.toString(j));
					xmllog.writeAttribute("DeltaX", Double.toString(delta[j]));
					xmllog.writeAttribute("Energy", Double.toString(energies[j]));
					xmllog.endEntity().flush();
				}
				xmllog.endEntity().flush();

				//Morse potential approximation
				xmllog.writeEntity("PolynomialFitting");
				PolynomialFunction fitted=PolynominalFit(degree);				
				if(fitted!=null){
					double[] coefs = fitted.getCoefficients();
					for (int j = 0; j < coefs.length; j++) {
						xmllog.writeEntity("Term").writeAttribute("Degree", Integer.toString(j));
						xmllog.writeAttribute("Coefficient", Double.toString(coefs[j]));
						xmllog.endEntity().flush();
					}
					
					double rmsE=0;
					xmllog.writeEntity("Validity");
					for (int j = 0; j < energies.length; j++) {
						double calcE=fitted.value(delta[j]);
						double deltaE=Math.abs(calcE-energies[j]);
						rmsE+=deltaE*deltaE;
						xmllog.writeEntity("Step").writeAttribute("id", Integer.toString(j));
						xmllog.writeAttribute("DeltaX", Double.toString(delta[j]));
						xmllog.writeAttribute("ObservedE", Double.toString(energies[j]));
						xmllog.writeAttribute("CalcE", Double.toString(calcE));
						xmllog.writeAttribute("DeltaE", Double.toString(deltaE));
						xmllog.endEntity().flush();
					}
					xmllog.writeEntity("RMS");
						xmllog.writeNormalText(Double.toString(rmsE));
					xmllog.endEntity();
					xmllog.endEntity().flush();

					///

					xmllog.writeEntity("Results");
						double reducedMass=cluster.getReducedMass(m);//in amu
						double forceConst=coefs[2]*2;//in Hartree.Angstrom^-2

						double kConvert=Math.sqrt(5.4857990943e-4)*1e8*7.297352569e-3;
						double anFreq=kConvert*Math.sqrt(forceConst/reducedMass)/(2.0*Math.PI);//in cm-1
						xmllog.writeAttribute("AnFreq",Double.toString(anFreq));
					xmllog.endEntity();
				}				
				xmllog.endEntity().flush();
			}

		xmllog.endEntity().flush();
	}

	String GenGaussianInputFromNormalModes(int m){
		double[] coord=cluster.getCoords();
		double[] direction=cluster.getNormalModeVectors()[m];


		String input="";
		delta=new double[nPoints];
		for(int i=0;i<nPoints;i++){
			//temporary scaling factor
			double beta=(2.0*i/(nPoints-1)-1)*deltaX;
			delta[i]=beta;

			String g03header="%chk="+sFileCheckPoint+String.format("-m%d.chk", m)+" \n%mem=1GB \n%nproc=2 \n#p "+basisSet;


			if(i!=0)	g03header="--Link1--\n"+	g03header + " Guess=Read ";
			input+=g03header+"\n";
			input+=String.format("\nM= %d i= %d beta= %1.5f\n\n",m,i,beta);

			input+=chargeAndMultiplicity+"\n";//Charge and multiplicity

			double[] newCoords=new double[coord.length];

			MTools.VEC_PLUS_VEC(newCoords, coord, direction, 1.0, beta);

			for(int j=0;j<cluster.getNAtoms();j++){
				input+=cluster.getAtomicSymbol(j)+" ";
				input+=String.format("%12.8f %12.8f %12.8f\n",newCoords[j*3+0],newCoords[j*3+1],newCoords[j*3+2]);
			}

			input+="\n";
		}
		
		return input;
	}

	PolynomialFunction PolynominalFit(int degree) {
			PolynomialFitter fitter = new PolynomialFitter(degree, new LevenbergMarquardtOptimizer());
			fitter.addObservedPoint(1, 0, cluster.getEnergy());
			for (int i = 0; i < delta.length; i++) {
				fitter.addObservedPoint(1, delta[i], energies[i]);
			}
			PolynomialFunction fitted=null;
			try {
				fitted = fitter.fit();
			} catch (OptimizationException ex) {
				Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
			}
			return fitted;

	}
	
}
