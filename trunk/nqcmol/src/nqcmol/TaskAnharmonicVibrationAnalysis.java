/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.potential.GaussianInterfacePotential;
import nqcmol.potential.Potential;
import nqcmol.tools.MTools;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math.optimization.OptimizationException;
import org.apache.commons.math.optimization.fitting.CurveFitter;
import org.apache.commons.math.optimization.fitting.ParametricRealFunction;
import org.apache.commons.math.optimization.fitting.PolynomialFitter;
import org.apache.commons.math.optimization.general.LevenbergMarquardtOptimizer;
import org.kohsuke.args4j.Option;


/**
 *
 * @author nqc
 */
public class TaskAnharmonicVibrationAnalysis extends Task {

	@Option(name="-a",usage="parameter for Gaussian. It should be in this format \"method@charge@multiplicity\". For example: \"b3lyp\6-31+G*@0@1\" (note the quote)",metaVar="STRING")
    String sFileParam="hf/31G@1@1";

	private String basisSet="";
	private String chargeAndMultiplicity="0 1";
	private String sFileCheckPoint="checkpoint.chk";

	Potential pot;
	GaussianInterfacePotential gau=null;

	@Override
	public String getName(){
		return "AnharmonicVibrationAnalysis";
	}

	@Override
	protected void Initialize() {
		try {
			super.Initialize();
			pot = MolExtra.SetupPotential("g03", sFileParam, "Hartree", "");
			assert (pot.getEquation().contains("Gaussian"));
			gau = (GaussianInterfacePotential) pot;
			xmllog.writeNormalText(gau.Info(1));
			basisSet = gau.getBasisSet();
			chargeAndMultiplicity = gau.getChargeAndMultiplicity();
			String filename = (new File(sFileIn)).getName();
			String ext = (filename.lastIndexOf(".") == -1) ? "" : filename.substring(filename.lastIndexOf(".") + 1, filename.length());
			sFileCheckPoint = filename.replace("." + ext, "");
		} catch (IOException ex) {
			Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
		}
	}


	@Override
	protected void Process() {
		try {
//			xmllog.writeEntity("Note");
//			xmllog.writeText(" Time is measured in seconds");
//			xmllog.endEntity();
			if (cluster.Read(fileIn, "g03")) {
				//Reading normal modes
				cluster.setTag((new File(sFileIn)).getName());
				xmllog.writeEntity("Cluster").writeAttribute("Tag", cluster.getTag());
				xmllog.writeAttribute("nAtom", Integer.toString(cluster.getNAtoms()));
				xmllog.writeAttribute("Energy", Double.toString(cluster.getEnergy()));
				if(cluster.getNormalModeVectors()!=null){
					
					if(iMode==-1){
						for(int i=0;i<cluster.getFreqs().length;i++){
							CalculateAnharmonicityFreq(i);
						}
					}else{
						int iFrom=iMode;
						int iTo= Math.min(iMode + nModes,cluster.getFreqs().length);
						for(int i=iFrom;i<iTo;i++){
							CalculateAnharmonicityFreq(i);
						}
					}
				}

				xmllog.endEntity().flush();
			}
			fileIn.close();
		} catch (IOException ex) {
			Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
		}
	}


	///===================
	@Option(name = "-np", usage = "number of processors using in Gaussian. Default is 1.", metaVar = "INTEGER")
	int numOfProcessors=1;

	@Option(name = "-num", usage = "number of points using in approximation. Default is 5.", metaVar = "INTEGER")
	int nPoints=5;
	double[] delta=null;
	double[] energies=null;

	@Option(name = "-m", usage = "starting mode to be calculated (count from 0), number of calculated modes is determined from -nmodes option. If -1 (default), all normal modes will be calculated.", metaVar = "INTEGER")
	int iMode=-1;

	@Option(name = "-nmodes", usage = "number of modes will be calculated (default is 1). The starting mode is determined from -m option", metaVar = "INTEGER")
	int nModes=1;

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
				long duration =  System.currentTimeMillis();
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
				//System.err.println(sGaussianInput);
				gau.setCluster(cluster);
				energies = gau.getEnergies(sGaussianInput);
				for (int j = 0; j < energies.length; j++) {
					xmllog.writeEntity("Step").writeAttribute("id", Integer.toString(j+1));
					xmllog.writeAttribute("DeltaX", Double.toString(delta[j]));
					xmllog.writeAttribute("Energy", Double.toString(energies[j]));
					xmllog.endEntity().flush();
				}
				xmllog.endEntity().flush();

				double[] params=null,calcE=null;
				double anFreq=0;

				//Polynomial approximation
				xmllog.writeEntity("PolynomialFitting");
					PolynomialFit(params,calcE,degree);
						double reducedMass=cluster.getReducedMass(m);//in amu
						double forceConst=params[2]*2;//in Hartree.Angstrom^-2
						double kConvert=Math.sqrt(5.4857990943e-4)*1e8*7.297352569e-3;
						anFreq=kConvert*Math.sqrt(forceConst/reducedMass)/(2.0*Math.PI);//in cm-1
					xmllog.writeAttribute("AnFreq",Double.toString(anFreq));

					//details of parameters
					for (int j = 0; j < params.length; j++) {
						xmllog.writeEntity("Term").writeAttribute("Degree", Integer.toString(j));
						xmllog.writeAttribute("Coefficient", Double.toString(params[j]));
						xmllog.endEntity().flush();
					}
					//writing validity of fitting: rmsE etc
					WriteValidity(calcE);
				xmllog.endEntity().flush();


				//Morse fitting
//				xmllog.writeEntity("MorsePotentialFitting");
//					MorseFit(params,calcE);
//					double De=params[0];
//					double alpha=params[1];
//					//results
//					xmllog.writeAttribute("AnFreq",Double.toString(anFreq));
//
//					//details of parameters
//					xmllog.writeAttribute("De", Double.toString(params[0]));
//					xmllog.writeAttribute("Alpha", Double.toString(params[1]));
//
//					//writing validity of fitting: rmsE etc
//					WriteValidity(calcE);
//
//				xmllog.endEntity().flush();

				xmllog.writeEntity("Duration");
					xmllog.writeNormalText(Double.toString((duration-System.currentTimeMillis())/1000.0));
				xmllog.endEntity().flush();
			}


		xmllog.endEntity().flush();
	}

	double WriteValidity(double[] calcE){
		double rmsE = 0;
		try {
			xmllog.writeEntity("Validity");
			for (int j = 0; j < energies.length; j++) {
				double deltaE = Math.abs(calcE[j] - energies[j]);
				rmsE += deltaE * deltaE;
			}
			rmsE = Math.sqrt(rmsE / energies.length);
			xmllog.writeAttribute("RMSE", Double.toString(rmsE));
			for (int j = 0; j < energies.length; j++) {
				double deltaE = Math.abs(calcE[j] - energies[j]);
				xmllog.writeEntity("Step").writeAttribute("id", Integer.toString(j));
				xmllog.writeAttribute("DeltaX", Double.toString(delta[j]));
				xmllog.writeAttribute("ObservedE", Double.toString(energies[j]));
				xmllog.writeAttribute("CalcE", Double.toString(calcE[j]));
				xmllog.writeAttribute("DeltaE", Double.toString(deltaE));
				xmllog.endEntity().flush();
			}
			xmllog.endEntity().flush();
		} catch (IOException ex) {
			Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
		}
		return rmsE;
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

			String g03header="%chk="+sFileCheckPoint+String.format("-m%d.chk", m);
			g03header+=" \n%mem=1GB \n%nproc="+String.format("%d", numOfProcessors)+" \n#p "+basisSet;
			//g03header+=" \n%nproc="+String.format("%d", numOfProcessors)+" \n#p "+basisSet;


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

	/**
	 * Fitting polynomial potential with observation (x,y)=(delta,energies)
	 * @param params coefficients of Morse potential, will be initialized and updated after fitting
	 * @param calcE values of energies calculated by new polynomial,will be initialized and updated after fitting
	 * @param degree maximum degree of polynomial
	 */
	void PolynomialFit(double[] params,double[] calcE,int degree) {
			PolynomialFitter fitter = new PolynomialFitter(degree, new LevenbergMarquardtOptimizer());
			fitter.addObservedPoint(1, 0, cluster.getEnergy());
			for (int i = 0; i < delta.length; i++) {
				fitter.addObservedPoint(1, delta[i], energies[i]);
			}
			PolynomialFunction fittedFunc=null;
			try {
				fittedFunc = fitter.fit();
				params=fittedFunc.getCoefficients();
				calcE=new double[energies.length];
				for (int j = 0; j < energies.length; j++) {
					calcE[j]=fittedFunc.value(delta[j]);
				}
			} catch (OptimizationException ex) {
				Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
			}
			
	}

	/**
	 * Fitting Morse potential with observation (x,y)=(delta,energies)
	 * @param params params of Morse potential, will be initialized and updated after fitting
	 * @param calcE values of energies calculated by Morse potential,will be initialized and updated after fitting
	 */
	void MorseFit(double[] params,double[] calcE) {
			CurveFitter fitter = new CurveFitter(new LevenbergMarquardtOptimizer());
			fitter.addObservedPoint(1, 0, cluster.getEnergy());
			for (int i = 0; i < delta.length; i++) {
				fitter.addObservedPoint(1, delta[i], energies[i]);
			}
			MorsePotential fittedFunc=new MorsePotential();
			double[] inputParam={1,1};

			try {
				params = fitter.fit(fittedFunc, inputParam);
				calcE=new double[energies.length];
				for (int j = 0; j < energies.length; j++) {
					calcE[j]=fittedFunc.value(delta[j],params);
				}
			} catch (OptimizationException ex) {
				Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
			} catch (FunctionEvaluationException ex) {
				Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
			} catch (IllegalArgumentException ex) {
				Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
			}
	}

	class MorsePotential implements ParametricRealFunction{

		@Override
		public double value(double x, double[] param) throws FunctionEvaluationException {
			double De=param[0];
			double a=param[1];
			double re=0;
			return De*Math.pow(1-Math.exp(a*(re-x)),2);
		}

		@Override
		public double[] gradient(double x, double[] param) throws FunctionEvaluationException {
			double De=param[0];
			double a=param[1];
			double re=0;
			double[] gradient=new double[2];

			gradient[0]=Math.pow(1-Math.exp(a*(re-x)),2);
			gradient[1]=2.0*De*(1-Math.exp(a*(re-x)))*(-(re-x)*Math.exp(a*(re-x)));
			return gradient;
		}
		
	}

	
}
