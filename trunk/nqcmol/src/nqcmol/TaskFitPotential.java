/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
import nqcmol.cluster.MolExtra;
import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.DifferentiableMultivariateVectorialFunction;
import org.apache.commons.math.analysis.MultivariateMatrixFunction;
import org.apache.commons.math.optimization.OptimizationException;
import org.apache.commons.math.optimization.VectorialPointValuePair;
import org.apache.commons.math.optimization.general.LevenbergMarquardtOptimizer;
import org.kohsuke.args4j.*;


import org.jgap.*;
import org.jgap.impl.*;

/**
 *
 * @author nqc
 */
public class TaskFitPotential extends TaskCalculate {		
	

	@Option(name="-ref",usage="To calculate binding energies. It MUST be in F_XYZ format and in the following order: neutral, protonated, deprotonated. ",metaVar="FILE")
	String sFileRef="";

    @Option(name="-random",usage="Randomize the input parameters withing the range")
	boolean isRandom=false;

    @Option(name="-nOpts",usage="Maximum number of iterations in local optimization. [0]",metaVar="INTEGER")
	int nOpts=0;
		
	final int cParamMax=100;

    @Option(name="-b",usage="Bound file. Each line should consist of at least tree columns: upper lower fixed.[]",metaVar="STRING")
	String sFileBound="";

	int nparam; //!< number of non-fixed parameters
	double[] param; //!< non-fixed parameters
	double[] ub;//!< non-fixed upper bounds, all fixed parameters are removed
	double[] lb;//!< non-fixed lower bounds, all fixed parameters are removed
	String[] label; //!< label of non-fixed parameters

	//int param_inp.length; //!< number of parameters (in total)
	double[] param_inp;//=new double[cParamMax]; //!< starting value of parameters
	double[] ub_inp;//=new double[cParamMax];//!< input upper and lower bounds
	double[] lb_inp;//=new double[cParamMax];
	int[] fixed;//=new int[cParamMax]; //!< determine wheter parameter are fixed
	String[] label_inp=new String[cParamMax]; //!< label of parameter

	Vector<Cluster> molRef= new Vector<Cluster>();//!< reference clusters
	Vector<Cluster> data=new Vector<Cluster>();//!< data/fitting clusters
	Vector<Double> weight=new Vector<Double>();//!< weighted numbers of each clusters

	double[] target; //!< contain target data, in that case is binding energy +rms

	Vector<String> datafile=new Vector<String>();
	Vector<Double> dataWeight=new Vector<Double>();


    int nIterations=0;
    int nEvaluations=0;
    double finalRMS=0;


    @Override
	public String getName() {
		return "FitPotential";
	}

    static final public String Option="fit";

    static final public String Descriptions="\t "+Option+" \t - "+ "Fit potential by using GA\n";
	
	void ProcessData(){
        double[] ener_ref={molRef.get(0).getEnergy(),molRef.get(1).getEnergy(),molRef.get(2).getEnergy()};
		target=new double[data.size()];
		for (int i=0; i< data.size(); i++){
			Cluster m=data.get(i);

			target[i]= CalcBindingEnergy(m,m.getEnergy(),ener_ref);
			//target[i]*=1000;
		}
	}

    private double CalcBindingEnergy(Cluster m,double ener,double[] ener_shift){
        double be=0;
        int nIons=m.getNonHydrogenNum();
        switch(m.getTotalCharge()){
            case  0:    be = ener - ener_shift[0]*(nIons); break;
            case  1:    be = ener - ener_shift[0]*(nIons-1) -     ener_shift[1]; break;
            case -1:    be = ener - ener_shift[0]*(nIons-1) -     ener_shift[2]; break;
            case  2:    be = ener - ener_shift[0]*(nIons-2) - 2.0*ener_shift[1]; break;
        }
        return be;
    }


	@Override
	protected void Initialize(){
	try {
		xmllog.writeAttribute("InputFile", sFileIn).writeAttribute("FormatInput", sFormatIn);
        xmllog.writeAttribute("OutputFile", sFileOut);
        xmllog.writeAttribute("FileBound", sFileBound);
        xmllog.writeAttribute("RandomParam", Boolean.toString(isRandom));
		xmllog.writeAttribute("Verbose", Integer.toString(verbose));

		xmllog.writeEntity("Initialize");
		pot = MolExtra.SetupPotential(sPotential, sFileParam, sUnit,sMethod);
		xmllog.writeNormalText(pot.Info(1));

		xmllog.writeEntity("ReadBound");
        ReadBound();
        xmllog.endEntity().flush();

		xmllog.writeEntity("ReadReference");
			Scanner scanRef = new Scanner(new File(sFileRef));
			molRef.clear();
			while (scanRef.hasNext()) {
				Cluster tmpMol = new Cluster();
				tmpMol.Read(scanRef, "xyz");

				tmpMol.CorrectOrder();
				if(tmpMol.getNAtoms()==0) break;
				molRef.add(tmpMol);

				xmllog.writeEntity("Cluster");
				xmllog.writeAttribute("id", Integer.toString(molRef.size()-1));
				xmllog.writeAttribute("nAtoms", Integer.toString(tmpMol.getNAtoms()));
				xmllog.writeAttribute("Energy", Double.toString(tmpMol.getEnergy()));
				xmllog.endEntity().flush();
			}
			scanRef.close();
		xmllog.endEntity().flush();

		xmllog.writeEntity("ReadData");
		for (int i = 0; i < datafile.size();i++) {
			Scanner scanData = new Scanner(new File(datafile.get(i)));
			int count=0;

			while (scanData.hasNext()) {
				Cluster tmpMol = new Cluster();
				tmpMol.Read(scanData, "xyz");
				if(tmpMol.getNAtoms()==0) break;
				tmpMol.CorrectOrder();

				data.add(tmpMol);
				weight.add(dataWeight.get(i));
				count++;
			}
			xmllog.writeEntity("Data");
			xmllog.writeAttribute("File", datafile.get(i));
			xmllog.writeAttribute("Size", Integer.toString(count));
			xmllog.writeAttribute("Weight", Double.toString(dataWeight.get(i)));
			xmllog.endEntity().flush();
			scanData.close();
		}
		xmllog.endEntity().flush();

		xmllog.writeEntity("TotalSize");
		xmllog.writeNormalText(Integer.toString(data.size()));
		xmllog.endEntity().flush();


		xmllog.writeEntity("ProcessData");
		ProcessData();
		xmllog.endEntity().flush();

	} catch (IOException ex) {
		logger.severe("Cannot write to xmllog");
	}

}


	@Override
	protected void Process(){
		try {
			xmllog.writeEntity("InitialEvaluation");
				Evaluate(param,4);
			xmllog.endEntity().flush();

			xmllog.writeEntity("Optimization");
				long duration=System.currentTimeMillis();

                FitUsingApache();
                //FitUsingGA();
				duration=(System.currentTimeMillis()-duration);
				
			xmllog.endEntity().flush();

            xmllog.writeEntity("ConvergenceInfo");
                    xmllog.writeAttribute("TotalRMS", Double.toString(finalRMS));
                    //xmllog.writeAttribute("MaxInterations", Integer.toString(fit.getMaxIterations()));
                    //xmllog.writeAttribute("MaxEvaluations", Integer.toString(fit.getMaxEvaluations()));
                    xmllog.writeAttribute("Iterations", Integer.toString(nIterations));
                    xmllog.writeAttribute("Evaluations", Integer.toString(nEvaluations));
                    xmllog.writeAttribute("Duration", Double.toString(duration/1000.0));
            xmllog.endEntity().flush();

			xmllog.writeEntity("FinalEvaluation");
				Evaluate(param, 4);
			xmllog.endEntity().flush();
			
			xmllog.writeEntity("FinalParameters");	
				WriteParam(param);

				double[] result=new double[param_inp.length];
				ConvertParam(param,result);

				
			xmllog.endEntity().flush();
		} catch (IllegalArgumentException ex) {
			Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
        }
	}


    protected void ReadBound() throws FileNotFoundException, IOException{
		Scanner fin = new Scanner(new File(sFileBound));
        param_inp=pot.getParam();
        ub_inp=new double[cParamMax];//!< input upper and lower bounds
        lb_inp=new double[cParamMax];
        fixed=new int[cParamMax]; //!< determine wheter parameter are fixed

		int nfixed = 0;

        for(int i=0;i<param_inp.length;i++){
			String line = fin.nextLine();
			if (line.isEmpty()) break;

			StringTokenizer tokenizer;
			tokenizer = new StringTokenizer(line, "\t ,;"); //read no of atoms
			String info;
			//info = tokenizer.nextToken();
			//param_inp[i] = Double.parseDouble(info);
			info = tokenizer.nextToken();			lb_inp[i] = Double.parseDouble(info);
			info = tokenizer.nextToken();			ub_inp[i] = Double.parseDouble(info);
			info = tokenizer.nextToken();			fixed[i] = Integer.parseInt(info);
			label_inp[i] = "";
			while (tokenizer.hasMoreTokens()) {
				label_inp[i] = tokenizer.nextToken();
			}
			//sscanf(line.c_str(),"%lf %lf %lf %d %line", &param_inp[param_inp.length], &lb_inp[param_inp.length], &ub_inp[param_inp.length], &fixed[param_inp.length]);

            if((isRandom)&&(fixed[i]==0)){
                param_inp[i]=lb_inp[i]+Math.random()*(ub_inp[i]-lb_inp[i]);
            }

			nfixed += fixed[i];
			xmllog.writeEntity("Param").writeAttribute("id", Integer.toString(i));
			xmllog.writeAttribute("Value", Double.toString(param_inp[i]));
			xmllog.writeAttribute("Lower", Double.toString(lb_inp[i]));
			xmllog.writeAttribute("Upper", Double.toString(ub_inp[i]));
			xmllog.writeAttribute("Fixed", Integer.toString(fixed[i]));
			xmllog.writeAttribute("Label", label_inp[i]);
			xmllog.endEntity().flush();
		}

		nparam = param_inp.length - nfixed;
		fin.close();
		xmllog.endEntity();

		xmllog.writeEntity("ActualParam");
		xmllog.writeAttribute("NParam", Integer.toString(nparam));
		xmllog.writeAttribute("NFixed", Integer.toString(nfixed));
		param = new double[nparam];
		ub = new double[nparam];
		lb = new double[nparam];
		label = new String[nparam];
		int t = 0;
		for (int i = 0; i < param_inp.length; i++) {
			if (fixed[i] == 0) {
				lb[t] = lb_inp[i];
				ub[t] = ub_inp[i];
				param[t] = param_inp[i];
				label[t] = label_inp[i];
				//cout<<i << ":  " << param_inp[i] << " " << lb_inp[i] << " " << ub_inp[i] << " " <<  fixed[i] << " "<<label_inp[i]<<endl;
				assert (lb[t] <= ub[t]);
				t++;
			}
		}
		xmllog.endEntity().flush();
    }

    protected void WriteParam(double[] param){      
            double[] result = new double[param_inp.length];
            ConvertParam(param, result);

            if (!sFileOut.isEmpty()) {
                  try {
                        fileOut = new FileWriter(new File(sFileOut));
                        pot.setParam(result);
                        pot.writeParam(fileOut);
                        fileOut.close();
                   } catch (IOException ex) {
                       logger.severe("Cannot write to "+sFileOut);
                    }
            }

            for(int i=0;i<result.length;i++){
					xmllog.writeEntity("Param").writeAttribute("id", Integer.toString(i));
					xmllog.writeAttribute("Value", Double.toString(result[i]));
					xmllog.writeAttribute("Lower", Double.toString(lb_inp[i]));
					xmllog.writeAttribute("Upper", Double.toString(ub_inp[i]));
					xmllog.writeAttribute("Fixed", Integer.toString(fixed[i]));
					xmllog.writeAttribute("Label", label_inp[i]);
					xmllog.endEntity().flush();
			}       
    }


    double Evaluate(double[] p,int verbose){
		double[] y=new double[data.size()];
		return Evaluate(p,y,verbose);
	}

	double Evaluate(double[] p,double[] y,int verbose){
        double[] alpha=new double[param_inp.length];
        ConvertParam(p,alpha);
        pot.setParam(alpha);

        if(verbose>=1){
                xmllog.writeEntity("Alpha");
                String s = "";
                for (int i = 0; i < alpha.length; i++) {
                    s += String.format("%f ", alpha[i]);
                }
                xmllog.writeNormalText(s);
                xmllog.endEntity().flush();
        }
        //calculate
        double[] ener_shift=new double[molRef.size()];
        double rms=0;

        for(int i=0; i<molRef.size(); i++){
            Cluster m=molRef.get(i);
            pot.setCluster(m);
            ener_shift[i] = pot.getEnergy(m.getCoords());
        }

        if (verbose>=2){
            xmllog.writeEntity("Reference");
            for (int i = 0; i < ener_shift.length; i++) {
                Cluster m = molRef.get(i);
                xmllog.writeEntity("Cluster");
                xmllog.writeAttribute("id", Integer.toString(i));
                xmllog.writeAttribute("nAtoms", Integer.toString(m.getNAtoms()));
                xmllog.writeAttribute("Energy", Double.toString(ener_shift[i]));
                xmllog.endEntity().flush();
            }
            xmllog.endEntity().flush();
        }

        //cout <<"\n\n@$    inputs# nAtoms Weight     Obsr.Data        Calc.Data    deltaE" << endl;
        for (int i=0; i<data.size(); i++){
            Cluster m=data.get(i);
            //calculate binding energy
            pot.setCluster(m);
            double ener = pot.getEnergy(m.getCoords());

            y[i] = CalcBindingEnergy(m,ener,ener_shift);
            //y[i]*=1000;

            double dE = Math.abs(y[i] - target[i])* Math.sqrt(weight.get(i));
            rms+=dE*dE;

            if(verbose>=3){
                xmllog.writeEntity("Eval");
                xmllog.writeAttribute("id", Integer.toString(i));
                xmllog.writeAttribute("nAtoms", Integer.toString(m.getNAtoms()));
                xmllog.writeAttribute("Weight", Double.toString(weight.get(i)));
                xmllog.writeAttribute("ObservedBE", Double.toString(target[i]));
                xmllog.writeAttribute("CalcE", Double.toString(ener));
                xmllog.writeAttribute("CalcBE", Double.toString(y[i]));
                xmllog.writeAttribute("DeltaE", Double.toString(dE));
                xmllog.endEntity().flush();
            }
        }

        rms=Math.sqrt(rms/data.size());

        if (verbose>=1){
                xmllog.writeEntity("TotalRMS");
                xmllog.writeNormalText(Double.toString(rms));
                xmllog.endEntity().flush();

        }

        return rms;
    }

	private void ConvertParam(double[] actual,double[] p){
		int c=0;
		for(int i=0; i<param_inp.length; i++)
			if (fixed[i]==0){
				p[i] = actual[c];
				c++;
			}else p[i] = param_inp[i];
	}

	@Override
	public void ParseArguments(String[] args) {
		parser = new CmdLineParser(this);
		try {
			parser.parseArgument(args);
			if(isHelp){
				parser.printUsage(System.out);
				return ;
			}

			for(int i=0;i<50;i++){
				boolean bSet=false;
				String opt=String.format("-i%d",i+1);
					for(int j=0;j<args.length-1;j++)
						if(args[j].contentEquals(opt)){
							bSet=true;
							datafile.add(args[j+1]);
						}

				if(bSet){
					double t=1;
					opt=String.format("-w%d",i+1);

					for(int j=0;j<args.length-1;j++)
						if(args[j].contentEquals(opt))
							t=Double.parseDouble(args[j+1]);

					dataWeight.add(t);
				}
			}

			if(!sFileOut.isEmpty()) fileOut=new FileWriter(new File(sFileOut));

		} catch (CmdLineException ex) {
			logger.severe(ex.getMessage());
		} catch (IOException ex) {
			logger.severe(ex.getMessage());
		}

	}


     //=================== for APACHE fitting
    public class APACHEObjectiveFunction implements DifferentiableMultivariateVectorialFunction {//used with apache commom optimization
            public class Jacobian implements MultivariateMatrixFunction {
				@Override
				public double[][] value(double[] a) throws FunctionEvaluationException, IllegalArgumentException {
					double[][] result=new double[data.size()][a.length];
					double ratio=1e-3;

					//one point jacobian
					double[] y0=APACHEObjectiveFunction.this.value(a);
					for(int i=0;i<a.length;i++){
						//double delta=Math.abs(a[i])*ratio;
						double delta=ratio;
						double tmp=a[i];
						a[i]+=delta;
						double[] y2=APACHEObjectiveFunction.this.value(a);
						a[i]=tmp;
						for(int j=0;j<result[i].length;j++)
							result[i][j]=(y2[j]-y0[j])/(delta);
					}

					//two points jacobian
//					for(int i=0;i<a.length;i++){
//						double delta=Math.abs(a[i])*ratio;
//						a[i]-=delta;
//						double[] y1=APACHEObjectiveFunction.this.value(a);
//						a[i]+=2*delta;
//						double[] y2=APACHEObjectiveFunction.this.value(a);
//						a[i]-=delta;
//						for(int j=0;j<result[i].length;j++)
//							result[i][j]=(y2[j]-y1[j])/(2*delta);
//					}


					return result;
				}
			}

		Jacobian jacobian=new Jacobian();

		@Override
		public MultivariateMatrixFunction jacobian() {
			return jacobian;
		}

		@Override
		public double[] value(double[] a) throws FunctionEvaluationException, IllegalArgumentException {
			double[] y=new double[data.size()];
            Evaluate(a, y, 0);
			//System.err.printf(" valuesize=%d\n", y.length);
			return y;
		}

	}

	protected void FitUsingApache(){
        try {
            APACHEObjectiveFunction obj = new APACHEObjectiveFunction();
            LevenbergMarquardtOptimizer fit = new LevenbergMarquardtOptimizer();
            double[] newweight = new double[weight.size()];
            for (int i = 0; i < newweight.length; i++) {
                newweight[i] = weight.get(i);
            }
            //System.err.printf("targetsize=%d \n",target.length);
            fit.setCostRelativeTolerance(1e-14);
            fit.setOrthoTolerance(1e-14);
            fit.setParRelativeTolerance(1e-14);
            fit.setMaxIterations(nOpts);
            fit.setMaxEvaluations(100*nOpts);
            //fitting here
            
            VectorialPointValuePair r = fit.optimize(obj, target, newweight, param);
            
            param = r.getPoint(); //
           
            nIterations = fit.getIterations();
            nEvaluations = fit.getEvaluations();
            finalRMS = fit.getRMS();
            
        } catch (FunctionEvaluationException ex) {            
            logger.severe(ex.getMessage());
        } catch (OptimizationException ex) {
            logger.warning(ex.getMessage());
        } catch (IllegalArgumentException ex) {
            logger.severe(ex.getMessage());
        }
    }

    /*/=================== for LMA fitting
    public class LMAObjectiveFunction extends LMAMultiDimFunction {
        @Override
        public double getY(double[] x, double[] a) {
                double[] alpha=new double[param_inp.length];
                ConvertParam(a,alpha);
                pot.setParam(alpha);

                //calculate
                double[] ener_shift=new double[molRef.size()];

                for(int i=0; i<molRef.size(); i++){
                        Cluster m=molRef.get(i);
                        pot.setCluster(m);
                        ener_shift[i] = pot.getEnergy(m.getCoords());
                }

                int i=0;//=(int)x;
                Cluster m=data.get(i);
                //calculate binding energy
                double ener = pot.getEnergy(m.getCoords());

                int index=MolExtra.indexForBE(m, molRef);
                double y = ener - ener_shift[0]*(m.getNonHydrogenNum()-1) - ener_shift[index];
                return y;
        }

        @Override
        public double getPartialDerivate(double[] x, double[] a, int parameterIndex) {
            double[] t=a.clone();
                double delta=1e-5;
                t[parameterIndex]-=delta;
                double y1=getY(x,t);
                t[parameterIndex]+=2*delta;
                double y2=getY(x,t);

                return (y2-y1)/(2*delta);
        }
	}

    protected void FitUsingLMA(){
        //MultivariateRealOptimizer fit = new NelderMead();
    	//RealPointValuePair r = null; //=new RealPointValuePair();
        //double finalRMS=Math.sqrt(fit.getChiSquare());

        LMA lma = new LMA(
			new LMAObjectiveFunction(),
			new double[] {1, 1, 1},
			new double[][] {
				// y x0 x1
				{0, 2, 6},
				{5, 10, 2},
				{7, 20, 4},
				{9, 30, 7},
				{12, 40, 6}
			}
		);
		lma.fit();

        nIterations=lma.iterationCount;
        finalRMS= Math.sqrt(lma.chi2);

		
//			"param0: " + lma.parameters[0] + ",\n" +
//			"param1: " + lma.parameters[1] + ",\n" +
//			"param2: " + lma.parameters[2]
//		
    }
    */

    //=================== for GA fitting
    public class GAFitnessFunction  extends FitnessFunction{
        @Override
        protected double evaluate(IChromosome indv) {
            double rms=0;
            double[] x=new double[indv.size()];
            for (int i = 0; i < indv.size(); i++) {
                    DoubleGene gene = (DoubleGene)indv.getGene(i);
                    x[i]=gene.doubleValue();
            }
            rms=1./Evaluate(x,0);
            return rms;
        }
    }

    protected void FitUsingGA(){      
        try {
            Configuration conf = new DefaultConfiguration();
            FitnessFunction myFunc =   new GAFitnessFunction();
            conf.setFitnessFunction(myFunc);
            
            conf.setPreservFittestIndividual(true);
            conf.setKeepPopulationSizeConstant(false);
            Gene[] sampleGenes = new Gene[param.length];
            for (int i = 0; i < param.length; i++) {
                    DoubleGene gene = new DoubleGene(conf, lb[i], ub[i]);
                    gene.setAllele(param[i]);
                    sampleGenes[i] = gene;
            }

            IChromosome sampleChromosome = new Chromosome(conf, sampleGenes);
            conf.setSampleChromosome(sampleChromosome);
            // Finally, we need to tell the Configuration object how many
            // Chromosomes we want in our population. The more Chromosomes,
            // the larger number of potential solutions (which is good for
            // finding the answer), but the longer it will take to evolve
            // the population (which could be seen as bad).
            // ------------------------------------------------------------
            conf.setPopulationSize(50);

             // Create random initial population of Chromosomes.
            // ------------------------------------------------
            Genotype population = Genotype.randomInitialGenotype(conf);
               //replace the first chromosome to the input
            population.getPopulation().setChromosome(0, sampleChromosome);

                         
            // Evolve the population. Since we don't know what the best answer
            // is going to be, we just evolve the max number of times.
            // ---------------------------------------------------------------
            int percentEvolution = (int)Math.max(1, nOpts / 100);
            for (int i = 0; i < nOpts; i++) {
              population.evolve();
              if (( i % percentEvolution == 0) || i==(nOpts-1)) {
                  IChromosome fittest = population.getFittestChromosome();
                  double fitness = fittest.getFitnessValue();

                  xmllog.writeEntity("Step").writeAttribute("id", Integer.toString(i));
					xmllog.writeAttribute("RMS", Double.toString(1.0/fitness));
				  xmllog.endEntity().flush();

                  //System.out.println(" size ="+ fittest.size());
                  double[] currentParam=new double[fittest.size()];
                  for (int k = 0; k < fittest.size(); k++) {
                    DoubleGene gene = (DoubleGene)fittest.getGene(k);
                    currentParam[k]=gene.doubleValue();
                   // System.out.println(" param ="+ currentParam[k]);
                  }
                  WriteParam(currentParam);
                }
            }
            

           IChromosome bestSolutionSoFar = population.getFittestChromosome();
           finalRMS=1.0/bestSolutionSoFar.getFitnessValue();
                     
           for (int i = 0; i < bestSolutionSoFar.size(); i++) {
                    DoubleGene gene = (DoubleGene)bestSolutionSoFar.getGene(i);
                    param[i]=gene.doubleValue();
            }
        } catch (InvalidConfigurationException ex) {
            Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
        }
              
    }
}

