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
import nqcmol.potential.Potential;
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

    @Option(name="-b",usage="Bound file. Each line should consist of at least three columns including [upper] [lower] [fixed].[]",metaVar="STRING")
	String sFileBound="";

    @Option(name="-weight",usage="Using linear weighted scheme.")
	boolean bWeighted=false;



    //@Option(name="-mpi",usage="Using MPI to run parallel or not.[false]",metaVar="STRING")
	boolean bMPI=false;

	

    @Override
	public String getName() {
		return "FitPotential";
	}

    static final public String Option="fit";

    static final public String Descriptions="\t "+Option+" \t - "+ "Fit potential by using GA\n";

    public static void main(String[] args) throws IOException  {
        new TaskFitPotential().Execute(args);
	}

//    @Override
//    public void Execute(String[] args) {
//         for(int i=0;i<args.length;i++){
//            if(args[i].contentEquals("-mpi")) bMPI=true;
//        }
//
//        if(!bMPI) super.Execute(args);
//        else{
//            String[] newargs=MPI.Init(args);
//            rank = MPI.COMM_WORLD.Rank();
//            ncpu = MPI.COMM_WORLD.Size();
//            ParseArguments(newargs);
//
//            try {
//                String sFileLog="main-t"+rank+".xml";
//                stdwriter=new BufferedWriter(new FileWriter(new File(sFileLog)));
//                xmllog = new XmlWriter(stdwriter);
//                System.out.println("Hi from <>"+rank);
//
//                xmllog.writeEntity(getName());
//                Initialize();
//                if(rank==MASTER) {//if u r MASTER, proceed
//                    Process();
//                }else{
//                    while(true){//if u r slaves, just waiting for the orders
//                        mpi.Status status=MPI.COMM_WORLD.Probe(MPI.ANY_SOURCE,MPI.ANY_TAG);
//                        if(status.tag==-1)                         break;
//                        if(status.tag==1){
//                            //MPI.COMM_WORLD.Recv(&a[0],mess.size, newtype,mess.destfrom,mess.tag);
//                            //return true;
//                        }
//                    }
//                }
//                Finalize();
//                xmllog.endEntity().close();
//                stdwriter.close();
//            } catch (IOException ex) {
//                logger.severe(ex.getMessage());
//            }
//            MPI.Finalize();
//        }
//    }


	
	@Override
	protected void Initialize(){
        try {
            xmllog.writeAttribute("InputFile", sFileIn).writeAttribute("FormatInput", sFormatIn);
            xmllog.writeAttribute("OutputFile", sFileOut).writeAttribute("FileBound", sFileBound);
            xmllog.writeAttribute("RandomParam", Boolean.toString(isRandom));
            xmllog.writeAttribute("Verbose", Integer.toString(verbose));

            xmllog.writeEntity("Initialize");
                pot = MolExtra.SetupPotential(sPotential, sFileParam, sUnit,sMethod);
                potpool=new Potential[ncpu];
                for(int i=0;i<ncpu;i++) potpool[i]=MolExtra.SetupPotential(sPotential, sFileParam, sUnit,sMethod);
                xmllog.writeNormalText(pot.XMLInfo(1));

                xmllog.writeEntity("ReadBound");
                    ReadBound();
                xmllog.endEntity().flush();

                xmllog.writeEntity("ReadReference");
                    ReadReference();
                xmllog.endEntity().flush();

                xmllog.writeEntity("ReadData");
                    ReadData();
                    xmllog.writeEntity("TotalSize");
                    xmllog.writeNormalText(Integer.toString(data.size()));
                    xmllog.endEntity().flush();
                xmllog.endEntity().flush();

            xmllog.endEntity().flush();

            } catch (IOException ex) {
                logger.severe(ex.getMessage());
            }
    }


	@Override
	protected void Process(){
			xmllog.writeEntity("InitialEvaluation");
				Evaluate(param,4);
			xmllog.endEntity().flush();

			xmllog.writeEntity("Optimization").writeAttribute("Method",sMethod );
				long duration=System.currentTimeMillis();

                nEvaluations=0;
                if(nOpts>0)
                    if(sMethod.contentEquals("GA"))
                        FitUsingGeneticAlgorithm();
                    else
                        FitUsingApacheLevenbergMarquardt();

				duration=(System.currentTimeMillis()-duration);

                 xmllog.writeEntity("ConvergenceInfo");
                    xmllog.writeAttribute("TotalRMS", Double.toString(finalRMS));
                    //xmllog.writeAttribute("MaxInterations", Integer.toString(fit.getMaxIterations()));
                    //xmllog.writeAttribute("MaxEvaluations", Integer.toString(fit.getMaxEvaluations()));
                    xmllog.writeAttribute("Iterations", Integer.toString(nIterations));
                    xmllog.writeAttribute("Evaluations", Integer.toString(nEvaluations));
                    xmllog.writeAttribute("Duration", Double.toString(duration/1000.0));
                xmllog.endEntity().flush();
			xmllog.endEntity().flush();


            if(finalParam!=null){
                xmllog.writeEntity("FinalEvaluation");
                    Evaluate(finalParam, 4);
                xmllog.endEntity().flush();

                xmllog.writeEntity("FinalParameters");
                    WriteParam(finalParam);
                xmllog.endEntity().flush();
            }
	}


    @Override
	public void ParseArguments(String[] args_) {
        this.args=args_;
		parser = new CmdLineParser(this);
		try {
			parser.parseArgument(args);
			if(isHelp){
				parser.printUsage(System.out);
				return ;
			}

			for(int i=0;i<50;i++){
				//boolean bSet=false;
				String opt=String.format("-i%d",i+1);
					for(int j=0;j<args.length-1;j++)
						if(args[j].contentEquals(opt)){
							//bSet=true;
							datafile.add(args[j+1]);
						}

//				if(bSet){
//					double t=1;
//					opt=String.format("-w%d",i+1);
//
//					for(int j=0;j<args.length-1;j++)
//						if(args[j].contentEquals(opt))
//							t=Double.parseDouble(args[j+1]);
//
//					dataWeight.add(t);
//				}
			}
		} catch (CmdLineException ex) {
			logger.severe(ex.getMessage());
		} 
	}


    //===============================
    final int MASTER=0;

    int rank = MASTER;
    int ncpu = 2;
    Potential[] potpool;


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
	String[] label_inp; //!< label of parameter

	Vector<Cluster> molRef= new Vector<Cluster>();//!< reference clusters
	Vector<Cluster> data=new Vector<Cluster>();//!< data/fitting clusters
	double[] weight;//!< weighted numbers of each clusters

	double[] target; //!< contain target data, in that case is binding energy +rms

	Vector<String> datafile=new Vector<String>();

    int nIterations=0;
    int nEvaluations=0;
    double[] finalParam;
    double finalRMS=0;

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

    protected void ReadReference() throws FileNotFoundException{
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
    }

    protected void ReadData() throws FileNotFoundException{
        ArrayList<ArrayList<Double> > arrayData=new ArrayList<ArrayList<Double>>();
        //reading the data
        for (int i = 0; i < datafile.size();i++) {
            ArrayList<Double> tmpData= new ArrayList<Double>();
			Scanner scanData = new Scanner(new File(datafile.get(i)));
			int count=0;

			while (scanData.hasNext()) {
				Cluster tmpMol = new Cluster();
				tmpMol.Read(scanData, "xyz");
				if(tmpMol.getNAtoms()==0) break;
				tmpMol.CorrectOrder();

				data.add(tmpMol);

                tmpData.add(tmpMol.getEnergy());
				count++;
			}
            arrayData.add(tmpData);

			xmllog.writeEntity("Data");
			xmllog.writeAttribute("File", datafile.get(i));
			xmllog.writeAttribute("Size", Integer.toString(count));
			xmllog.endEntity().flush();
			scanData.close();
		}
        
        //calculate target energies
        double[] ener_ref={molRef.get(0).getEnergy(),molRef.get(1).getEnergy(),molRef.get(2).getEnergy()};
		target=new double[data.size()];
		for (int i=0; i< data.size(); i++){
			Cluster m=data.get(i);
			target[i]= CalcBindingEnergy(m,m.getEnergy(),ener_ref);
			//target[i]*=1000;
		}
        
        
        //calculate weighted numbers
        weight = new double[data.size()];
        if(bWeighted){
            int i=0;
            for (ArrayList<Double> d : arrayData) {
                //finding the smallest ones
                double minE=Double.MAX_VALUE;
                double maxE=-Double.MAX_VALUE;
                for(Double k : d){
                    minE=Math.min(minE,k);
                    maxE=Math.max(maxE,k);
                }

                //System.out.println(" Min = " + minE + " max = "+maxE);
                double c=maxE+0.1*(maxE-minE)/d.size();
                double normalizationFactor=0;
               
                for(Double k : d){
                   normalizationFactor+=(c-k);
                }
                

                //scale weighted numbers
                for(Double k : d){
                    weight[i]=(c-k)/normalizationFactor*d.size()/weight.length;
                    i++;
                }

            }
        }else //if not weighted
            for(int i=0;i<weight.length;i++){
                weight[i]=1.0/weight.length;
            }
    }   

    protected void ReadBound() throws FileNotFoundException, IOException{		
        param_inp=pot.getParam().clone();
        ub_inp=new double[param_inp.length];//!< input upper and lower bounds
        lb_inp=new double[param_inp.length];
        fixed=new int[param_inp.length]; //!< determine wheter parameter are fixed
        label_inp=new String[param_inp.length];

		int nfixed = 0;

        File file=new File(sFileBound);
        
        if(file.canRead()){//if bound file exists, read it
            Scanner fin = new Scanner(file);
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
            fin.close();
        }else{ //if not, predict from input value
            double fluc=0.2;
            if(!sFileBound.isEmpty()) fluc=Double.parseDouble(sFileBound);
            for(int i=0;i<param_inp.length;i++){
                lb_inp[i]=param_inp[i]-fluc*Math.abs(param_inp[i]);
                ub_inp[i]=param_inp[i]+fluc*Math.abs(param_inp[i]);
                fixed[i]=0;
                label_inp[i] = "";
            }

        }

		nparam = param_inp.length - nfixed;		

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
                        FileWriter fileOut = new FileWriter(new File(sFileOut));
                        pot.setParam(result);
                        pot.writeParam(fileOut);
                        fileOut.close();
                   } catch (IOException ex) {
                       logger.log(Level.SEVERE, "Cannot write to {0}", sFileOut);
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

//    class ThreadEvaluate extends Thread{
//        private double[] calcEner;
//        private int iFrom=0;
//        private int iTo=0;
//        private Potential pot;
//
//        public void setRange(int iFrom,int iTo){
//            this.iFrom=iFrom;
//            this.iTo=iTo;
//        }
//
//        public void setCalcEner(double[] a){
//            this.calcEner=a;
//        }
//
//        public void setPotential(Potential pot){
//            this.pot=pot;
//        }
//
//
//        @Override
//        public void run() {
//			for (int i=iFrom; i<iTo; i++){
//				Cluster m=data.get(i);
//				//calculate binding energy
//				this.pot.setCluster(m);
//				this.calcEner[i] = this.pot.getEnergy(m.getCoords());
//                //System.out.println(Thread.currentThread().getName()+" id = "+i+" e="+this.calcEner[i]);//
//			}
//        }
//
//    }

    double Evaluate(double[] p,int verbose){
		double[] y=new double[data.size()];
		return Evaluate(p,y,verbose);
	}

    double Evaluate(double[] p,double[] y,int verbose){
        double[] alpha=new double[param_inp.length];
        ConvertParam(p,alpha);
        pot.setParam(alpha);

        if(verbose>=2){
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

        if (verbose>=3){
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

        //calculate y here
        double[] calcEner=new double[data.size()];
        //calculate y here
//        ThreadEvaluate[] threads=new ThreadEvaluate[ncpu];
//
//        for(int i=0;i<threads.length;i++){
//            System.out.println("Iniatialize " + i);
//            threads[i] = new ThreadEvaluate();
//            threads[i].setCalcEner(calcEner);
//            threads[i].setPotential(potpool[i]);
//            int nCalcs=(int)(data.size()/ncpu)+1;
//            if(i!=threads.length-1) threads[i].setRange(nCalcs*i,nCalcs*(i+1));
//            else threads[i].setRange(nCalcs*i,data.size());
//        }
//
//        for(int i=0;i<threads.length;i++){
//            threads[i].start();
//        }
//
//        for (int i = 0; i < threads.length; i++) {
//           System.out.println("Join " + i);
//           try {
//              threads[i].join();
//           } catch (InterruptedException ignore) {
//           }
//       }


        for (int i=0; i<data.size(); i++){
            Cluster m=data.get(i);
            //calculate binding energy
            pot.setCluster(m);
            calcEner[i] = pot.getEnergy(m.getCoords());
            double ener=calcEner[i];
            y[i] = CalcBindingEnergy(m,ener,ener_shift);

            double dE = Math.abs(y[i] - target[i]);
            rms+=dE*dE*weight[i];

            if(verbose>=4){
                xmllog.writeEntity("Eval");
                xmllog.writeAttribute("id", Integer.toString(i));
                xmllog.writeAttribute("nAtoms", Integer.toString(m.getNAtoms()));
                xmllog.writeAttribute("Weight", Double.toString(weight[i]));
                xmllog.writeAttribute("ObservedBE", Double.toString(target[i]));
                xmllog.writeAttribute("CalcE", Double.toString(ener));
                xmllog.writeAttribute("CalcBE", Double.toString(y[i]));
                xmllog.writeAttribute("DeltaE", Double.toString(dE));
                xmllog.endEntity().flush();
            }
        }

        rms=Math.sqrt(rms);

        if (verbose>=1){
             xmllog.writeEntity("Step").writeAttribute("id", Integer.toString(nEvaluations));
					xmllog.writeAttribute("RMS", Double.toString(rms));
		     xmllog.endEntity().flush();
        }

        nEvaluations++;
        return rms;
    }

	double Evaluate_old(double[] p,double[] y,int verbose){
        double[] alpha=new double[param_inp.length];
        ConvertParam(p,alpha);
        pot.setParam(alpha);

        if(verbose>=2){
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

        if (verbose>=3){
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

        //omp for private(m) reduction(+:rms)
        for (int i=0; i<data.size(); i++){
            Cluster m=data.get(i);
            //calculate binding energy
            pot.setCluster(m);
            double ener = pot.getEnergy(m.getCoords());

            y[i] = CalcBindingEnergy(m,ener,ener_shift);
            //y[i]*=1000;

            double dE = Math.abs(y[i] - target[i])* Math.sqrt(weight[i]);
            rms+=dE*dE;

            if(verbose>=4){
                xmllog.writeEntity("Eval");
                xmllog.writeAttribute("id", Integer.toString(i));
                xmllog.writeAttribute("nAtoms", Integer.toString(m.getNAtoms()));
                xmllog.writeAttribute("Weight", Double.toString(weight[i]));
                xmllog.writeAttribute("ObservedBE", Double.toString(target[i]));
                xmllog.writeAttribute("CalcE", Double.toString(ener));
                xmllog.writeAttribute("CalcBE", Double.toString(y[i]));
                xmllog.writeAttribute("DeltaE", Double.toString(dE));
                xmllog.endEntity().flush();
            }
        }

        rms=Math.sqrt(rms/data.size());

        //convert unit
        rms=Potential.ConvertUnit(rms,"Hartree","kcal/mol");

        if (verbose>=1){
             xmllog.writeEntity("Step").writeAttribute("id", Integer.toString(nEvaluations));
					xmllog.writeAttribute("RMS", Double.toString(rms));
		     xmllog.endEntity().flush();
//                xmllog.writeEntity("TotalRMS");
//                xmllog.writeNormalText(Double.toString(rms));
//                xmllog.endEntity().flush();
        }

        nEvaluations++;
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
						a[i]+=delta;//disturb a bit
						double[] y2=APACHEObjectiveFunction.this.value(a);
						a[i]=tmp;//restore to original value
						for(int j=0;j<data.size();j++)
							result[j][i]=(y2[j]-y0[j])/(delta);
					}

					//two points jacobian
//                    double[][] result=new double[1][a.length];
//                    double ratio=1e-3;
//					for(int i=0;i<a.length;i++){
//						double delta=Math.abs(a[i])*ratio;
//						a[i]-=delta;
//						double[] y1=APACHEObjectiveFunction.this.value(a);
//						a[i]+=2*delta;
//						double[] y2=APACHEObjectiveFunction.this.value(a);
//						a[i]-=delta;
//						result[0][i]=(y2[0]-y1[0])/(2*delta);
//					}


					return result;
				}
			}

		Jacobian jacobian=new Jacobian();

        double bestRMS=1e9;
        double[] bestParam;

        public double getBestRMS(){
            return bestRMS;
        }

        public double[] getBestParam(){
            return bestParam;
        }

		@Override
		public MultivariateMatrixFunction jacobian() {
			return jacobian;
		}

		@Override
		public double[] value(double[] a) throws FunctionEvaluationException, IllegalArgumentException {
			double[] y=new double[data.size()];
            double rms;
            if(nEvaluations%100==0)
                rms=Evaluate(a, y, 1);
            else
                rms=Evaluate(a, y, 0);

            if(rms<bestRMS){
                bestParam=a;
                bestRMS=rms;
            }
            //double[] y= new double[1];
           // y[0]=Evaluate(a,0);
			//System.err.printf(" valuesize=%d\n", y.length);
			return y;
		}

	}

	protected void FitUsingApacheLevenbergMarquardt() {
        try {
            APACHEObjectiveFunction obj = new APACHEObjectiveFunction();
            LevenbergMarquardtOptimizer fit = new LevenbergMarquardtOptimizer();
            //System.err.printf("targetsize=%d \n",target.length);
            fit.setCostRelativeTolerance(1e-14);
            fit.setOrthoTolerance(1e-14);
            fit.setParRelativeTolerance(1e-14);
            fit.setMaxIterations(nOpts);
            fit.setMaxEvaluations(100 * nOpts);
            //fitting here
            VectorialPointValuePair r = fit.optimize(obj, target, weight, param);
            finalParam = obj.getBestParam(); //
            finalRMS = obj.getBestRMS();
            nIterations = fit.getIterations();
        } catch (FunctionEvaluationException ex) {
            Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
        } catch (OptimizationException ex) {
            Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IllegalArgumentException ex) {
            Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
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

    protected void FitUsingGeneticAlgorithm(){
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
                  IChromosome bestSolutionSoFar = population.getFittestChromosome();
                  double fitness = bestSolutionSoFar.getFitnessValue();

                  xmllog.writeEntity("Step").writeAttribute("id", Integer.toString(i));
					xmllog.writeAttribute("RMS", Double.toString(1.0/fitness));
				  xmllog.endEntity().flush();

                  //System.out.println(" size ="+ fittest.size());
                  double[] currentBestParam=new double[bestSolutionSoFar.size()];
                  for (int k = 0; k < bestSolutionSoFar.size(); k++) {
                    DoubleGene gene = (DoubleGene)bestSolutionSoFar.getGene(k);
                    currentBestParam[k]=gene.doubleValue();
                   // System.out.println(" param ="+ currentParam[k]);
                  }
                  WriteParam(currentBestParam);
                }
            }
            

           IChromosome bestSolutionSoFar = population.getFittestChromosome();
           finalRMS=1.0/bestSolutionSoFar.getFitnessValue();

           finalParam=new double[bestSolutionSoFar.size()];
           for (int i = 0; i < bestSolutionSoFar.size(); i++) {
                    DoubleGene gene = (DoubleGene)bestSolutionSoFar.getGene(i);
                    finalParam[i]=gene.doubleValue();
            }
        } catch (InvalidConfigurationException ex) {
            Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
        }
              
    }
}

