/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import jaolho.data.lma.LMA;
import jaolho.data.lma.LMAFunction;
import java.io.*;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.DifferentiableMultivariateVectorialFunction;
import org.apache.commons.math.analysis.MultivariateMatrixFunction;
import org.kohsuke.args4j.CmdLineException;
import org.kohsuke.args4j.CmdLineParser;
import org.kohsuke.args4j.Option;



/**
 *
 * @author nqc
 */
public class TaskFitPotential extends TaskCalculate {	
	public class LMAObjectiveFunction extends LMAFunction {
		@Override
		public double getY(double x, double[] a) {
			double[] alpha=new double[nparam_inp];
			ConvertParam(a,alpha);
			pot.setParam(alpha);

			//calculate
			double[] ener_shift=new double[molRef.size()];
			
			for(int i=0; i<molRef.size(); i++){
				Cluster m=molRef.get(i);
				pot.setCluster(m);
				ener_shift[i] = pot.getEnergy(m.getCoords());
			}

			int i=(int)x;
			Cluster m=data.get(i);
			//calculate binding energy
			double ener = pot.getEnergy(m.getCoords());

			int index=MolExtra.indexForBE(m, molRef);
			double y = ener - ener_shift[0]*(m.getNonHydrogenNum()-1) - ener_shift[index];
			return y;
		}

		@Override
		public double getPartialDerivate(double x, double[] a, int parameterIndex) {
			double[] t=a.clone();
			double delta=1e-5;
			t[parameterIndex]-=delta;
			double y1=getY(x,t);
			t[parameterIndex]+=2*delta;
			double y2=getY(x,t);

			return (y2-y1)/(2*delta);
		}
	}

	public class APACHEObjectiveFunction implements DifferentiableMultivariateVectorialFunction {//used with apache commom optimization

		@Override
		public MultivariateMatrixFunction jacobian() {
			throw new UnsupportedOperationException("Not supported yet.");
		}

		@Override
		public double[] value(double[] a) throws FunctionEvaluationException, IllegalArgumentException {
			double[] alpha=new double[nparam_inp];
			ConvertParam(a,alpha);
			pot.setParam(alpha);

			//calculate
			double[] ener_shift=new double[molRef.size()];

			for(int i=0; i<molRef.size(); i++){
				Cluster m=molRef.get(i);
				pot.setCluster(m);
				ener_shift[i] = pot.getEnergy(m.getCoords());
			}

			double[] y=new double[data.size()];
			for (int i=0; i<data.size(); i++){
				Cluster m=data.get(i);
				//calculate binding energy
				pot.setCluster(m);
				double ener = pot.getEnergy(m.getCoords());

				int index=MolExtra.indexForBE(m, molRef);
				y[i] = ener - ener_shift[0]*(m.getNonHydrogenNum()-1) - ener_shift[index];
			}
			return y;
		}

	}

	@Option(name="-ref",usage="to calculate binding energies. It MUST be in F_XYZ format and in the following order: neutral, protonated, deprotonated. ",metaVar="FILE")
	String sFileRef="";

	
		
		final int cParamMax=100;

		int nparam; //!< number of actual parameters
		double[] param; //!< actual parameters
		double[] ub;//!< actual bounds, all fixed parameters are removed
		double[] lb;
		String[] label; //!< label of parameter

		int nparam_inp; //!< number of parameters (in total)
		double[] param_inp=new double[cParamMax]; //!< starting value of parameters
		double[] ub_inp=new double[cParamMax];//!< input upper and lower bounds
		double[] lb_inp=new double[cParamMax];
		int[] fixed=new int[cParamMax]; //!< determine wheter parameter are fixed
		int nfixed; //!< number of fixed parameters
		String[] label_inp=new String[cParamMax]; //!< label of parameter

		Vector<Cluster> molRef= new Vector<Cluster>();
		Vector<Cluster> data=new Vector<Cluster>();
		Vector<Double> weight=new Vector<Double>();

		double[][] observation; //!< consist of observation data, in that case is binding energy +rms

		Vector<String> datafile=new Vector<String>();
		Vector<Double> dataWeight=new Vector<Double>();

	LMAObjectiveFunction fit=new LMAObjectiveFunction();

	
	void ProcessData(){
		observation=new double[2][data.size()];
		for (int i=0; i< data.size(); i++){
			Cluster m=data.get(i);
			int index=MolExtra.indexForBE(m,molRef);
			observation[0][i]=i;
			observation[1][i]= m.getEnergy() - (m.getNonHydrogenNum()-1)*molRef.get(0).getEnergy() - molRef.get(index).getEnergy();
		}
	}

	double Evaluate(double[] p,int debug){
			double[] alpha=new double[nparam_inp];
			ConvertParam(p,alpha);
			pot.setParam(alpha);


			if(debug>=1){
				try {
					xmllog.writeEntity("Alpha");
					String s = "";
					for (int i = 0; i < alpha.length; i++) {
						s += String.format("%f ", alpha[i]);
					}
					xmllog.writeNormalText(s);
					xmllog.endEntity().flush();
				} catch (IOException ex) {
					Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
				}

			}
			//calculate
			double[] ener_shift=new double[molRef.size()];
			double rms=0;

			for(int i=0; i<molRef.size(); i++){
				Cluster m=molRef.get(i);
				pot.setCluster(m);
				ener_shift[i] = pot.getEnergy(m.getCoords());
			}

			if (debug>=2){
				try {
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
				} catch (IOException ex) {
					Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
				}
			}

			double[] y=new double[data.size()];
			//cout <<"\n\n@$    inputs# nAtoms Weight     Obsr.Data        Calc.Data    deltaE" << endl;
			for (int i=0; i<data.size(); i++){
				Cluster m=data.get(i);
				//calculate binding energy
				pot.setCluster(m);
				double ener = pot.getEnergy(m.getCoords());

				int index=MolExtra.indexForBE(m, molRef);
				y[i] = ener - ener_shift[0]*(m.getNonHydrogenNum()-1) - ener_shift[index];

				double dE = Math.abs(y[i] - observation[1][i])* Math.sqrt(weight.get(i));
				rms+=dE*dE;

				if(debug>=3){
					try {
						xmllog.writeEntity("Eval");
						xmllog.writeAttribute("id", Integer.toString(i));
						xmllog.writeAttribute("nAtoms", Integer.toString(m.getNAtoms()));
						xmllog.writeAttribute("Weight", Double.toString(weight.get(i)));
						xmllog.writeAttribute("ObserveBE", Double.toString(observation[1][i]));
						xmllog.writeAttribute("CalcEnergy", Double.toString(ener));
						xmllog.writeAttribute("CalcBE", Double.toString(y[i]));
						xmllog.writeAttribute("DeltaE", Double.toString(dE));
						xmllog.endEntity().flush();
					} catch (IOException ex) {
						Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
					}
				}
			}

			rms=Math.sqrt(rms/data.size());

			if (debug>=1){
				try {
					xmllog.writeEntity("TotalRMS");
					xmllog.writeNormalText(Double.toString(rms));
					xmllog.endEntity().flush();
				} catch (IOException ex) {
					Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
				}
			}

			return rms;
		}

	void ConvertParam(double[] actual,double[] p){
		int c=0;
		for(int i=0; i<nparam_inp; i++)
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
				//parser.printUsage(System.out);
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
			Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
		} catch (IOException ex) {
			Logger.getLogger(Task.class.getName()).log(Level.SEVERE, null, ex);
		}

	}


	@Override
	public String getName() {
		return "FitPotential";
	}

	@Override
	protected void Initialize(){
	try {
		xmllog.writeAttribute("InputFile", sFileIn).writeAttribute("FormatInput", sFormatIn);
		xmllog.writeAttribute("Verbose", Integer.toString(verbose));

		xmllog.writeEntity("Initialize");
		pot = MolExtra.SetupPotential(sPotential, sFileParam, sUnit,sMethod);
		xmllog.writeNormalText(pot.Info(1));

		xmllog.writeEntity("InitialParameters");
		Scanner fin = new Scanner(new File(sFileParam));
		nparam_inp = nfixed = 0;
		int i = 0;
		while (fin.hasNext()) {
			String line = fin.nextLine();
			if (line.isEmpty()) {
				break;
			}
			StringTokenizer tokenizer;
			tokenizer = new StringTokenizer(line, "\t ,;"); //read no of atoms
			String info;
			info = tokenizer.nextToken();
			param_inp[i] = Double.parseDouble(info);
			info = tokenizer.nextToken();
			lb_inp[i] = Double.parseDouble(info);
			info = tokenizer.nextToken();
			ub_inp[i] = Double.parseDouble(info);
			info = tokenizer.nextToken();
			fixed[i] = Integer.parseInt(info);
			label_inp[i] = "";
			while (tokenizer.hasMoreTokens()) {
				label_inp[i] = tokenizer.nextToken() + " ";
			}
			//sscanf(line.c_str(),"%lf %lf %lf %d %line", &param_inp[nparam_inp], &lb_inp[nparam_inp], &ub_inp[nparam_inp], &fixed[nparam_inp]);
			nfixed += fixed[i];
			xmllog.writeEntity("Param").writeAttribute("id", Integer.toString(i));
			xmllog.writeAttribute("Value", Double.toString(param_inp[i]));
			xmllog.writeAttribute("Lower", Double.toString(lb_inp[i]));
			xmllog.writeAttribute("Upper", Double.toString(ub_inp[i]));
			xmllog.writeAttribute("Fixed", Integer.toString(fixed[i]));
			xmllog.writeAttribute("Label", label_inp[i]);
			xmllog.endEntity().flush();
			i++;
		}
		nparam_inp = i;
		nparam = nparam_inp - nfixed;
		fin.close();
		xmllog.endEntity();

		xmllog.writeEntity("SetActualParameters");
		xmllog.writeAttribute("NParam", Integer.toString(nparam));
		xmllog.writeAttribute("NFixed", Integer.toString(nfixed));
		param = new double[nparam];
		ub = new double[nparam];
		lb = new double[nparam];
		label = new String[nparam];
		int t = 0;
		for (i = 0; i < nparam_inp; i++) {
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
		for (i = 0; i < datafile.size();i++) {
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

		xmllog.endEntity().flush();


	} catch (IOException ex) {
		Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
	}

}
	
	@Override
	protected void Process(){
		try {
			xmllog.writeEntity("InitialEvaluation");
				Evaluate(param, 4);
			xmllog.endEntity().flush();

			xmllog.writeEntity("Optimization");
				//for(int i=0;i<10;i++)
					//System.out.printf("%d %f\n",i,fit.getY(i,param));
					
				//MultivariateRealOptimizer opt = new NelderMead();

				//RealPointValuePair r = null; //=new RealPointValuePair();

				long duration=System.currentTimeMillis();
				LMA lma = new LMA(
					fit,
					param,
					observation
				);
				lma.fit();
				//opt.setMaxIterations(nOpts);
				//r = opt.optimize(this, GoalType.MINIMIZE,param);
				duration=(System.currentTimeMillis()-duration);
				//double newEnergy = r.getValue();


				
				xmllog.writeAttribute("TotalRMS", Double.toString(Math.sqrt(lma.getRelativeChi2())));
				//xmllog.writeAttribute("NEvals", Integer.toString(opt.getEvaluations()));
				xmllog.writeAttribute("Duration", Double.toString(duration/1000.0));
			xmllog.endEntity().flush();

			xmllog.writeEntity("FinalEvaluation");
				//Evaluate(param, 4);
			xmllog.endEntity().flush();
			
			xmllog.writeEntity("FinalParameters");				
				if(!sFileOut.isEmpty()){
					fileOut=new FileWriter(new File(sFileOut));
				}

//				double[] result=new double[nparam_inp];
//				ConvertParam(r.getPoint(),result);
//
//				for(int i=0;i<result.length;i++){
//					if(fileOut!=null)
//						fileOut.write(String.format("%f\t%f\t%f\t%d\t%s\n",result[i],lb_inp[i],ub_inp[i],fixed[i],label_inp[i]));
//
//					xmllog.writeEntity("Param").writeAttribute("id", Integer.toString(i));
//					xmllog.writeAttribute("Value", Double.toString(result[i]));
//					xmllog.writeAttribute("Lower", Double.toString(lb_inp[i]));
//					xmllog.writeAttribute("Upper", Double.toString(ub_inp[i]));
//					xmllog.writeAttribute("Fixed", Integer.toString(fixed[i]));
//					xmllog.writeAttribute("Label", label_inp[i]);
//					xmllog.endEntity().flush();
//				}
				if(fileOut!=null)
					fileOut.close();
			xmllog.endEntity().flush();

		} catch (IllegalArgumentException ex) {
			Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
		} catch (IOException ex) {
			Logger.getLogger(TaskFitPotential.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	
}
