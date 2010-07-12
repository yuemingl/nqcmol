/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.potential;

import java.io.IOException;
import java.io.Writer;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.Cluster;
import nqcmol.tools.MTools;
import nqcmol.tools.XmlWriter;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.*;
import org.apache.commons.math.optimization.*;
import org.apache.commons.math.optimization.general.*;

/**
 *
 * @author nqc
 */
public class Potential {
	
	public Potential() {
		isValidSetup=false;
	}

	public Potential(Cluster cluster) {
		setCluster(cluster);
		isValidSetup=false;
	}


	// general variables
	protected boolean isValidSetup=false;

	/**
	 * Get the value of isValidSetup
	 *
	 * @return the value of isValidSetup
	 */
	public boolean isValidSetup() {
		return isValidSetup;
	}


   /**
	* Molecule to be evaluated or minimized
	*/
	protected Cluster cluster;

	/**
	 * Get the value of cluster
	 *
	 * @return the value of cluster
	 */
	public Cluster getCluster() {
		return cluster;
	}

	/**
	 * Set cluster to the internal one. Note that, some potential CAN change cluster to suitable format.
	 * @param cluster new value of cluster
	 */
	public void setCluster(Cluster cluster_) {
		this.cluster = cluster_;
		nCoords=cluster.getNcoords();
	}


	int nCoords=0;
	
	public int getNCoords(){
		return cluster.getNcoords();
	}
    // minimization variables

	protected String nativeUnit = "au";

	/**
	 * Get the value of nativeUnit
	 *
	 * @return the value of nativeUnit
	 */
	public String getNativeUnit() {
		return nativeUnit;
	}


    /**
	 * The unit (kcal/mol, kJ/mol, ...) in which the energy is expressed
     */
	protected String unit = "au";

	/**
	 * Get the value of unit
	 *
	 * @return the value of unit
	 */
	public String getUnit() {
		return unit;
	}

	/**
	 * Set unit
	 *
	 * @return the value of unit
	 */
	public void setUnit(String s) {
		if(!s.isEmpty())	unit=s;
	}

    /**
	 * Does this force field have analytical gradients.
     * @return True if all analytical gradients are implemented.
     */
    public boolean HasAnalyticalGradients() { return false; }

	/**
	 * Does this force field have analytical hessian (2nd derivatives).
     * @return True if all analytical gradients are implemented.
     */
    public boolean HasAnalyticalHessian() { return false; }

	protected double[] param=null;

	/**
	 * Get the value of param
	 *
	 * @return the value of param
	 */
	public double[] getParam() {
		return param;
	}

	/**
	 * Set the value of param
	 *
	 * @param param new value of param
	 */
	public void setParam(double[] param) {
		this.param = param.clone();
		//System.arraycopy(param,0, this.param,0,param.length);
	}

	/**
	 * Set the value of param
	 *
	 * @param param new value of param
	 */
	public void setParam(String filename) {
	}


	/**
	 * Get the value of equation
	 *
	 * @return the value of equation
	 */
	public String getEquation() {
		String equation = "Not Specified";
		return equation;
	}

    /////////////////////////////////////////////////////////////////////////
    // Public methods for evaluation, gradients and optimize              //
    /////////////////////////////////////////////////////////////////////////

	/**
	 * setup new cluster  and then calculate energy
	 * @return calculated energy
	 */
	public double getEnergy(Cluster cluster_) {
		setCluster(cluster_);
		double ener=Energy_(cluster.getCoords());
		cluster.setEnergy(ener);
		return ener;
	}

	/**
	 * return current energy of  cluster
	 * @param isEnergyUpdate: if true, it will call Energy function to update energy
	 * @return
	 */
	public double getEnergy(boolean isEnergyUpdate){
		if(isEnergyUpdate){
			double ener=Energy_(cluster.getCoords());
			cluster.setEnergy(ener);
		}
		return cluster.getEnergy();
	}

	/**
	 * return energy. Does not update neither coordinates and energy of cluster
	 * @return energy
	 */
	public double getEnergy(final double[] coords_){
		return Energy_(coords_);
	}

	public double[][] getGradient(Cluster cluster_){
		setCluster(cluster_);
		Gradient_(cluster.getCoords(),cluster.getGradient());
		return cluster.getHessian();
	}

	/**
	 * assume that cluster is set already, update cluster coordinates and calculate gradients
	 * @param coords_ new coordination
	 * @return reference of cluster gradient;
	 */
	public double[] getGradient(final double[] coords_){
		cluster.setCoords(coords_);
		Gradient_(cluster.getCoords(),cluster.getGradient());
		return cluster.getGradient();
	}

	/**
	 * return gradients. Does not update neither coordinates and gradients of cluster
	 */
	public void getGradient(final double[] coords_,double[] gradient_){
		Gradient_(coords_,gradient_);
	}

	/**
	 * return current energy of  cluster
	 * @param isGradUpdate: if true, it will call Gradient to update gradient
	 * @return reference of cluster gradient
	 */
	public double[] getGradient(boolean isGradUpdate){
		if(isGradUpdate){
			Gradient_(cluster.getCoords(),cluster.getGradient());
		}
		return cluster.getGradient();
	}

	public double[][] getHessian(boolean isHessianUpdate){
		if(isHessianUpdate) Hessian_(cluster.getCoords(),cluster.getHessian());
		return cluster.getHessian();
	}

	public double[][] getHessian(Cluster cluster_){
		setCluster(cluster_);
		Hessian_(cluster.getCoords(),cluster.getHessian());
		return cluster.getHessian();
	}

	public double[][] getNumericalHessian(Cluster cluster_,double delta){
		setCluster(cluster_);
		NumericalHessian_(cluster.getCoords(),cluster.getHessian(),delta);
		return cluster.getHessian();
	}

	public double getDeltaEnergy(final double[] coords_,final boolean[] isUpdate){//!< require setup mol first
		return DeltaEnergy_(coords_,isUpdate) ;
	}

	/**
	 * Perform optimization on the cluster set previouly by setCluster
	 * @return true if succesfully converged
	 */
	public boolean Optimize(){//!< require setup mol first
		//setting the variable
		nEvals=0;

//		if(debug>=2){
//			cout<<"Dimension = "<<Dim_();
//			cout<<" Energy="<<evaluate_(x)<<" Method="<<optMethod<<endl;
//		}

		double[] x=cluster.getCoords();
		optimizedE=cluster.getEnergy();

		if(optMethod.contentEquals("DFPMIN")){//quasi newton
			iConverged=DFPmin(x);
			
		}else if(optMethod.contentEquals("CG")){//conjugate gradient
			DifferentiableMultivariateRealOptimizer opt=new NonLinearConjugateGradientOptimizer(ConjugateGradientFormula.FLETCHER_REEVES);
			DiffFunc func=new DiffFunc();
			opt.setMaxEvaluations(nMaxEvals);
			opt.setConvergenceChecker(new Convergence());

			RealPointValuePair r = null;
			try {
				r=opt.optimize(func, GoalType.MINIMIZE, x);
			} catch (FunctionEvaluationException ex) {
				//Logger.getLogger(Potential.class.getName()).log(Level.SEVERE, null, ex);
			} catch (OptimizationException ex) {
				//Logger.getLogger(Potential.class.getName()).log(Level.SEVERE, null, ex);
			} catch (IllegalArgumentException ex) {
				Logger.getLogger(Potential.class.getName()).log(Level.SEVERE, null, ex);
			}

			optimizedE=r.getValue();

			x=r.getPoint();
            RmsGrad=this.getRMSGrad(x);
		}else if(optMethod.contentEquals("NETMEAD")){//conjugate gradient

		}


		cluster.setCoords(x);
		cluster.setEnergy(optimizedE);
        cluster.setRmsGrad(RmsGrad);

		return (iConverged!=0);
	}

	/**
	 * Perform optimization
	 * @param cluster_ : input cluster, will be optimized
	 * @return true if succesfully converged
	 */
	public boolean Optimize(Cluster cluster_){
		this.setCluster(cluster_);
		return Optimize();
	}


    /**
	 * Validate the analytical gradients by comparing them to numerical ones.
	 * @param deltaX small displacment for numerical gradients.
	 * @return a string in XML format.
    */
    public String ValidateGradient(double deltaX){
		String sVerbose="";

		try {
			double maxerr = 100;
			//XmlWriter.test1();
			Writer writer = new java.io.StringWriter();
			XmlWriter xmlwriter = new XmlWriter(writer);
			xmlwriter.writeEntity("ValidateGradient");
			xmlwriter.writeAttribute("Equation",getEquation());
			xmlwriter.writeAttribute("DeltaX",Double.toString(deltaX));
			xmlwriter.writeEntity("Note");
			xmlwriter.writeText(" Relative Error is percentage");
			xmlwriter.endEntity();

//			if (verbose > 0) {
//				System.out.printf("\nValidate Gradients of %s\n\n", getEquation());
//				System.out.print("ATOM IDX                    NUMERICAL GRADIENT                        ANALYTICAL GRADIENT                     REL. ERROR ()   \n");
//				System.out.print("-------------------------------------------------------------------------------------------------------------------------------\n");
//			}
			if (HasAnalyticalGradients()) {
				Gradient_(cluster.getCoords(), cluster.getGradient());
			}
			double[] anagrad = cluster.getGradient();
			double[] numgrad = new double[cluster.getNcoords()];
			NumericalGradient_(cluster.getCoords(), numgrad, deltaX);
			for (int i = 0; i < cluster.getNAtoms(); i++) {
				double[] err = new double[3];
				for (int k = 0; k < 3; k++) {
					if (numgrad[i * 3 + k] != 0) {
						err[k] = 100 * Math.abs((anagrad[i * 3 + k] - numgrad[i * 3 + k]) / numgrad[i * 3 + k]);
						maxerr = Math.max(err[k],maxerr);
					} else {
						if (anagrad[i * 3 + k] != 0)
                            err[k] = 100.0;
                        else err[k]=0;
					}
				}

				//String info=String.format("%2d       (%12.6f, %12.6f, %12.6f)  (%12.6f, %12.6f, %12.6f)  (%6.3f, %6.3f, %6.3f)\n", i, numgrad[i * 3], numgrad[i * 3 + 1], numgrad[i * 3 + 2], anagrad[i * 3], anagrad[i * 3 + 1], anagrad[i * 3 + 2], err[0], err[1], err[2]);

				xmlwriter.writeEntity("Test");
					xmlwriter.writeAttribute("AtomID",Integer.toString(i));
					xmlwriter.writeEntity("NumGrad");
						xmlwriter.writeAttribute("x",Double.toString(numgrad[i*3+0]));
						xmlwriter.writeAttribute("y",Double.toString(numgrad[i*3+1]));
						xmlwriter.writeAttribute("z",Double.toString(numgrad[i*3+2]));
					xmlwriter.endEntity();
					xmlwriter.writeEntity("AnaGrad");
						xmlwriter.writeAttribute("x",Double.toString(anagrad[i*3+0]));
						xmlwriter.writeAttribute("y",Double.toString(anagrad[i*3+1]));
						xmlwriter.writeAttribute("z",Double.toString(anagrad[i*3+2]));
					xmlwriter.endEntity();
					xmlwriter.writeEntity("RelError");
						xmlwriter.writeAttribute("x",Double.toString(err[0]));
						xmlwriter.writeAttribute("y",Double.toString(err[1]));
						xmlwriter.writeAttribute("z",Double.toString(err[2]));
					xmlwriter.endEntity();
				xmlwriter.endEntity();

//				if (verbose > 0) {
//					System.out.printf("%2d       (%12.6f, %12.6f, %12.6f)  (%12.6f, %12.6f, %12.6f)  (%6.3f, %6.3f, %6.3f)\n", i, numgrad[i * 3], numgrad[i * 3 + 1], numgrad[i * 3 + 2], anagrad[i * 3], anagrad[i * 3 + 1], anagrad[i * 3 + 2], err[0], err[1], err[2]);
//				}
			}
			xmlwriter.writeEntity("MaxError");
			xmlwriter.writeText(Double.toString(maxerr));
			xmlwriter.endEntity();
//			if (verbose > 0) {
//				System.out.printf(" Maximum error = %f \n", maxerr);
//				System.out.printf("-------------------------------------------------------------------------------------------------------------------------------\n");
//			}

			xmlwriter.endEntity();
			xmlwriter.close();

			sVerbose=writer.toString();

		} catch (IOException ex) {
			Logger.getLogger(Potential.class.getName()).log(Level.SEVERE, null, ex);
		} 
		return sVerbose;
	}

	/**
	 * @param arg0: unit name (see UnitTable)
	 * @return index (integer value) of unit in UnitTable
	 */

	public static int getIndexOfUnit(String arg0){
		for(int i=0;i<UnitTable.length;i++){
				if(arg0.contentEquals(UnitTable[i])) return i;
			}
		return 0;
	}

	protected static String[] UnitTable={"Hartree","kcal/mol","kJ/mol","eV","au"};

	/**
	 * convert energy value from (from_) unit to (to_)
	 * @param e input energy
	 * @param from original unit
	 * @param to target unit
	 * @return new value of energy in (to_) unit
	 */
	public static double ConvertUnit(double e,String from,  String to){
		if(from.contentEquals(to))	return e;
		int iFrom=getIndexOfUnit(from);
		int iTo=getIndexOfUnit(to);

		final double HART_2_KCAL=627.509391;
		double HART_2_KJ=2625.5;
		double HART_2_EV=27.2113845;
		double KCAL_2_HART=1.593601649e-3;
		double KCAL_2_KJ=4.184;
		double KCAL_2_EV=0.043364107;
		double KJ_2_HART=3.808798324e-4;
		double KJ_2_KCAL=0.2390057;
		double KJ_2_EV=0.010364267;
		double EV_2_HART=0.036749324;
		double EV_2_KCAL=23.06054626;
		double EV_2_KJ=96.48535156;

		double[][] ConvertTable={ {1.0,HART_2_KCAL, HART_2_KJ, HART_2_EV } ,
								  {KCAL_2_HART, 1.0, KCAL_2_KJ, KCAL_2_EV },
					 			  {KJ_2_HART,KJ_2_KCAL,1.0, KJ_2_EV },
								  {EV_2_HART,EV_2_KCAL,EV_2_KJ,1.0} };

	   return e*ConvertTable[iFrom][iTo];
	}
	
	/////////////////////////////////////////////////////////////////////////
    // Energy Evaluation : supposed to be overwriten                       //
    /////////////////////////////////////////////////////////////////////////

	protected double Energy_(final double[] x){ return 1.0;};
	protected double DeltaEnergy_(double[] x,boolean[] isUpdate){ return 1.0;};
	protected void Gradient_(final double[] x, double[] grad){
		NumericalGradient_(x,grad,1e-4);
	};

	protected void Hessian_(final double[] x, double[][] hessian){
		NumericalHessian_(x,hessian,1e-5);
	};

	protected void NumericalGradient_(final double[] xt, double[] grad,double delta){
		double e1,e2;
		double[] x=xt.clone();

	//for(int i=0;i<_ncoords;i++) x[i]=xt[i];
	//System.out.printf(" No coords = "<<_ncoords<<endl;
		double rms=0;
		for(int i=0;i<cluster.getNcoords();i++){
			x[i]+=delta;	e1=Energy_(x);	x[i]-=delta;
			x[i]-=delta;	e2=Energy_(x);	x[i]+=delta;
			//System.out.printf(i<<" = "<<e1<<" "<<e2<<" = "<<(e1-e2)/(2.0*delta)<<endl;
			
				grad[i]=(e1-e2)/(2.0*delta);			

		}
	}

	protected void NumericalHessian_(final double[] _p, double[][] hessian,double delta){
		assert delta>0;
		assert _p.length>0;

		double[] xDelta = _p.clone();
//		double[] xDelta=new double[_p.length];
//		for(int i=0;i<_p.length;i++) xDelta[i]=_p[i];
		
		if(HasAnalyticalGradients()){
			double[][] d1=new double[_p.length][_p.length];
			double[][] d2=new double[_p.length][_p.length];
			for(int i=0;i<_p.length;i++){ //using central finite difference
				xDelta[i]+=delta;
				Gradient_(xDelta,d1[i]);
//					System.err.printf("x : ");			MTools.PrintArray(xDelta);
//					System.err.printf("D1[%d] : ",i);	MTools.PrintArray(d1[i]);
				xDelta[i]-=delta;
				
				xDelta[i]-=delta;
				Gradient_(xDelta,d2[i]);
//					System.err.printf("x : ");			MTools.PrintArray(xDelta);
//					System.err.printf("D2[%d] : ",i);	MTools.PrintArray(d2[i]);
				xDelta[i]+=delta;				
				
			}

			for(int i=0;i<_p.length;i++)
				for(int j=i;j<_p.length;j++){
					//_gg[l][i]=(dfold1[l]-dfold2[l])/(2.0*delta);
					hessian[i][j]=hessian[j][i]=(d1[i][j]-d2[i][j] + d1[j][i]-d2[j][i] )/(4.0*delta);
				}
//				System.err.printf("X \n");
//				MTools.PrintArray(xDelta);
//				System.err.printf("Hessian \n");
//				MTools.PrintArray(hessian);

//			double[] gradientAtXplusDelta = new double[_p.length];
//			double[] gradientAtXminusDelta = new double[_p.length];
//			for (int i = 0; i < _p.length; i++) {
//				xDelta[i]+=delta;	Gradient_(xDelta,gradientAtXplusDelta);		xDelta[i]-=delta;
//				xDelta[i]-=delta;	Gradient_(xDelta,gradientAtXminusDelta);	xDelta[i]+=delta;
//				for (int j = 0; j < _p.length; j++) {
//					hessian[i][j] = (gradientAtXplusDelta[j] - gradientAtXminusDelta[j]) / (2 * delta);
//				}
//			}
		}else{ //using forward finite difference
			double fx=Energy_(xDelta);
			double[] f1=new double[_p.length];
			double[][] d1=new double[_p.length][_p.length];
			for(int i=0;i<_p.length;i++){
				xDelta[i]+=delta;
				f1[i]=Energy_(_p);

				for(int j=i;j<_p.length;j++){
					xDelta[j]+=delta;
					d1[i][j]= Energy_(_p);
					xDelta[j]-=delta;
				}
				_p[i]-=delta;
			}
			for(int i=0;i<_p.length;i++)
				for(int j=i;j<_p.length;j++){
						hessian[i][j]= hessian[j][i]=(d1[i][j]-f1[i]-f1[j] +fx)/(delta*delta);
				}
		}
	}

    /////////////////////////////////////////////////////////////////////////
    // Logging                                                             //
    /////////////////////////////////////////////////////////////////////////
	/**
	 * Return info of potential in xml format
	 * @param verbose =1: equation name and unit, 2= setting info, 3= status
	 */
	public String Info(int verbose){//!< print setting parameters
		String info = "";
		try {
			//!< print setting parameters
			Writer writer = new java.io.StringWriter();
			XmlWriter xmlwriter = new XmlWriter(writer);
			xmlwriter.writeEntity("Potential");
			if ((verbose & 1) !=0) {
				xmlwriter.writeAttribute("Equation", getEquation()).writeAttribute("Unit",getUnit());
			}

			if( (verbose & 2) != 0){
				xmlwriter.writeAttribute("MaxEvals", Integer.toString(nMaxEvals));
				xmlwriter.writeAttribute("EnergyTolerance", Double.toString(EnergyTol));
				xmlwriter.writeAttribute("GradTolerance", Double.toString(GradTol));
			}
			xmlwriter.endEntity();
			xmlwriter.close();
			info = writer.toString();

		} catch (IOException ex) {
			Logger.getLogger(Potential.class.getName()).log(Level.SEVERE, null, ex);
		} 
		return info;
	}
        
	//	String basisset; //!< for Gaussian or other DFT calculations

	/////////////////////////////////////////////////////////////////////////
    // Energy Minimization                                                 //
    /////////////////////////////////////////////////////////////////////////
//
//	enum CONVERGENCE_TYPE{
//		ON_ENERGY,
//		ON_RMSGRAD,
//		ON_MAXEVALS,
//		ON_ERROR
//	};

	int linesearch; //!< LineSearch type
    /**
	 * Set the LineSearchType. The default type is LineSearchType::Simple.
     *  @param type The LineSearchType to be used in SteepestDescent and ConjugateGradients.
     */
    void setLineSearchType(int type){
      linesearch = type;
    }
	
    /**
	 * Get the LineSearchType.
     * @return The current LineSearchType.
     */
    int getLineSearchType(){
      return linesearch;
    }

	

	int nMaxEvals=10000; //!< number of maximum evaluations
	public int getMaxEvals(){ return nMaxEvals;}
	public void setMaxEvals(int a){ nMaxEvals=a;};

	double EnergyTol=1e-8;//!< funtional tolerance criteria in optimization
	public double getEnergyTol(){return EnergyTol;}
	public void setEnergyTol(double a){ EnergyTol=a;}

	double GradTol=1e-6;//!< gradient tolerance criteria in optimization
	public double getGradientTol(){return GradTol;}
	public void setGradientTol(double a){ GradTol=a;}

	double MaxStepSize=0.1;//!< maximum step size
	public double getMaxStepSize(){return MaxStepSize;}
	public void setMaxStepSize(double a){ MaxStepSize=a;}

	String optMethod="DFPMIN"; //!< optimization method

	public String getOptMethod(){ return optMethod;}
	public void setOptMethod(String a){ optMethod=a.toUpperCase(); }

	//parameters returned by Optimize
	int nEvals=0; //!< number of evaluations performed
	public int getNEvals(){ return nEvals;}
	public void setNEvals(int a){ nEvals=a; }

	protected double optimizedE = 0;
	
	int iConverged;//!< check whether converged or not
	double nSeconds=0; //!< time(seconds) for optimization
	public double getSeconds(){ return nSeconds;}

	double MaxGrad,RmsGrad;//!< Maximum gradient component and Root Mean Square of gradient, will be updated in optimization

    /**
     * return RMS of gradient of coordinates x. A real gradient is performed
     * @param x coordinates
     * @return
     */
	public double getRMSGrad(double[] x){//!< 
		if(x.length==0) return 1;
		double[] gradient=new double[x.length];
		getGradient(x,gradient);		
		return MTools.CalculateRMS(gradient);
	}

	public double getRMSGrad(){//!< return RMS of gradient
		return RmsGrad;
	}

	public double getMaxGrad(){return MaxGrad;}

	final double EPS=1.0e-12;
	//Convergence criterion on xfinal values.
	final double TOLX=(4*EPS);
	//Ensures sufficient decrease in function value.
	final double ALF=1.0e-4;

	/**
	 * Line search for DFPMIN. During calculation, if nEvals > nMaxEvals, it will terminate. 
	 * @param xinit Starting point
	 * @param finit Starting energy
	 * @param grad gradients
	 * @param direction direction
	 * @param xfinal Resultant point after search
	 * @return Resultant fitness
	 */
	private double LineSearch(double[] xinit, double finit, double[] grad, double[] direction, double[] xfinal){
		double a,alam,alam2=0,alamin,b,disc,f2=0,rhs1,rhs2,tmplam;

		double ffinal=finit;

		double sum=0;
		for (int i=0;i<nCoords;i++) sum += direction[i]*direction[i];
		sum=Math.sqrt(sum);

		if (sum > MaxStepSize)
			for (int i=0;i<nCoords;i++) direction[i] *= MaxStepSize/sum; //Scale if attempted step is too big.

		double slope=0;
		for (int i=0;i<nCoords;i++)	slope += grad[i]*direction[i];
		//if (slope >= 0.0) {};//cout<<"Roundoff problem in lnsrch.\nCoords";


		double test=0.0; //Compute min.

		for (int i=0;i<nCoords;i++) {
			double temp=Math.abs(direction[i])/Math.max(Math.abs(xinit[i]),1.0);
			if (temp > test) test=temp;
		}
		alamin=TOLX/test;
		alam=1.0;  //Always try full Newton step first.

		while ( nEvals < nMaxEvals ) { //Start of iteration loop.
			for (int i=0;i<nCoords;i++) xfinal[i]=xinit[i]+alam*direction[i];

			ffinal=getEnergy(xfinal);


			if (alam < alamin) { //Convergence on xfinal. For zero finding,the calling program should verify the convergence.
				for (int i=0;i<nCoords;i++) xfinal[i]=xinit[i];
				return 1;
			} else
				if (ffinal <= finit+ALF*alam*slope) return ffinal; //Convergence on energy
				else { //Backtrack.
					if (alam == 1.0)
						tmplam = -slope/(2.0*(ffinal-finit-slope)); //First time.
					else {// Subsequent backtracks.
						rhs1 = ffinal- finit - alam*slope;
						rhs2=f2-finit - alam2*slope;
						a=(rhs1/(alam*alam)-rhs2/(alam2*alam2))/(alam-alam2);
						b=(-alam2*rhs1/(alam*alam)+alam*rhs2/(alam2*alam2))/(alam-alam2);
						if (a == 0.0) tmplam = -slope/(2.0*b);
						else {
							disc=b*b-3.0*a*slope;
							if (disc < 0.0) tmplam=0.5*alam;
							else if (b <= 0.0)
									tmplam=(-b+Math.sqrt(disc))/(3.0*a);
								 else
									 tmplam=-slope/(b+Math.sqrt(disc));
						}
						if (tmplam > 0.5*alam)
							tmplam=0.5*alam;
					}
				}
			alam2=alam;
			f2 = ffinal;
			alam=Math.max(tmplam,0.1*alam);
		}

		return ffinal;
	}


	/**
	 * DFPmin optimizer. Will update optimizedE, MaxRMS, RMSGrad and nEvals
	 * @param x Starting point. Will be optimized
	 * @return Resultant fitness
	 */
	private int DFPmin(double p[]){
		int i,its,j,iter;
		double fac,fad,fae,fp,sumdg,sumxi,test=0;
		optimizedE=0;

		double[] dg=new double [nCoords];
		double[] g=new double [nCoords];
		double[] hdg=new double [nCoords];
		double[][] hessin=new double [nCoords][nCoords];
		double[] pnew=new double [nCoords];
		double[] direction=new double [nCoords];

		fp=getEnergy(p);  //Calculate starting function value and gradient,
		optimizedE=fp;
		//System.err.printf(" ffinal = %f \n",ffinal);
		getGradient(p,g);	

		double sum=0;
		for (i=0;i<nCoords;i++) { //and initialize the inverse Hessian to the unit matrix.
			for (j=0;j<nCoords;j++) hessin[i][j]=0.0;
			hessin[i][i]=1.0;
			direction[i] = -g[i]; //Initial line direction.
			sum += p[i]*p[i];
		}
//		try {
//			FileWriter fileOut = new FileWriter(new File("opt.xyz"));

			int step=0;
		while ( nEvals < nMaxEvals ){ //Main loop over the iterations.
			optimizedE=LineSearch(p,fp,g,direction,pnew);
			
			
			//System.err.printf(" fret = %f \n",optimizedE);

			//The new function evaluation occurs in lnsrch; save the function value in fp for the
			//next line search. It is usually safe to ignore the value of check.
			fp = optimizedE;
			for (i=0;i<nCoords;i++) {
				direction[i]=pnew[i]-p[i]; //Update the line direction,
				p[i]=pnew[i];  // and the current point.
			}


			test=0.0; //Test for convergence on xfinal.
			for (i=0;i<nCoords;i++) {
				double temp=Math.abs(direction[i]);
				if (temp > test) test=temp;
			}


			for (i=0;i<nCoords;i++) dg[i]=g[i]; //Save the old gradient,
			getGradient(p,g);  //and get the new gradient.

			RmsGrad=0; //Test for convergence on zero gradient
			MaxGrad=0; //test for maximum force.
			for (i=0;i<nCoords;i++) {
				RmsGrad+=g[i]*g[i];
				MaxGrad=Math.max(Math.abs(g[i]),MaxGrad);
			}
			RmsGrad=Math.sqrt(RmsGrad/nCoords);

//			cluster.setCoords(pnew);
//			cluster.setTag(step);
//			cluster.setEnergy(fp);
//			cluster.setRmsGrad(RmsGrad);
//			cluster.Write(fileOut,"xyz");

			//if (test < EnergyTol){ System.out.print(" Die because of energy tol\n");  return 1;}
			if ((RmsGrad < GradTol) && (MaxGrad < GradTol*4.0)){  return 2;}

			for (i=0;i<nCoords;i++) dg[i]=g[i]-dg[i]; //Compute difference of gradients,
			for (i=0;i<nCoords;i++) { //and difference times current matrix.
				hdg[i]=0.0;
				for (j=0;j<nCoords;j++) hdg[i] += hessin[i][j]*dg[j];
			}
			fac=fae=sumdg=sumxi=0.0; //Calculate dot products for the denominators.
			for (i=0;i<nCoords;i++) {
				fac += dg[i]*direction[i];
				fae += dg[i]*hdg[i];
				sumdg += Math.pow(dg[i],2);
				sumxi += Math.pow(direction[i],2);
			}
			if (fac > Math.sqrt(EPS*sumdg*sumxi)) { //Skip update if fac not sufficiently positive.
				fac=1.0/fac;
				fad=1.0/fae;
		//The vector that makes BFGS different from DFP:
				for (i=0;i<nCoords;i++) dg[i]=fac*direction[i]-fad*hdg[i];
				for (i=0;i<nCoords;i++) { //The BFGS updating formula:
					for (j=i;j<nCoords;j++) {
						hessin[i][j] += fac*direction[i]*direction[j] -fad*hdg[i]*hdg[j]+fae*dg[i]*dg[j];
						hessin[j][i]=hessin[i][j];
					}
				}
			}
			for (i=0;i<nCoords;i++) { //Now calculate the next direction to go,
				direction[i]=0.0;
				for (j=0;j<nCoords;j++) direction[i] -= hessin[i][j]*g[j];
			}
			step++;
		}// and go back for another iteration.

//			} catch (IOException ex) {
//			Logger.getLogger(Potential.class.getName()).log(Level.SEVERE, null, ex);
//		}

		return 0;//if overrun
	}

	//============================== Functions which are integrated with apache common opt
	class Convergence implements RealConvergenceChecker{

		@Override
		public boolean converged(int arg0, RealPointValuePair arg1, RealPointValuePair arg2) {
			if(arg2.getValue() - arg1.getValue() < EnergyTol ) return true;

			double rms=getRMSGrad(arg2.getPoint());
			if(rms < GradTol) return true;
			return false;
		}

	}

	class DiffFunc implements DifferentiableMultivariateRealFunction  {
		class Gradient implements MultivariateVectorialFunction{
			@Override
			public double[] value(double[] arg0) throws FunctionEvaluationException, IllegalArgumentException {
				return getGradient(arg0);
			}
		}

		Gradient gradient=new Gradient();

		@Override
		public MultivariateRealFunction partialDerivative(int arg0) {
			throw new UnsupportedOperationException("Not supported yet.");
		}

		@Override
		public MultivariateVectorialFunction gradient() {
			return  gradient;
		}

		@Override
		public double value(double[] arg0) throws FunctionEvaluationException, IllegalArgumentException {
			return getEnergy(arg0);
		}


	}
}
