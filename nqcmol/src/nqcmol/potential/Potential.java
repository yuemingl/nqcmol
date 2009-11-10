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
import nqcmol.xml.XmlWriter;

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
	protected Cluster cluster = new Cluster();

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
	public void setCluster(Cluster cluster) {
		this.cluster = (Cluster) cluster.clone();
		nCoords=cluster.getNcoords();
	}


	int nCoords=0;
	
	public int getNCoords(){
		return cluster.getNcoords();
	}
    // minimization variables
	

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
	 * Does this force field have analytical gradients.
     * @return True if all analytical gradients are implemented.
     */
    public boolean HasAnalyticalGradients() { return false; }

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
	 * setup mol  and then calculate energy
	 * @param isMolUpdate if true, new energy will be updated to mol
	 * @return calculated energy
	 */
	public double getEnergy(Cluster cluster_,boolean isMolUpdate) {
		setCluster(cluster_);
		double ener=Energy_(cluster.getCoords());
		cluster.setEnergy(ener);
		if(isMolUpdate) cluster_.setEnergy(ener);
		return ener;
	}

	/**
	 * return current energy of internal cluster
	 * @param isEnergyUpdate: if true, it will call Energy function to update
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
	 * return current energy. Do not update the internal cluster.  Require setup mol first.
	 *
	 * @return energy
	 */
	public double getEnergy(final double[] coords_){
		return Energy_(coords_) ;
	}

	
	public double[] getGradient(final double[] coords_){//!< require setup mol first, if update==true, recalculate the gradients
		cluster.setCoords(coords_);
		Gradient_(cluster.getCoords(),cluster.getGradient());
		return cluster.getGradient();
	}

	public void getGradient(final double[] coords_,double[] gradient_){//!< require setup mol first, if update==true, recalculate the gradients
		Gradient_(coords_,gradient_);
	}

	public double[] getGradient(boolean isGradUpdate){//!< require setup mol first, if update==true, recalculate the gradients
		if(isGradUpdate){
			Gradient_(cluster.getCoords(),cluster.getGradient());
		}
		return cluster.getGradient();
	}


	public double getDeltaEnergy(final double[] coords_,final boolean[] isUpdate){//!< require setup mol first
		return DeltaEnergy_(coords_,isUpdate) ;
	}


	public boolean Optimize(){//!< require setup mol first
		//setting the variable
		nEvals=0;

//		if(debug>=2){
//			cout<<"Dimension = "<<Dim_();
//			cout<<" Energy="<<evaluate_(x)<<" Method="<<opt_med<<endl;
//		}

//	if(opt_med==DFPMIN){
//		double *xb=new double[dim];
//		for(int i=0;i<dim;i++)	xb[i]=x[i];
//		iConverged=dfpmin(xb,dim,nMaxEval, &final_eval);
//		for(int i=0;i<dim;i++)	final_x[i]=xb[i];
//		delete[] xb;
//		return neval;
//	}
	//cout<<(h.checkConvg()?"Converged":"Not Yet")<<endl;
	//cout<<"Time = "<<h.getSeconds()<<endl;
	//cout<<"nEval = "<<h.getNEvals()<<endl;

	//start optmization

//
//		switch(opt_med){
//			default:
//				//iConverged=optimize_QNewton(nlp);
//				iConverged=optimize_nlf1<OptQNewton>(nlp);
//			break;
//			case CONJUGATEGRAD :
//				//iConverged=optimize_CG(nlp);
//				iConverged=optimize_nlf1<OptCG>(nlp);
//			break;
//			case LBFGS:
//				iConverged=optimize_nlf1<OptLBFGS>(nlp);
//			break;
//			case FDNEWTON:
//				iConverged=optimize_nlf1<OptFDNewton>(nlp);
//			break;
//
//
//		}
		Double fret=new Double(0);
		iConverged=DFPmin(cluster.getCoords());

			System.err.printf(" Final fret = %f \n",fret);
		cluster.setEnergy(optimizedE);

		return (iConverged!=0);
	}

	public boolean Optimize(Cluster cluster_,boolean isMolUpdate){//!< setup mol  and then optimize mol, if update==true, new energy, structures will be updated to mol
		this.setCluster(cluster_);
		Optimize();
		if(isMolUpdate){
			cluster_.setCoords(this.cluster.getCoords());
			cluster_.setEnergy(this.cluster.getEnergy());
		}
		return true;
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
						err[k] = 100.0;
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


	protected static String[] UnitTable={"Hartree","kcal/mol","kJ/mol","eV","au"};

	public static double ConvertUnit(double e,String from,  String to){//!< convert e from (from_) unit to (to_)
		int iFrom=0;
		int iTo=0;
		for(int i=0;i<UnitTable.length;i++){
			if(from.contentEquals(UnitTable[i])) iFrom=i;
			if(to.contentEquals(UnitTable[i])) iTo=i;
		}

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
			if(!HasAnalyticalGradients()) NumericalGradient_(x,grad,1e-4);
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

	protected void NumericalHessian_(final double[] _p, double[][] _hessian,double delta){
		
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

	int opt_med; //!< optimization method

	int nMaxEvals=100; //!< number of maximum evaluations
	int getMaxEvals(){ return nMaxEvals;}
	void setMaxEvals(int a){ nMaxEvals=a;};

	double EnergyTol=1e-8;//!< funtional tolerance criteria in optimization
	double getEnergyTol(){return EnergyTol;}
	void setEnergyTol(double a){ EnergyTol=a;}

	double GradTol=1e-4;//!< gradient tolerance criteria in optimization
	double getGradientTol(){return GradTol;}
	void setGradientTol(double a){ GradTol=a;}

	double MaxStepSize=0.01;//!< maximum step size
	double getMaxStepSize(){return MaxStepSize;}
	void setMaxStepSize(double a){ MaxStepSize=a;}

	int getOptMethod(){ return opt_med;}
	void setOptMethod(int a){ opt_med=a; }



	//parameters returned by Optimize
	int nEvals=0; //!< number of evaluations performed
	int getNEvals(){ return nEvals;}
	void setNEvals(int a){ nEvals=a; }

	protected double optimizedE = 0;

	/**
	 * Get the value of finalE
	 *
	 * @return the value of finalE
	 */
	public double getOptimizedEnergy() {
		return optimizedE;
	}

	/**
	 * Get the value of OptimizedCoords
	 *
	 * @return the value of OptimizedCoords
	 */
	public double[] getOptimizedCoords() {
		return cluster.getCoords();
	}



	
	int iConverged;//!< check whether converged or not
	double nSeconds; //!< time(seconds) for optimization
	double getSeconds(){ return nSeconds;}

	double MaxGrad,RmsGrad;//!< Maximum gradient component and Root Mean Square of gradient, will be updated in optimization
	double getRMSGrad(){//!< return RMS of graident
		return RmsGrad;
	}

	double getMaxGrad(){return MaxGrad;}

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

			if (test < EnergyTol)  return 1;
			if ((RmsGrad < GradTol) && (MaxGrad < GradTol*4.0)) return 2;

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
		}// and go back for another iteration.

		return 0;//if overrun
	}

}
