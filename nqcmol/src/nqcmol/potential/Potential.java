/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.potential;

import nqcmol.Cluster;

/**
 *
 * @author nqc
 */
public class Potential {
	enum OPT_METHODS{ DFPMIN, QNEWTON, CONJUGATEGRAD, LBFGS, FDNEWTON };

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
	 * Set the value of cluster
	 *
	 * @param cluster new value of cluster
	 */
	public void setCluster(Cluster cluster) {
		this.cluster = new Cluster();
	}


    // minimization variables
	int nMaxEval; //!< number of maximum evaluations
	double FcnTol;//!< funtional tolerance criteria in optimization
	double GradTol;//!< gradient tolerance criteria in optimization
	double MaxStep;//!< maximum step size

	int _linesearch; //!< LineSearch type
	int opt_med; //!< optimization method

	double MaxGrad,RmsGrad;//!< Maximum gradient component and Root Mean Square of gradient, will be updated in optimization
	int iConverged;//!< check whether converged or not
	int neval; //!< number of evaluations performed
	double nsecs; //!< time(seconds) for optimization


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
	 * Set the value of unit
	 *
	 * @param unit new value of unit
	 */
	public void setUnit(String unit) {
		this.unit = unit;
	}

    /* Does this force field have analytical gradients defined for all
     * calculation components (bonds, angles, non-bonded, etc.)
     * If this is true, code should default to using OBFF_ANALYTICAL_GRADIENT
     * for SteepestDescent() or ConjugateGradients().
     * \return True if all analytical gradients are implemented.
     */
    public boolean HasAnalyticalGradients() { return false; }

	protected double[] param;

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

	int getMaxEval(){ return nMaxEval;}
	void setMaxEval(int a){ nMaxEval=a;};

	double getFcnTol(){return FcnTol;}
	void setFcnTol(double a){ FcnTol=a;}

	double getGradTol(){return GradTol;}
	void setGradTol(double a){ GradTol=a;}

	double getMaxStep(){return MaxStep;}
	void setMaxStep(double a){ MaxStep=a;}

	int getOptMethod(){ return opt_med;}
	void setOptMethod(int a){ opt_med=a; }

	int getNEvals(){ return neval;}
	double getSeconds(){ return nsecs;}

	double getRMSGrad(){//!< return RMS of graident
		return RmsGrad;
	}

	double getMaxGrad(){return MaxGrad;}


    /////////////////////////////////////////////////////////////////////////
    // Energy Evaluation                                                   //
    /////////////////////////////////////////////////////////////////////////

	/**
	 * setup mol  and then calculate energy
	 * @param isMolUpdate if true, new energy will be updated to mol
	 * @return calculated energy
	 */
	double getEnergy(Cluster cluster_,boolean isMolUpdate) {
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

	double getEnergy(boolean isEnergyUpdate){
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
	double getEnergy(final double[] coords_){
		return Energy_(coords_) ;
	}

	
	double[] getGradient(final double[] coords_){//!< require setup mol first, if update==true, recalculate the gradients
		cluster.setCoords(coords_);
		Gradient_(cluster.getCoords(),cluster.getGradient());
		return cluster.getGradient();
	}

	void getGradient(final double[] coords_,double[] gradient_){//!< require setup mol first, if update==true, recalculate the gradients
		Gradient_(coords_,gradient_);
	}

	double[] getGradient(boolean isGradUpdate){//!< require setup mol first, if update==true, recalculate the gradients
		if(isGradUpdate){
			Gradient_(cluster.getCoords(),cluster.getGradient());
		}
		return cluster.getGradient();
	}


	double getDeltaEnergy(final double[] coords_,final boolean[] isUpdate){//!< require setup mol first
		return DeltaEnergy_(coords_,isUpdate) ;
	}


	boolean Optimize(){//!< require setup mol first

	}

	boolean Optimize(Cluster cluster_,boolean isMolUpdate){//!< setup mol  and then optimize mol, if update==true, new energy, structures will be updated to mol
		this.setCluster(cluster_);
		Optimize();
		if(isMolUpdate){
			cluster_.setCoords(this.cluster.getCoords());
			cluster_.setEnergy(this.cluster.getEnergy());
		}
		return true;
	}

	static public String[] UnitTable={"Hartree","kcal/mol","kJ/mol","eV","au"};

	static double ConvertUnit(double e,String from,  String to){//!< convert e from (from_) unit to (to_)
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
	
	protected double Energy_(final double[] x){ return 1.0;};
	protected double DeltaEnergy_(double[] x,boolean[] isUpdate){ return 1.0;};
	protected void Gradient_(final double[] x, double[] grad){
			if(!HasAnalyticalGradients()) NumericalGradient_(x,grad,1e-4);
		};

	protected void NumericalGradient_(final double[] xt, double[] grad,double delta){
		double e1,e2;
		double[] x=xt.clone();

	//for(int i=0;i<_ncoords;i++) x[i]=xt[i];
	//cout<<" No coords = "<<_ncoords<<endl;
		double rms=0;
		for(int i=0;i<cluster.getNcoords();i++){
			x[i]+=delta;	e1=Energy_(x);	x[i]-=delta;
			x[i]-=delta;	e2=Energy_(x);	x[i]+=delta;
			//cout<<i<<" = "<<e1<<" "<<e2<<" = "<<(e1-e2)/(2.0*delta)<<endl;
			
				grad[i]=(e1-e2)/(2.0*delta);			
	}
	}

	protected void NumericalHessian_(final double[] _p, double[][] _hessian,double delta){
		
	}

    /////////////////////////////////////////////////////////////////////////
    // Logging                                                             //
    /////////////////////////////////////////////////////////////////////////
	//! \name Methods for logging
	//@{
	public string Info(int verbose=0);//!< print setting parameters

	
    //@}


    /////////////////////////////////////////////////////////////////////////
    // Energy Minimization                                                 //
    /////////////////////////////////////////////////////////////////////////

    //! \name Methods for energy minimization
    //@{
    /*! Set the LineSearchType. The default type is LineSearchType::Simple.
     *  \param type The LineSearchType to be used in SteepestDescent and ConjugateGradients.
     */
    void setLineSearchType(int type){
      _linesearch = type;
    }
    /*! Get the LineSearchType.
     *  \return The current LineSearchType.
     */
    int getLineSearchType(){
      return _linesearch;
    }

    /////////////////////////////////////////////////////////////////////////
    // Validation                                                          //
    /////////////////////////////////////////////////////////////////////////

    //! \name Methods for forcefield validation
    /*!
      Validate the analytical gradients by comparing them to numerical ones. This function has to
      be implemented force field specific. (debugging)
    */
    double ValidateGradients(int verbose=1);
    //@}

		String basisset; //!< for Gaussian or other DFT calculations

	private	void lnsrch(int n, double[] xold, double fold, double g[], double p[], double x[],double *fret, double stpmax, int *check);
	private	int dfpmin(double p[], int n, int iterMax, double *fret);
}
