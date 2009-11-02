/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.openscience.cdk.modeling.forcefield;

import javax.vecmath.GMatrix;
import javax.vecmath.GVector;
import org.openscience.cdk.interfaces.IAtomContainer;


/**
 *
 * @author nqc
 */
public class OSS2Function implements IPotentialFunction {
	String energyFunctionShape = " Lennard-Jones energy ";
	double energy = 0;
	GVector energyGradient = null;
	GVector order2ndErrorApproximateGradient = null;
	GMatrix energyHessian = null;
	double[] forHessian = null;
	GMatrix order2ndErrorApproximateHessian = null;
	double[] forOrder2ndErrorApproximateHessian = null;
	int functionEvaluationNumber = 0;

	//private LoggingTool logger;


	/**
	 *  Constructor 
	 *
	 */
	public OSS2Function()  {
		//logger.debug("LJEnergyFunction Constructor");
		//logger = new LoggingTool(this);
		double[] t={0.9614,	1.81733805,	1.73089134,	0.29575999,	332.2852922, 2.84683483,	1.0956864,	0.00000204,	0.00575558,	2.6425374,	1.1274385, 3.3255655,	0.07847734,	40.4858734,	1.3290955,	-41.7260608, 1.3509491, 0.0629396,	6.77909903,	1.8117836,	-0.0420073,	0.16488265,	-0.02509795, -0.37814525, -0.31667595, -0.0114672, 0.06061075, 0.528281, 1.1617627, 0.0765232, -0.21208835, -0.1025385, -0.0762216, -0.228694, 0, -0.029092, 6.25, 0.00377233,	6.25, 1.444 };

 		for(int i=0;i<40;i++)	 alpha[i]=t[i];

		ParseParameters();
	}

	public OSS2Function(IAtomContainer mol)  {
		//logger.debug("LJEnergyFunction Constructor");
		//logger = new LoggingTool(this);
		nAtom=mol.getAtomCount();
		nO=(nAtom%3 == 2)?(nAtom/3+1):(nAtom/3);
		x=new double[nAtom][3];
		q=new double[nAtom];
		r=new double[nAtom][nAtom];
		r2=new double[nAtom][nAtom];
		isUpdate=new boolean[nAtom];
		Vpairwise=new double[nAtom][nAtom];

		Vq=new double[nAtom][nAtom];
		VHOH=new double[nO][nAtom][nAtom];

		Scd=new double[nAtom][nAtom];
		Sdd=new double[nO][nO];
		D=new double[nO*3][nO*3];
		T=new double[nO*3][nO*3];
		E=new double[nO*3];
		mu=new double[nO*3];

		for(int i=0; i<nAtom; i++ ){
			q[i] = ((i<nO)? -2 : 1);
			isUpdate[i]=true;		
		}

		double[] t={0.9614,	1.81733805,	1.73089134,	0.29575999,	332.2852922, 2.84683483,	1.0956864,	0.00000204,	0.00575558,	2.6425374,	1.1274385, 3.3255655,	0.07847734,	40.4858734,	1.3290955,	-41.7260608, 1.3509491, 0.0629396,	6.77909903,	1.8117836,	-0.0420073,	0.16488265,	-0.02509795, -0.37814525, -0.31667595, -0.0114672, 0.06061075, 0.528281, 1.1617627, 0.0765232, -0.21208835, -0.1025385, -0.0762216, -0.228694, 0, -0.029092, 6.25, 0.00377233,	6.25, 1.444 };

 		for(int i=0;i<40;i++)	 alpha[i]=t[i];

		ParseParameters();
	}


	@Override
	public String toString(){
		return "OSS2";
	}


	/**
	 *  Evaluate the LJEnergyFunction for a given molecule
	 *
	 *@param  molecule  Current molecule.
	 *@return        LennardJones energy function value.
	 */
	public double energyFunctionOfAMolecule(IAtomContainer molecule) {

		GVector coords3d = ForceFieldTools.getCoordinates3xNVector(molecule);

		energy = energyFunction(coords3d);

		return energy;
	}



	/**
	 *  Evaluate the LennardJones energy function for a given 3xN point
	 *
	 *@param  coords3d  Current molecule 3xN coordinates.
	 *@return        LennardJones energy function value.
	 */
	public double energyFunction(GVector _p) {
		energy=0;

		//System.out.print(" Size of vec = "+ Integer.toString(coords3d.getSize()));

		for(int i=0;i<nAtom;i++) if(isUpdate[i]) {
			x[i][0]= _p.getElement(i*3+0);
			x[i][1]= _p.getElement(i*3+1);
			x[i][2]= _p.getElement(i*3+2);
		// cout<<i<<" d "<<x[i][0]<<" - "<<x[i][1]<<" -  "<<x[i][2]<<endl;
		}

		CalcDistance();
		CalcPairWiseInteraction();
		CalcThreeBodyInteraction();
		CalcChargeChargeInteraction();

		CalcCutoffFunc();
		CalcDipoleMoment();


		energy=Vpairwise_total+VHOH_total;

		System.out.printf( "VOO: %f  VOH: %f VHOH: %f VHH:%f \n",VOO,VOH,VHOH_total,VHH);
		System.out.printf( "Vq: %f \n",Vq_total);
//		System.out.printf( "Self-dipole: " <<setprecision(12)<< Vsd <<endl;
//		System.out.printf( "Dipole-dipole interaction: " << setprecision(12)<<Vdd << endl;
//		System.out.printf( "Charge-dipole interaction: " << setprecision(12)<<Vcd << endl;
//		System.out.printf( "Vel: " << setprecision(12)<<Vel << endl;
		System.out.printf( "Energy: %f \n",energy);


		return energy;
	}


	/**
	 *  Evaluate the gradient for the LennardJones energy function in a given 3xN point
	 *
	 *@param  coords3d  Current molecule coordinates.
	 */
	public void setEnergyGradient(GVector _p) {
		//setOrder2ndErrorApproximateEnergyGradient(coords3d);

		//logger.debug("coords3d : " + coords3d);
		energyGradient = new GVector(_p.getSize());

		
		//logger.debug("bs.getGradientLennardJonesSumEB() = " + bs.getGradientLennardJonesSumEB());
		//logger.debug("ab.getGradientLennardJonesSumEA() = " + ab.getGradientLennardJonesSumEA());
		int nAtom=_p.getSize()/3;
		for (int i=0; i < nAtom ; i++) {
//			energyGradient.setElement(i,
			double[] a={0,0,0};
			for(int j=0;j<nAtom;j++)	if(i!=j){
				double rij2=  Math.pow( (_p.getElement(i*3+0)-_p.getElement(j*3+0) ),2)
							+ Math.pow( (_p.getElement(i*3+1)-_p.getElement(j*3+1) ),2)
							+ Math.pow( (_p.getElement(i*3+2)-_p.getElement(j*3+2) ),2);
				double rij6 = Math.pow(rij2,3);

				double rij14= rij6*rij6*rij2;
				double f=24.0*(rij6-2.0)/(rij14);
				for(int k=0;k<3;k++)	a[k]+=f*(_p.getElement(i*3+k)-_p.getElement(j*3+k));
			}
			for(int k=0;k<3;k++)	energyGradient.setElement(i*3+k,a[k]);
		}		
	}


	/**
	 *  Get the gradient for the potential energy function in a given point
	 *
	 *@return        Gradient value
	 */
	public GVector getEnergyGradient() {
		return energyGradient;
		//return order2ndErrorApproximateGradient;
	}


	/**
	 *  Evaluate a 2nd order approximation of the gradient
	 *
	 *
	 *@param  coord3d  Current molecule coordinates.
	 */
	public void set2ndOrderErrorApproximateGradient(GVector point) {
		order2ndErrorApproximateGradient = new GVector(point.getSize());
		double sigma = Math.pow(0.000000000000001,0.33);
		GVector xplusSigma = new GVector(point.getSize());
		GVector xminusSigma = new GVector(point.getSize());

		for (int m = 0; m < order2ndErrorApproximateGradient.getSize(); m++) {
			xplusSigma.set(point);
			xplusSigma.setElement(m,point.getElement(m) + sigma);
			xminusSigma.set(point);
			xminusSigma.setElement(m,point.getElement(m) - sigma);
			order2ndErrorApproximateGradient.setElement(m,(energyFunction(xplusSigma) - energyFunction(xminusSigma)) / (2 * sigma));
		}

		//logger.debug("order2ndErrorApproximateGradient : " + order2ndErrorApproximateGradient);
	}


	/**
	 *  Get the 2nd order error approximate gradient of the angle bending term.
	 *
	 *
	 *@return           Angle bending 2nd order error approximate gradient value.
	 */
	public GVector get2ndOrderErrorApproximateGradient() {
		return order2ndErrorApproximateGradient;
	}


	/**
	 *  Evaluate the hessian for the potential energy function in a given point
	 *
	 *@param  point  Current molecule coordinates.
	 */
	public void setEnergyHessian(GVector point) {
		/*double[] forHessian = {2,0,0,0,4,0,0,0,1};
		energyHessian = new GMatrix(3, 3, forHessian);
		*/
		set2ndOrderErrorApproximateHessian(point);
	}


	/**
	 *  Get the hessian for the potential energy function in a given point
	 *
	 *@return        Hessian value
	 */
	public GMatrix getEnergyHessian() {
		//return energyHessian;
		return order2ndErrorApproximateHessian;
	}


	/**
	 *  Get the hessian of the potential energy function in a given point.
	 *
	 *@return        Hessian energy value in the wished point.
	 */
	 public double[] getForEnergyHessian() {
		 //return forHessian;
		 return forOrder2ndErrorApproximateHessian;
	 }


	/**
	 *  Evaluate a 2nd order approximation of the Hessian
	 *  given the atoms coordinates.
	 *
	 *@param  coord3d  Current molecule coordinates.
	 */
	public void set2ndOrderErrorApproximateHessian(GVector coord3d) {
		forOrder2ndErrorApproximateHessian = new double[coord3d.getSize() * coord3d.getSize()];
		double sigma = Math.pow(0.000000000000001,0.33);
		GVector xplusSigma = new GVector(coord3d.getSize());
		GVector xminusSigma = new GVector(coord3d.getSize());
		GVector gradientAtXplusSigma = new GVector(coord3d.getSize());
		GVector gradientAtXminusSigma = new GVector(coord3d.getSize());

		int forHessianIndex;
		for (int i = 0; i < coord3d.getSize(); i++) {
			xplusSigma.set(coord3d);
			xplusSigma.setElement(i,coord3d.getElement(i) + sigma);
			setEnergyGradient(xplusSigma);
			gradientAtXplusSigma.set(this.getEnergyGradient());
			xminusSigma.set(coord3d);
			xminusSigma.setElement(i,coord3d.getElement(i) - sigma);
			setEnergyGradient(xminusSigma);
			gradientAtXminusSigma.set(this.getEnergyGradient());
			for (int j = 0; j < coord3d.getSize(); j++) {
				forHessianIndex = i*coord3d.getSize()+j;
				forOrder2ndErrorApproximateHessian[forHessianIndex] = (gradientAtXplusSigma.getElement(j) - gradientAtXminusSigma.getElement(j)) / (2 * sigma);
				//(energyFunction(xplusSigma) - 2 * fx + energyFunction(xminusSigma)) / Math.pow(sigma,2);
			}
		}
		forOrder2ndErrorApproximateHessian[8] = 1;
		order2ndErrorApproximateHessian = new GMatrix(coord3d.getSize(), coord3d.getSize());
		order2ndErrorApproximateHessian.set(forOrder2ndErrorApproximateHessian);
		//logger.debug("order2ndErrorApproximateHessian : " + order2ndErrorApproximateHessian);
	}


	/**
	 *  Get the 2nd order error approximate Hessian
	 *
	 *
	 *@return           2nd order error approximate Hessian value.
	 */
	public GMatrix get2ndOrderErrorApproximateHessian() {
		return order2ndErrorApproximateHessian;
	}


	//========================== implement OSS2
	int nAtom,nO;
	double[] q;
	double[][] r,x,r2;
	boolean[] isUpdate;

	double[][] Vpairwise;
	double Vpairwise_total,VHOH_total,Vq_total;

	double VOO,VOH,VHH;

	double[][][] VHOH;
	double[][] Vq;
	double[] E,mu;

	double[][] Scd,Sdd,D,T;


	double p_r0;
	double p_theta0;
	double[] p_a=new double[3];
	double[] p_b=new double[3];
	double[] p_c=new double[3];
	double[] p_h=new double[6];
	double[] p_o=new double[8];
	double[] p_k=new double[17];
	double[] p_m=new double[4];
	double p_alpha;

	double[] alpha=new double[40];

	void ParseParameters(){
		 p_r0= alpha[0];
		p_theta0= alpha[1];
 p_a[1]= alpha[2];
 p_a[2]= alpha[3];
 p_b[1]= alpha[4];
 p_b[2]= alpha[5];
 p_c[1]= alpha[6];
 p_c[2]= alpha[7];

 p_h[1]= alpha[8];
 p_h[2]= alpha[9];
 p_h[3]= alpha[10];
 p_h[4]= alpha[11];
 p_h[5]= alpha[12];

 p_o[1]= alpha[13];
 p_o[2]= alpha[14];
 p_o[3]= alpha[15];
 p_o[4]= alpha[16];
 p_o[5]= alpha[17];
 p_o[6]= alpha[18];
 p_o[7]= alpha[19];

 p_k[1]= alpha[20];
 p_k[2]= alpha[21];
 p_k[3]= alpha[22];
 p_k[4]= alpha[23];
 p_k[5]= alpha[24];
 p_k[6]= alpha[25];
 p_k[7]= alpha[26];
 p_k[8]= alpha[27];
 p_k[9]= alpha[28];
 p_k[10]= alpha[29];
 p_k[11]= alpha[30];
 p_k[12]= alpha[31];
 p_k[13]= alpha[32];
 p_k[14]= alpha[33];
 p_k[15]= alpha[34];
 p_k[16]= alpha[35];

 p_m[1]= alpha[36];
 p_m[2]= alpha[37];
 p_m[3]= alpha[38];

 p_alpha= alpha[39];
		//for(int i=0;i<40;i++) System.out.printf(" a[%d] =%f\n", i,alpha[i]);
	}


	double SQR(double x){
		return x*x;
	}
	
	void CalcDistance(){
		for(int i=0; i<nAtom; i++ ){
			for(int j=i+1; j< nAtom;j++)
				if(isUpdate[i] || isUpdate[j]){
					r2[i][j] = 0;
					for(int k=0;k<3 ; k++) r2[i][j] += SQR(x[i][k] - x[j][k]);
					r2[j][i]=r2[i][j];
					r[i][j] = r[j][i] = Math.sqrt(r2[i][j]);
			}
		}
	}

	void CalcPairWiseInteraction(){
		VOO=VOH=VHH=0;
		for(int i=0; i<nO; i++ )
			for(int j=i+1; j< nO;j++){
				if(isUpdate[i] || isUpdate[j]){
					Vpairwise[i][j] = p_o[1]*Math.exp(-p_o[2] * r[i][j])
						 + p_o[3]*Math.exp(-p_o[4] * r[i][j])
						 + p_o[5]*Math.exp(-p_o[6] * SQR(r[i][j] - p_o[7]))
						 + Math.exp(-30.0*(r[i][j]-1.8));
					//System.out.printf("Come here %f \n",Vpairwise[i][j]);
				}
				VOO+=Vpairwise[i][j];
			}

		//------------ H - H
		for(int i=nO; i<nAtom; i++ ){
			for(int j=i+1; j< nAtom;j++){
				if(isUpdate[i] || isUpdate[j]){
					Vpairwise[i][j] = Math.exp(-100.0*(r[i][j]-1.1));
				}
				VHH+=Vpairwise[i][j];
			}
		}

		double h5tmp = Math.pow(1.0-p_h[5],2) / (Math.pow(1.0-p_h[5],2) + p_h[5]*p_h[5]);
		
		for(int i=0; i<nO; i++ ){
			for(int j=nO; j< nAtom;j++){
				if(isUpdate[i] || isUpdate[j]){
					Vpairwise[i][j] = p_h[1] * Math.pow(1 - h5tmp*Math.exp(-p_h[3]*(r[i][j]-p_h[2]))
							- (1-h5tmp)*Math.exp(-p_h[4]*(r[i][j]-p_h[2])), 2) - p_h[1]
							+ Math.exp(-70.0*(r[i][j]-0.3));
				}
				VOH+=Vpairwise[i][j];
			}
		}

		Vpairwise_total=0;
		for(int i=0; i<nAtom; i++ )
			for(int j= 0; j< nAtom;j++)
				Vpairwise_total+=Vpairwise[i][j];
	}	
	
	void CalcThreeBodyInteraction(){
		double r11, r21, r12, r22, t, dt, dt2,fcutoff,temp;
		for(int i=0;i<nO;i++){
		for(int j=nO;j< nAtom ; j++){
		for(int k=j+1; k< nAtom ; k++){
			if(isUpdate[i] || isUpdate[j] || isUpdate[k]){					
					r11 = r[i][j] - p_r0;
					r21 = r[i][k] - p_r0;
					r12 = SQR(r11);
					r22 = SQR(r21);
					if(Math.abs(p_m[1]*(r12+r22))<40){ // improved by chinh
						// theta & shifted theta
						t = Math.acos( (r2[i][j] + r2[i][k] - r2[j][k]) / (2.0*r[i][j]*r[i][k]));
						dt = t - p_theta0;
						dt2 = SQR(dt);

						fcutoff = Math.exp(-(p_m[1]*(r12+r22) +  p_m[2]*dt2 + p_m[3]*(r12+r22)*dt2));

						temp=0;

						temp = p_k[1]+ p_k[2]*(r11+r21) + p_k[3]*dt
						+ p_k[4]*(r12+r22) + p_k[5]*r11*r21 + p_k[6]*dt2 + p_k[7]*(r11+r21)*dt
						+ p_k[8]*(r12*r11+r22*r21) + p_k[9]*(r12*r21+r11*r22) + p_k[10]*dt*dt2 + p_k[11]*(r12+r22)*dt
						+ p_k[12]*r11*r21*dt + p_k[13]*(r11+r21)*dt2
						+ p_k[14]*(SQR(r12) + SQR(r22)) + p_k[15]*r12*r22 + p_k[16]*SQR(dt2);


						//cout<<i<<" "<<j<<" "<<k<<" "<<fcutoff<<" "<<temp<<endl;
						VHOH[i][j][k] =temp*fcutoff;
				}
			}
		}
		}
		}

		VHOH_total=0;
		for(int i=0;i<nO;i++)
		for(int j=nO;j< nAtom ; j++)
		for(int k=j+1; k< nAtom ; k++) VHOH_total=VHOH[i][j][k];

	}

	void CalcChargeChargeInteraction(){
		Vq_total=0;
		for(int i=0; i<nAtom; i++ ){
			for(int j=i+1; j< nAtom;j++){
				if(isUpdate[i] || isUpdate[j]){
					Vq[i][j] = q[i]*q[j] / r[i][j];
				}
				Vq_total+=Vq[i][j];
			}
		}
	}

	void CalcCutoffFunc(){
		for(int i=0; i<nO;i++) {
			for(int j=nO; j< nAtom; j++) if(isUpdate[i] || isUpdate[j]){
				Scd[i][j] = Scd[j][i] = 1 / (r[i][j]*(r2[i][j] + p_a[1]*Math.exp(-p_a[2]*r[i][j])));
			}

			for(int j=i+1;j< nO; j++) if(isUpdate[i] || isUpdate[j]){
				Scd[i][j] = Scd[j][i] = 1 / (r[i][j]*(r2[i][j] + p_b[1]*Math.exp(-p_b[2]*r[i][j])));
				Sdd[i][j] = Sdd[j][i] = 1 / (r[i][j]*(r2[i][j] + p_c[1]*Math.exp(-p_c[2]*r[i][j])));
			}
		}

		double[] rij=new double[3];

		for(int i=0; i<nO;i++)
		for(int j=i+1;j< nO; j++) if(isUpdate[i] || isUpdate[j]){
			for(int k=0;k<3;k++){
				rij[k] = x[j][k] - x[i][k];
			}
			for(int k=0;k<3;k++){
				for(int l=0;l<3;l++){
					if (k==l)
						T[3*i+k][3*j+l] = 1 - 3*rij[k]*rij[l] / r2[i][j];
					else
						T[3*i+k][3*j+l] = -3*rij[k]*rij[l] / r2[i][j];

					D[3*i+k][3*j+l]  = T[3*i+k][3*j+l] * Sdd[i][j];

					T[3*j+l][3*i+k] = T[3*i+k][3*j+l];
					D[3*j+l][3*i+k] = D[3*i+k][3*j+l];
				}
			}
		}
	}

	void CalcDipoleMoment(){
		
		// <part 4: matrix E>
	double[] tmp=new double[3];
	double[] Etmp=new double[3*nO];
	double[] rij=new double[3];

	for(int i=0; i< nO;i++){
		for(int k=0;k<3;k++) tmp[k] = 0;

		for(int j=0;j< i; j++) if(isUpdate[i] || isUpdate[j])
			for(int k=0;k<3;k++) {
				rij[k] = x[j][k] - x[i][k];
				tmp[k] += -2 * Scd[i][j] * rij[k];
			}
		

		for(int j=i+1;j< nO; j++) if(isUpdate[i] || isUpdate[j])
			for(int k=0;k<3;k++){
				rij[k] = x[j][k] - x[i][k];
				tmp[k] += -2 * Scd[i][j] * rij[k];
			}

		for(int j=nO;j< nAtom; j++) if(isUpdate[i] || isUpdate[j])	{
			for(int k=0;k<3;k++){
				rij[k] = x[j][k] - x[i][k];
				//cout << Scd[i][j] << endl;
				tmp[k] += 1 * Scd[i][j] * rij[k];
			}
		

		for(int k=0;k<3;k++){
			Etmp[3*i+k]=E[3*i+k] = -p_alpha*tmp[k];

		}
	}
	// </part 4: matrix E>
	}

	double[][] Dtmp=new double[3*nO][3*nO];

	// <part5: mu>
	for(int i=0; i< 3*nO;i++){
		for(int j=i;j< 3*nO; j++){
			if (i == j)
				Dtmp[i][j] = 1 + p_alpha*D[i][j];
			else
				Dtmp[i][j] = Dtmp[j][i]= p_alpha*D[i][j];
		}
	}
	//============== Here, we need to solve Dtmp*mu = Etmp for mu ==/
	System.out.printf(" mu before solving: "); 
	for(int i=0; i< 3*nO;i++){ System.out.printf("%f ", mu[i]);}
	System.out.println("");



	#if defined(OSS2_ITER)
		double ftol=1e-5;
		#ifdef OSS2_DEBUG
		cout<<"\n++++++++++++++ Solve by Iterative \n";
		#endif
		IR_Gauss_Seidel(Dtmp,mu,Etmp,3*nO,10000,ftol);
	#elif defined(OSS2_GSL)
		#ifdef OSS2_DEBUG
		cout<<"\n++++++++++++++ Solve by by GSL library \n";
		#endif
		GSL_LinearSolver(Dtmp,mu,Etmp,3*nO);
	#elif defined(OSS2_LAPACK)
		#ifdef OSS2_DEBUG
		cout<<"\n++++++++++++++ Solve by by Lapack \n";
		#endif
		Lapack_LinearSolver(Dtmp,mu,Etmp,3*nO);
	#else
		#ifdef OSS2_DEBUG
		cout<<"\n++++++++++++++ Solve by by Gaussian elimination \n";
		#endif
		LinearSolver(Dtmp,mu,Etmp,3*nO);
	#endif

	#ifdef OSS2_DEBUG
	cout<<" mu after solving: "; RP(i, 3*nO) cout <<setprecision(12)<< mu[i] << " ";cout<<endl;
	double residue=0;
	RP(i, 3*nO){
		double tmp = -Etmp[i];
		RP(j, 3*nO)	tmp += Dtmp[i][j] * mu[j];
		//cout<<setprecision(4)<<Dtmp[i][j]<<" ";
		residue=(fabs(tmp)>residue)?fabs(tmp):residue;
	}
	cout<<" Tolerance D*mu - E = "<<residue<<endl;
	#endif
	//===================

	// cout << endl;
	// </part5: mu>
	}
}

