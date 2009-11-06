/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.potential;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.Cluster;
import org.apache.commons.math.linear.*;

/**
 *
 * @author nqc
 */
public class OSS2Potential extends Potential{
	public OSS2Potential(){};

	@Override
	public String getEquation(){
		String equation="OSS2";
		return equation;
	}

	
	@Override
	public boolean HasAnalyticalGradients() {
		return true;
	}

	@Override
	public String getUnit() {
		return "Hartree";
	}

	@Override
	public boolean isValidSetup() {
		return true;
	}

	@Override
	public void setParam(double[] param) {
		assert alpha.length<param.length;
		System.arraycopy(param,0,alpha,0,alpha.length);
	}

	@Override
	public void setParam(String filename) {
		try {
			Scanner scanner = new Scanner(new File(filename));
			for(int i=0;i<alpha.length;i++)
				alpha[i]=scanner.nextDouble();
		} catch (FileNotFoundException ex) {
			Logger.getLogger(OSS2Potential.class.getName()).log(Level.SEVERE, null, ex);
		}

	}



	@Override
	public void setCluster(Cluster cluster_) {
		super.setCluster(cluster_);
		InitializePrivateVariables(cluster_.getNAtoms());
	}

	@Override
	protected double Energy_(double[] _p){
		double energy=0;
		System.out.print(" Size of vec = "+ Integer.toString(_p.length));
		InitializePrivateVariables(_p.length/3);

		//System.out.print(" Size of vec = "+ Integer.toString(coords3d.getSize()));

		for(int i=0;i<nAtom;i++) {
			x[i][0]= _p[i*3+0];
			x[i][1]= _p[i*3+1];
			x[i][2]= _p[i*3+2];
			// cout<<i<<" d "<<x[i][0]<<" - "<<x[i][1]<<" -  "<<x[i][2]);
		}


		CalcDistance();


		double VOO=OOInteraction();
		double VOH=OHInteraction();
		double VHH=HHInteraction();
		double VHOH=ThreeBodyInteraction();

		CalcCutoffFunc();
		CalcElectricalField();
		CalcDipoleMoment();

		double Vq=ChargeChargeInteraction();
		double Vcd=ChargeDipoleInteraction();
		double Vdd=DipoleDipoleInteraction();
		double Vsd=SelfDipoleInteraction();


		double Vel=(Vdd + Vcd + Vq + Vsd) * rescale;
		energy=VOO+VOH+VHH+VHOH + Vel;

//		System.out.print("E : ");
//		for(int i=0;i<nO*3;i++)	System.out.printf( " %f ",E[i]);
//		System.out.println("");
//		System.out.printf( "VOO: %f  VOH: %f VHOH: %f VHH:%f \n",VOO,VOH,VHOH,VHH);
//		System.out.printf( "Vq: %f \n",Vq);
//		System.out.printf( "Self-dipole: %f \n" , Vsd );
//		System.out.printf( "Dipole-dipole interaction: %f \n" ,Vdd);
//		System.out.printf( "Charge-dipole interaction: %f \n" ,Vcd);
//		System.out.printf( "Vel: %f \n",Vel);
//		System.out.printf( "Energy: %f \n",energy);

		neval++;
		return energy;
	}

	@Override
	protected void Gradient_(double[] _p, double[] grad) {
		nAtom=_p.length/3;
		for (int i=0; i < nAtom ; i++) {
//			energyGradient.setElement(i,
			for(int k=0;k<3;k++)	grad[i*3+k]=0;
			for(int j=0;j<nAtom;j++)	if(i!=j){
				double rij2=  Math.pow( (_p[i*3+0]-_p[j*3+0] ),2)
							+ Math.pow( (_p[i*3+1]-_p[j*3+1] ),2)
							+ Math.pow( (_p[i*3+2]-_p[j*3+2] ),2);
				double rij6 = Math.pow(rij2,3);

				double rij14= rij6*rij6*rij2;
				double f=24.0*(rij6-2.0)/(rij14);
				for(int k=0;k<3;k++)	grad[i*3+k]+=f*(_p[i*3+k]-_p[j*3+k]);
			}

		}
	}

	private void InitializePrivateVariables(int nAtoms_){
		if(nAtom!=nAtoms_){
			nAtom=nAtoms_;
			nO=(nAtom%3 == 2)?(nAtom/3+1):(nAtom/3);
			x=new double[nAtom][3];
			q=new double[nAtom];
			r=new double[nAtom][nAtom];
			r2=new double[nAtom][nAtom];

			Scd=new double[nAtom][nAtom];
			Sdd=new double[nO][nO];
			D=new double[nO*3][nO*3];
			T=new double[nO*3][nO*3];
			E=new double[nO*3];
			mu=new double[nO*3];

			for(int i=0; i<nAtom; i++ ){
				q[i] = ((i<nO)? -2 : 1);
			}
		}
	}

	//========================== implement OSS2
	int nAtom,nO;
	double[] q;
	double[][] r,x,r2;

	double[] E,mu;

	double[][] Scd,Sdd,D,T;

	final double rescale=0.529177249;

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

	double[] alpha={0.9614,	1.81733805,	1.73089134,	0.29575999,	332.2852922, 2.84683483,	1.0956864,	0.00000204,	0.00575558,	2.6425374,	1.1274385, 3.3255655,	0.07847734,	40.4858734,	1.3290955,	-41.7260608, 1.3509491, 0.0629396,	6.77909903,	1.8117836,	-0.0420073,	0.16488265,	-0.02509795, -0.37814525, -0.31667595, -0.0114672, 0.06061075, 0.528281, 1.1617627, 0.0765232, -0.21208835, -0.1025385, -0.0762216, -0.228694, 0, -0.029092, 6.25, 0.00377233,	6.25, 1.444 };


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
			for(int j=i+1; j< nAtom;j++){
					r2[i][j] = 0;
					for(int k=0;k<3 ; k++) r2[i][j] += SQR(x[i][k] - x[j][k]);
					r2[j][i]=r2[i][j];
					r[i][j] = r[j][i] = Math.sqrt(r2[i][j]);
			}
		}
	}

	double OOInteraction(){
		double VOO=0;
		for(int i=0; i<nO; i++ )
			for(int j=i+1; j< nO;j++){
				{
					VOO+= p_o[1]*Math.exp(-p_o[2] * r[i][j])
						 + p_o[3]*Math.exp(-p_o[4] * r[i][j])
						 + p_o[5]*Math.exp(-p_o[6] * SQR(r[i][j] - p_o[7]))
						 + Math.exp(-30.0*(r[i][j]-1.8));
					//System.out.printf("Come here %f \n",Vpairwise[i][j]);
				}
			}
		return VOO;
	}

	double HHInteraction(){
		double VHH=0;
		for(int i=nO; i<nAtom; i++ ){
			for(int j=i+1; j< nAtom;j++){
				VHH+= Math.exp(-100.0*(r[i][j]-1.1));
			}
		}
		return VHH;
	}

	double OHInteraction(){
		double VOH=0;
		double h5tmp = Math.pow(1.0-p_h[5],2) / (Math.pow(1.0-p_h[5],2) + p_h[5]*p_h[5]);

		for(int i=0; i<nO; i++ ){
			for(int j=nO; j< nAtom;j++){
				VOH+= p_h[1] * Math.pow(1 - h5tmp*Math.exp(-p_h[3]*(r[i][j]-p_h[2]))
							- (1-h5tmp)*Math.exp(-p_h[4]*(r[i][j]-p_h[2])), 2) - p_h[1]
							+ Math.exp(-70.0*(r[i][j]-0.3));

			}
		}
		return VOH;
	}

	double ThreeBodyInteraction(){
		double VHOH=0;
		double r11, r21, r12, r22, t, dt, dt2,fcutoff,temp;
		for(int i=0;i<nO;i++){
			for(int j=nO;j< nAtom ; j++){
				for(int k=j+1; k< nAtom ; k++){
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


								//cout<<i<<" "<<j<<" "<<k<<" "<<fcutoff<<" "<<temp);
								VHOH+=temp*fcutoff;
						}
				}
			}
		}
		return VHOH;
	}

	double ChargeChargeInteraction(){
		double Vq=0;
		for(int i=0; i<nAtom; i++ ){
			for(int j=i+1; j< nAtom;j++){
				Vq+= q[i]*q[j] / r[i][j];
			}
		}
		return Vq;
	}

	void CalcCutoffFunc(){
		for(int i=0; i<nO;i++) {
			for(int j=nO; j< nAtom; j++) {
				Scd[i][j] = Scd[j][i] = 1 / (r[i][j]*(r2[i][j] + p_a[1]*Math.exp(-p_a[2]*r[i][j])));
			}

			for(int j=i+1;j< nO; j++) {
				Scd[i][j] = Scd[j][i] = 1 / (r[i][j]*(r2[i][j] + p_b[1]*Math.exp(-p_b[2]*r[i][j])));
				Sdd[i][j] = Sdd[j][i] = 1 / (r[i][j]*(r2[i][j] + p_c[1]*Math.exp(-p_c[2]*r[i][j])));
			}
		}

		double[] rij=new double[3];

		for(int i=0; i<nO;i++)
		for(int j=i+1;j< nO; j++) {
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

	void CalcElectricalField(){
		// <part 4: matrix E>
		double[] tmp=new double[3];

		for(int i=0; i< nO;i++){
			for(int k=0;k<3;k++) tmp[k] = 0;

			for(int j=0;j< nAtom; j++) if(i!=j)
				for(int k=0;k<3;k++) {
					double rijk = x[j][k] - x[i][k];
					tmp[k] += q[j] * Scd[i][j] * rijk;
				}

			for(int k=0;k<3;k++){
				E[3*i+k] = -p_alpha*tmp[k];
			}
		}
	}

	void CalcDipoleMoment(){
		double[][] Dtmp=new double[3*nO][3*nO];
		double[] Etmp=new double[3*nO];


		for(int i=0; i< 3*nO;i++){
			for(int j=i;j< 3*nO; j++){
				if(i==j)
					Dtmp[i][j] = Dtmp[j][i]=  ( 1 + p_alpha*D[i][j]);
				else
					Dtmp[i][j] = Dtmp[j][i]=  ( p_alpha*D[i][j]);
			}
			Etmp[i]=E[i];
		}


		//============== Here, we need to solve Dtmp*mu = Etmp for mu ==/
//		System.out.print(" mu before solving: ");
//		for(int i=0; i< 3*nO;i++){ System.out.printf("%f ", mu[i]);}
//		System.out.println("");

		RealMatrix A =new Array2DRowRealMatrix(Dtmp,false);
		DecompositionSolver solver = new LUDecompositionImpl(A).getSolver();
        RealVector b = new ArrayRealVector(Etmp,false);
		RealVector solution = solver.solve(b);

		mu=solution.getData();

//		System.out.print(" mu after solving: ");
//		for(int i=0; i< 3*nO;i++){ System.out.printf("%f ", mu[i]);}
//		System.out.println("");
	}

	double ChargeDipoleInteraction(){
		double Vcd = 0;
		for(int i=0;i< 3*nO; i++){
			Vcd -= mu[i] * E[i] / p_alpha;
		}
		return Vcd;
	}

	double DipoleDipoleInteraction(){
		double Vdd = 0;
		for(int i=0;i< 3*nO; i++){
			for(int j=i+1;j< 3*nO; j++){
				Vdd += mu[i]*D[i][j]*mu[j];
			}
		}
		return Vdd;
	}

	double SelfDipoleInteraction(){
		double Vsd = 0;
		for(int i=0;i< 3*nO; i++) Vsd += mu[i]*mu[i];
			Vsd /= 2*p_alpha;
		return Vsd;
	}
}
