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
import nqcmol.tools.MTools;

/**
 *
 * @author nqc
 */
public class OSS2Potential extends Potential{
	public OSS2Potential(){
		nativeUnit="Hartree";
		unit=nativeUnit;
		ParseParameters();
	};

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
	public boolean isValidSetup() {
		return true;
	}

	@Override
	public void setParam(double[] param) {
		assert alpha.length<param.length;
		System.arraycopy(param,0,alpha,0,alpha.length);
		ParseParameters();
	}

	@Override
	public void setParam(String filename) {
		if(!filename.isEmpty()){
		try {
			Scanner scanner = new Scanner(new File(filename));
			for(int i=0;i<alpha.length;i++){
				alpha[i]=scanner.nextDouble();
				scanner.nextLine();
			}
		} catch (FileNotFoundException ex) {
			Logger.getLogger(OSS2Potential.class.getName()).log(Level.SEVERE, null, ex);
		}
		}
		ParseParameters();
	}

	@Override
	public void setCluster(Cluster cluster_) {
		super.setCluster(cluster_);
		AllocatePrivateVariables(cluster_.getNAtoms());
		cluster.CorrectOrder();
	}

	@Override
	protected double Energy_(double[] _p){
		double energy=0;
		if(_p.length<=3) return 0;
		//System.out.print(" Size of vec = "+ Integer.toString(_p.length));
		AllocatePrivateVariables(_p.length/3);

		CalcDistanceAndCharge(_p,false);

		double VOO=OOInteraction(false);
		double VOH=OHInteraction(false);
		double VHH=HHInteraction(false);
		double VHOH=ThreeBodyInteraction(false);

		CalcCutoffFunc(false);
		CalcElectricalField();
		CalcDipoleMoment();

		double Vq=0;		
		double Vcd=0;
		double Vdd=0;		
		double Vsd=0;
		Vq=ChargeChargeInteraction(false);
		Vcd=ChargeDipoleInteraction(false);
		Vdd=DipoleDipoleInteraction(false);
		Vsd=SelfDipoleInteraction();


		double Vel=(Vdd + Vcd + Vq + Vsd) ;//* rescale;
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


		energy=ConvertUnit(energy,nativeUnit,unit);

		nEvals++;
		return energy;
	}

	@Override
	protected void Gradient_(double[] _p, double[] gradient) {
		double energy=0;
		//System.out.print(" Size of vec = "+ Integer.toString(_p.length));
		if(_p.length<=3) return ;
		
		AllocatePrivateVariables(_p.length/3);

		CalcDistanceAndCharge(_p,true);

		double VOO=OOInteraction(true);
		double VOH=OHInteraction(true);
		double VHH=HHInteraction(true);
		double VHOH=ThreeBodyInteraction(true);

		CalcCutoffFunc(true);
		CalcElectricalField();
		CalcDipoleMoment();

		double Vq=0;
		double Vcd=0;
		double Vdd=0;		
		double Vsd=0;
		Vq=ChargeChargeInteraction(true);
		Vcd=ChargeDipoleInteraction(true);
		Vdd=DipoleDipoleInteraction(true);
//		Vsd=SelfDipoleInteraction();
//
//		double Vel=(Vdd + Vcd + Vq + Vsd);// * rescale;
//		energy=VOO+VOH+VHH+VHOH + Vel;

		for(int i=0;i<nAtom;i++)
			for(int k=0;k<3;k++)
				gradient[i*3+k]=ConvertUnit(this.grad[i][k],nativeUnit,unit);


//		System.err.print("E : ");
//		for(int i=0;i<nO*3;i++)	System.err.printf( " %f ",E[i]);
//		System.err.println("");
//		System.err.printf( "VOO: %f  VOH: %f VHOH: %f VHH:%f \n",VOO,VOH,VHOH,VHH);
//		System.err.printf( "Vq: %f \n",Vq);
//		System.err.printf( "Self-dipole: %f \n" , Vsd );
//		System.err.printf( "Dipole-dipole interaction: %f \n" ,Vdd);
//		System.err.printf( "Charge-dipole interaction: %f \n" ,Vcd);
//		System.err.printf( "Vel: %f \n",Vel);
//		System.err.printf( "Energy: %f \n",energy);

		nEvals+=3;
	}

	

	//========================== implement OSS2
	int nAtom,nO;
	double[][] r,x,r2;
	double[] q,E,mu;
	double[][] Scd,Sdd,D,T;

	final double rescale=0.529177249;

	//for gradient
	double[][][] rv,dr;
	double[][] grad;
	double[][] dScd,dSdd;

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


	private void ParseParameters(){
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

	private void AllocatePrivateVariables(int nAtoms_){
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

			//if(true){ ///for gradient
				grad=new double[nAtom][3];
				rv=new double[nAtom][nAtom][3];
				dr=new double[nAtom][nAtom][3];
				dScd=new double[nAtom][nAtom];
				dSdd=new double[nAtom][nAtom];
				//clear gradient vector
				
			//}
		}

		for(int i=0;i<nAtom;i++)
			for(int k=0;k<3;k++)
					grad[i][k]=0;
	}

	void CalcDistanceAndCharge(double[] _p,boolean isGrad){
		for(int i=0;i<nAtom;i++) {
			x[i][0]= _p[i*3+0];
			x[i][1]= _p[i*3+1];
			x[i][2]= _p[i*3+2];
			//System.out.printf(" x[%d] = %f %f %f \n",i,x[i][0],x[i][1],x[i][2]);
		}
		
		for(int i=0; i<nAtom; i++ ){
			for(int j=i+1; j< nAtom;j++){
					r2[i][j] = 0;
					for(int k=0;k<3 ; k++) r2[i][j] += SQR(x[i][k] - x[j][k]);
					r2[j][i]=r2[i][j];
					r[i][j] = r[j][i] = Math.sqrt(r2[i][j]);


					if(isGrad) ///for gradient
					for(int k=0;k<3;k++){
							rv[i][j][k] = x[i][k] - x[j][k];
							rv[j][i][k] = -rv[i][j][k];
							dr[i][j][k] = rv[i][j][k]/r[i][j];
							dr[j][i][k] =-dr[i][j][k];
						}
			}
		}

		for(int i=0; i<nAtom; i++ ){
				q[i] = ((i<nO)? -2 : 1);
		}
	}

	double OOInteraction(boolean isGrad){
		double VOO=0;
		for(int i=0; i<nO; i++ )
			for(int j=i+1; j< nO;j++){				
					double pt1 = p_o[1]*Math.exp(-p_o[2] * r[i][j]);
					double pt2 = p_o[3]*Math.exp(-p_o[4] * r[i][j]);
					double pt3 = p_o[5]*Math.exp(-p_o[6] * SQR(r[i][j] - p_o[7]));

					VOO += pt1 + pt2 + pt3;

					double pt4 =  Math.exp(-30.0*(r[i][j]-1.8));
					VOO += pt4;

					if(isGrad){ /// for gradient
						// dp = pt1' + pt2' + pt3'
						double dp = (-p_o[2]*pt1 - p_o[4]*pt2 - p_o[6]*2.0*(r[i][j] - p_o[7])*pt3);
						dp+=-30.0*pt4; //prevent unphysically short O-O

						for(int k=0;k<3;k++){
							double dpk = dp * dr[i][j][k];
							grad[i][k] += dpk;
							grad[j][k] -= dpk;
						}
					}
					//System.out.printf("Come here %f \n",VOO);
			}
		return VOO;
	}

	double HHInteraction(boolean isGrad){
		double VHH=0;
		for(int i=nO; i<nAtom; i++ ){
			for(int j=i+1; j< nAtom;j++){
				double temp1=Math.exp(-100.0*(r[i][j]-1.1));//prevent unphysically short O-H
				VHH+=temp1;

				if(isGrad){
				for(int k=0;k<3;k++){
					double temp2=-100.0*temp1*dr[i][j][k]; //prevent unphysically short O-H
					grad[i][k]+=  temp2;
					grad[j][k]+= -temp2;
					}
				}
			}
		}
		return VHH;
	}

	double OHInteraction(boolean isGrad){
		double VOH=0;
		double h5tmp = Math.pow(1.0-p_h[5],2) / (Math.pow(1.0-p_h[5],2) + p_h[5]*p_h[5]);

		for(int i=0; i<nO; i++ ){
			for(int j=nO; j< nAtom;j++){
				//VOH+= p_h[1] * Math.pow(1 - h5tmp*Math.exp(-p_h[3]*(r[i][j]-p_h[2]))
				//			- (1-h5tmp)*Math.exp(-p_h[4]*(r[i][j]-p_h[2])), 2) - p_h[1]
				//			+ Math.exp(-70.0*(r[i][j]-0.3));

				double pt1 = h5tmp*Math.exp(-p_h[3]*(r[i][j]-p_h[2]));
				double pt2 = (1-h5tmp)*Math.exp(-p_h[4]*(r[i][j]-p_h[2]));

				VOH += p_h[1] * SQR(1 - pt1 - pt2) - p_h[1];

				double pt5=Math.exp(-70.0*(r[i][j]-0.3)); //prevent unphysically short 
				VOH += pt5;//prevent unphysically short O-H

				if(isGrad){
					// pt3 = pt1' + pt2'
					double pt3= -p_h[3] * pt1 - p_h[4]*pt2;
					double dp = -p_h[1] * 2 * (1-pt1-pt2) * pt3;
					dp+=-70.0*pt5; //prevent unphysically short
					for(int k=0;k<3;k++){
						double dpk = dp * dr[i][j][k];
						grad[i][k] += dpk;
						grad[j][k] -= dpk;
					}
				}

			}
		}
		return VOH;
	}

	double ThreeBodyInteraction(boolean isGrad){
		double VHOH=0;
		double r11, r21, r12, r22, theta, cost, dt, dt2,fcutoff,temp;
		for(int i=0;i<nO;i++){
			for(int j=nO;j< nAtom ; j++){
				for(int k=j+1; k< nAtom ; k++){
							r11 = r[i][j] - p_r0;
							r21 = r[i][k] - p_r0;
							r12 = SQR(r11);
							r22 = SQR(r21);
							if(Math.abs(p_m[1]*(r12+r22))<40){ // improved by chinh
								// theta & shifted theta
								cost = (r2[i][j] + r2[i][k] - r2[j][k]) / (2.0*r[i][j]*r[i][k]);
								theta = Math.acos( cost);
								dt = theta - p_theta0;
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

								if(isGrad){ ///for gradient
									// fcutoff'
									double dfcutoff_r1=-2.0*(p_m[1]*r11 + p_m[3]*r11*dt2)*fcutoff;
									double dfcutoff_r2=-2.0*(p_m[1]*r21 + p_m[3]*r21*dt2)*fcutoff;
									double dfcutoff_t =-2.0*(p_m[2]+p_m[3]*(r12+r22))*dt*fcutoff;

									// temp'
									double dtemp_r1 = p_k[2] + p_k[4]*2.0*r11 + p_k[5]*r21 + p_k[7]*dt
									+ p_k[8]*3.0*r12 + p_k[9]*(2.0*r11*r21+r22) + p_k[11]*2.0*r11*dt
									+ p_k[12]*r21*dt + p_k[13]*dt2 +p_k[14]*4.0*r11*r12 + p_k[15]*2.0*r11*r22;
									double dtemp_r2 = p_k[2] + p_k[4]*2.0*r21 + p_k[5]*r11 + p_k[7]*dt
									+ p_k[8]*3.0*r22 + p_k[9]*(2.0*r21*r11+r12) + p_k[11]*2.0*r21*dt
									+ p_k[12]*r11*dt + p_k[13]*dt2 +p_k[14]*4.0*r21*r22 + p_k[15]*2.0*r21*r12;
									double dtemp_t  = p_k[3] + p_k[6]*2.0*dt + p_k[7]*(r11+r21) + p_k[10]*3.0*dt2 + p_k[11]*(r12+r22)
									+ p_k[12]*r11*r21 + p_k[13]*(r11+r21)*2.0*dt + p_k[16]*4.0*dt*dt2;

									// VHOH'
									double dVHOH_r1 = fcutoff*dtemp_r1 + dfcutoff_r1*temp;
									double dVHOH_r2 = fcutoff*dtemp_r2 + dfcutoff_r2*temp;
									double dVHOH_t  = fcutoff*dtemp_t  + dfcutoff_t *temp;

									// r11', r21' and theta'
									// r11' = dr[i][j] && r21' = dr[i][k]
									double sint = Math.sqrt(1-SQR(cost));

									for(int t=0;t<3;t++){
										double dt_rj = (-1.0/sint) * (cost*dr[i][j][t] - dr[i][k][t])/r[i][j];
										double dt_rk = (-1.0/sint) * (cost*dr[i][k][t] - dr[i][j][t])/r[i][k];
										double dt_ri = -(dt_rj + dt_rk);

										grad[i][t]+= dVHOH_r1*dr[i][j][t] + dVHOH_r2*dr[i][k][t] + dVHOH_t * dt_ri;
										grad[j][t]+= dVHOH_r1*dr[j][i][t]  						 + dVHOH_t * dt_rj;
										grad[k][t]+=                        dVHOH_r2*dr[k][i][t] + dVHOH_t * dt_rk;
									}
								}
						}
				}
			}
		}
		return VHOH;
	}

	double ChargeChargeInteraction(boolean isGrad){
		double Vq=0;
		for(int i=0; i<nAtom; i++ ){
			for(int j=i+1; j< nAtom;j++){
				double temp1=q[i]*q[j] / r[i][j] *rescale;
				Vq+= temp1;
				if(isGrad){ //for gradient
					for(int k=0;k<3;k++){
						double temp2=(-q[i]*q[j]/r2[i][j])*dr[i][j][k]*rescale;
						grad[i][k]+= temp2;
						grad[j][k]+=-temp2;
					}
				}
			}
		}
		return Vq;
	}

	void CalcCutoffFunc(boolean isGrad){
		for(int i=0; i<nO;i++) {
			for(int j=nO; j< nAtom; j++) {
				Scd[i][j] = Scd[j][i] = 1 / (r[i][j]*(r2[i][j] + p_a[1]*Math.exp(-p_a[2]*r[i][j])));

				if(isGrad){ //for gradient
					dScd[i][j]= dScd[j][i] = -SQR(Scd[i][j])*(3.0*r2[i][j] + p_a[1]*Math.exp(-p_a[2]*r[i][j])*(1.0-p_a[2]*r[i][j]));
				}
			}

			for(int j=i+1;j< nO; j++) {
				Scd[i][j] = Scd[j][i] = 1 / (r[i][j]*(r2[i][j] + p_b[1]*Math.exp(-p_b[2]*r[i][j])));
				Sdd[i][j] = Sdd[j][i] = 1 / (r[i][j]*(r2[i][j] + p_c[1]*Math.exp(-p_c[2]*r[i][j])));

				if(isGrad){ //for gradient
					dScd[i][j] = dScd[j][i] =-SQR(Scd[i][j])*(3.0*r2[i][j] + p_b[1]*Math.exp(-p_b[2]*r[i][j])*(1.0-p_b[2]*r[i][j]));
					dSdd[i][j] = dSdd[j][i] =-SQR(Sdd[i][j])*(3.0*r2[i][j] + p_c[1]*Math.exp(-p_c[2]*r[i][j])*(1.0-p_c[2]*r[i][j]));
				}
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

		MTools.LinearSolver_apache(Dtmp,Etmp,mu);

//		System.out.print(" mu after solving: ");
//		for(int i=0; i< 3*nO;i++){ System.out.printf("%f ", mu[i]);}
//		System.out.println("");
	}

	double ChargeDipoleInteraction(boolean isGrad){
		double Vcd = 0;
		for(int i=0;i< 3*nO; i++){
			Vcd -= mu[i] * E[i] / p_alpha *rescale;
		}

		if(isGrad){
				for(int i=0;i<nAtom;i++){
					for(int j=0;j<nO;j++){
						if(i!=j){
							for(int t=0;t<3;t++){
								double tmp = 0;
								for(int k=0;k<3;k++) tmp += mu[3*j+k] * rv[i][j][k];

								tmp =q[i]*tmp*dScd[i][j]*dr[i][j][t];
								tmp+=q[i]*mu[3*j+t]*Scd[i][j];

								tmp*=rescale;

								grad[i][t]+=tmp;
								grad[j][t]-=tmp;
							}
						}
				}
			}
		}
		return Vcd;
	}

	double DipoleDipoleInteraction(boolean isGrad){
		double Vdd = 0;
		for(int i=0;i< 3*nO; i++){
			for(int j=i+1;j< 3*nO; j++){
				Vdd += mu[i]*D[i][j]*mu[j] *rescale;
			}
		}

		if(isGrad){
			for(int i=0;i<nO;i++){
				for(int j=0;j<nO;j++){
					if(i!=j){
						double pt1 = 0;
						for(int k=0;k<3;k++)
							for(int l=0;l<3;l++)
								pt1 += mu[3*i+k]*T[3*i+k][3*j+l]*mu[3*j+l];

						double pt2 = 0;
						double pt3 = 0;

						for(int k=0;k<3;k++) {
							pt2 += mu[i*3+k] * dr[i][j][k];
							pt3 += mu[j*3+k] * dr[i][j][k];
						}

						for(int t=0;t<3;t++){
							double temp = pt1*dSdd[i][j]*dr[i][j][t];
							temp += 3.0*(2.0*dr[i][j][t]*pt2*pt3
							- mu[j*3+t]*pt2
							- mu[i*3+t]*pt3)/r[i][j]*Sdd[i][j];

							grad[i][t]+=temp*rescale;
						}
					}
				}
			}
		}
		return Vdd;
	}

	double SelfDipoleInteraction(){
		double Vsd = 0;
		for(int i=0;i< 3*nO; i++) Vsd += mu[i]*mu[i];
		Vsd /= 2*p_alpha;
		Vsd*=rescale;
		return Vsd;
	}

	//================== for gradient
}
