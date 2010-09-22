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
import nqcmol.cluster.Cluster;
import nqcmol.tools.MTools;

/**
 *
 * @author nqc
 */
public class OSS2RPotential extends Potential{
	public OSS2RPotential(){
		nativeUnit="Hartree";
		unit=nativeUnit;
		ParseParameters();
	};

	@Override
	public String getEquation(){
		String equation="OSS2R";
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
			Logger.getLogger(OSS2RPotential.class.getName()).log(Level.SEVERE, null, ex);
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
		double VHOH=BendingAngleInteraction(false);

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
		double VHOH=BendingAngleInteraction(true);

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
	double[] p_h=new double[3];
	double[] p_o=new double[5];
	double[] p_k=new double[6];
	double[] p_m=new double[4];
	double p_alpha;
    double p_qO,p_qH;

	double[] alpha={0.96140000, 1.81733811,
                    1.73089135, 0.29575998, 336.05523682, 2.75564718, 0.20642018, 0.00000000,
                    0.00618378, 2.64253736,
                    40.48587418, 1.32909548,
                    -0.04034958, 0.16488265, -0.02451205, -0.37863523, -0.35223219,
                    6.25000000, 0.00000000, 6.25000000,
                    1.44400000, -2, 1};

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
		

		 p_o[1]= alpha[10];
		 p_o[2]= alpha[11];
		 

		 p_k[1]= alpha[12];
		 p_k[2]= alpha[13];
		 p_k[3]= alpha[14];
		 p_k[4]= alpha[15];
		 p_k[5]= alpha[16];
		
		 p_m[1]= alpha[17];
		 p_m[2]= alpha[18];
		 p_m[3]= alpha[19];

		 p_alpha= alpha[20];
         p_qO= alpha[21];
         p_qH = alpha[22];
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
				q[i] = ((i<nO)? p_qO : p_qH);
		}
	}

	double OOInteraction(boolean isGrad){
		double VOO=0;
		for(int i=0; i<nO; i++ )
			for(int j=i+1; j< nO;j++){
                    double rij6 = Math.pow(r2[i][j],3);
                    double rij12 = Math.pow(rij6,2);

                    double pt2 =  p_o[1]/rij12 ;
                    double pt3 =- p_o[2]/rij6;
					double pt4 =  Math.exp(-30.0*(r[i][j]-1.8));
					VOO += pt2 + pt3 + pt4;

					if(isGrad){ /// for gradient
						// dp = pt1' + pt2' + pt3'
						double dp = -6*(2*pt2+ pt3)/r[i][j];
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
		for(int i=0; i<nO; i++ ){
			for(int j=nO; j< nAtom;j++){
				//VOH+= p_h[1] * Math.pow(1 - h5tmp*Math.exp(-p_h[3]*(r[i][j]-p_h[2]))
				//			- (1-h5tmp)*Math.exp(-p_h[4]*(r[i][j]-p_h[2])), 2) - p_h[1]
				//			+ Math.exp(-70.0*(r[i][j]-0.3));

				double pt1 = Math.exp(-p_h[2]*(r[i][j]-p_r0));

				VOH += p_h[1] *Math.pow(1-pt1,2); //p_De*pow(1.0 - exp(-p_Be*(r[i][j]-p_r0)),2.0);

				double pt5=Math.exp(-70.0*(r[i][j]-0.3)); //prevent unphysically short 
				VOH += pt5;//prevent unphysically short O-H

				if(isGrad){
					// pt3 = pt1' + pt2'
					double dp= p_h[1] * 2 * (1-pt1) * p_h[2]* pt1;
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

	double BendingAngleInteraction(boolean isGrad){
		double Ebend=0;
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

								temp = p_k[1]+ p_k[2]*dt + p_k[3]*dt2 + p_k[4]*dt*dt2 + p_k[5]*SQR(dt2);


								//cout<<i<<" "<<j<<" "<<k<<" "<<fcutoff<<" "<<temp);
								Ebend+=temp*fcutoff;

								if(isGrad){ ///for gradient
									// fcutoff'
									//double dfcutoff_r1=-2.0*(p_m[1]*r11 + p_m[3]*r11*dt2)*fcutoff;
									//double dfcutoff_r2=-2.0*(p_m[1]*r21 + p_m[3]*r21*dt2)*fcutoff;
									double dfcutoff_t =-2.0*(p_m[2]+p_m[3]*(r12+r22))*dt*fcutoff;

									// temp'									
									double dtemp_t  = p_k[2] + p_k[3]*2.0*dt +p_k[4]*3.0*dt2 + p_k[5]*4.0*dt*dt2;

									// VHOH'
									double dEbend_t  = fcutoff*dtemp_t  + dfcutoff_t *temp;

									// r11', r21' and theta'
									// r11' = dr[i][j] && r21' = dr[i][k]
									double sint = Math.sqrt(1-SQR(cost));

									for(int t=0;t<3;t++){
										double dt_rj = (-1.0/sint) * (cost*dr[i][j][t] - dr[i][k][t])/r[i][j];
										double dt_rk = (-1.0/sint) * (cost*dr[i][k][t] - dr[i][j][t])/r[i][k];
										double dt_ri = -(dt_rj + dt_rk);

										grad[i][t]+= dEbend_t * dt_ri;
										grad[j][t]+= dEbend_t * dt_rj;
										grad[k][t]+= dEbend_t * dt_rk;
									}
								}
                            }
				}
			}
		}
		return Ebend;
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
