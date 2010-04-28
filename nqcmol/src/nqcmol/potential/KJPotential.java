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
public class KJPotential extends Potential{
	public KJPotential(){
		nativeUnit="kcal/mol";
		unit=nativeUnit;
		ParseParameters();
	};

	@Override
	public String getEquation(){
		String equation="KJ";
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
			Logger.getLogger(KJPotential.class.getName()).log(Level.SEVERE, null, ex);
		}
		}
		ParseParameters();
	}

    private void ParseParameters(){
		p_r0= alpha[0];
		p_theta0= alpha[1];
		p_rSigma= alpha[2];
		p_q= alpha[3];
		//for(int i=0;i<alpha.length;i++) System.out.printf(" a[%d] =%f\n", i,alpha[i]);
	}


	@Override
	public void setCluster(Cluster cluster_) {
		super.setCluster(cluster_);
		//AllocatePrivateVariables(cluster_.getNAtoms());
		cluster.CorrectOrder();
        nO=cluster.getNonHydrogenNum();
		for(int i=0;i<nO;i++){
			int count=0;
			for(int j=nO;j<cluster.getNAtoms();j++)
				if((cluster.distance(i,j)<1.3)&&(nO+2*i+count<cluster.getNAtoms()-1)){
					cluster.SwapAtom(nO+2*i+count,j);
					count++;
				}
		}
	}

	@Override
	protected double Energy_(double[] _p){
		double energy=0;
		if(_p.length<=3) return 0;

        AddSigmaCharge(1);
		//System.out.print(" Size of vec = "+ Integer.toString(_p.length));
//		AllocatePrivateVariables(_p.length/3);
//
//		CalcDistanceAndCharge(_p,false);
//
//		double VOO=OOInteraction(false);
//		double VOH=OHInteraction(false);
//		double VHH=HHInteraction(false);
		double V2body=TwoBodyInteraction(false);
//
//		CalcCutoffFunc(false);
//		CalcElectricalField();
//		CalcDipoleMoment();
//
		double Vq=0;
//		double Vcd=0;
//		double Vdd=0;
//		double Vsd=0;
		Vq=ChargeChargeInteraction(false);
//		Vcd=ChargeDipoleInteraction(false);
//		Vdd=DipoleDipoleInteraction(false);
//		Vsd=SelfDipoleInteraction();
//
//
//		double Vel=(Vdd + Vcd + Vq + Vsd) ;//* rescale;
		energy=V2body + Vq;

//		System.out.print("E : ");
//		for(int i=0;i<nO*3;i++)	System.out.printf( " %f ",E[i]);
//		System.out.println("");
		System.out.printf( "VOO: %f \n",V2body);
		System.out.printf( "Vq: %f \n",Vq);
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
		
//		AllocatePrivateVariables(_p.length/3);
//
//		CalcDistanceAndCharge(_p,true);
//
//		double VOO=OOInteraction(true);
//		double VOH=OHInteraction(true);
//		double VHH=HHInteraction(true);
//		double VHOH=ThreeBodyInteraction(true);
//
//		CalcCutoffFunc(true);
//		CalcElectricalField();
//		CalcDipoleMoment();
//
//		double Vq=0;
//		double Vcd=0;
//		double Vdd=0;
//		double Vsd=0;
//		Vq=ChargeChargeInteraction(true);
//		Vcd=ChargeDipoleInteraction(true);
//		Vdd=DipoleDipoleInteraction(true);
////		Vsd=SelfDipoleInteraction();
////
////		double Vel=(Vdd + Vcd + Vq + Vsd);// * rescale;
////		energy=VOO+VOH+VHH+VHOH + Vel;
//
//		for(int i=0;i<nAtom;i++)
//			for(int k=0;k<3;k++)
//				gradient[i*3+k]=ConvertUnit(this.grad[i][k],nativeUnit,unit);


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

	

	//========================== implement KJ
    int nAtoms,nO;
    double[][] grandCoords;
    double[][] xO;
    double[][] xH1,xH2;
    double[][] xS;
    double[] q;
    double[][] E;

    final double rescale=0.529177249*627.509;

    double p_r0=0.957;
	double p_theta0=104.5;
    double p_rSigma=0.138;
    double p_q=0.6228;
    double p_epsilon=0.42;
    double p_sigma=3.17;


    double[] alpha={p_r0,p_theta0,p_rSigma,p_q,p_epsilon,p_sigma};

    void AddSigmaCharge(int debug){
        cluster.Write(System.out, "xyz");
        xO=new double[cluster.getNonHydrogenNum()][4];
        xH1=new double[cluster.getNonHydrogenNum()][4];
        xH2=new double[cluster.getNonHydrogenNum()][4];
        xS=new double[cluster.getNonHydrogenNum()][4];

        double[] x=cluster.getCoords();
        for(int i=0;i<nO;i++){
            //assign coords and charge of oxygen
            for(int j=0;j<3;j++)    xO[i][j]=x[i*3+j];
            xO[i][3]=2.0*p_q;

            int iH1=nO+i*2;
            int iH2=nO+i*2+1;

            //calculate OH1 and OH2 unit vectors
            double[] u1=new double[3];  for(int j=0;j<3;j++) u1[j]=x[iH1*3+j]-xO[i][j];  MTools.NORMALIZE(u1);
            double[] u2=new double[3];  for(int j=0;j<3;j++) u2[j]=x[iH2*3+j]-xO[i][j];  MTools.NORMALIZE(u2);


            double c=MTools.DOTPRODUCT(u1, u2);//angle of H-O-H

            double cp=Math.cos(Math.toRadians(p_theta0)); //angle of H-O-H expected in KJ

            //calculate new OH1, OH2 unit vectors in KJ
            double[] p1=new double[3]; 
            double[] p2=new double[3];

            double a1=0.5*(Math.sqrt((1+cp)/(1+c)) + Math.sqrt((1-cp)/(1-c)) );
            double a2=0.5*(Math.sqrt((1+cp)/(1+c)) - Math.sqrt((1-cp)/(1-c)) );

            MTools.VEC_PLUS_VEC(p1, u1, u2, a1,  a2);

            a1=0.5*(Math.sqrt((1+cp)/(1+c)) - Math.sqrt((1-cp)/(1-c)) );
            a2=0.5*(Math.sqrt((1+cp)/(1+c)) + Math.sqrt((1-cp)/(1-c)) );
            MTools.VEC_PLUS_VEC(p2, u1, u2, a1, a2);

            //calculate direction of O-sigma, aka symetric axis of water molecules
            double[] s=new double[3];
            MTools.VEC_PLUS_VEC(s,u1,u2,1,1);
            MTools.NORMALIZE(s);

            //calculate coords of H1 and H2 by adding vO to r0*(OH1, OH2) unit vectors
             for(int j=0;j<3;j++){
                 xH1[i][j]=xO[i][j] + p_r0*p1[j];
                 xH2[i][j]=xO[i][j] + p_r0*p2[j];
                 xS[i][j] =xO[i][j] + p_rSigma*s[j];
             }

            xH1[i][3]=xH2[i][3]=p_q;
            xS[i][3]=-4.0*p_q;
        }

        //subtitute to grandCoords
        grandCoords=new double[xO.length+xH1.length+xH2.length+xS.length][4];
        for(int j=0;j<4;j++){
            int c=0;
            for(int i=0;i<xO.length;i++)   grandCoords[c+i][j]=xO[i][j];
            c+=xO.length;
            for(int i=0;i<xH1.length;i++)   grandCoords[c+i][j]=xH1[i][j];
            c+=xH1.length;
            for(int i=0;i<xH2.length;i++)   grandCoords[c+i][j]=xH2[i][j];
            c+=xH2.length;
            for(int i=0;i<xS.length;i++)   grandCoords[c+i][j]=xS[i][j];

        }

        if(debug>=1){
            System.out.printf("%d\n\n",4*nO);
            for(int i=0;i<nO;i++){
                System.out.printf("O %12.6f %12.6f %12.6f  %6.3f\n", xO[i][0], xO[i][1], xO[i][2], xO[i][3]);
                System.out.printf("H %12.6f %12.6f %12.6f  %6.3f\n", xH1[i][0],xH1[i][1],xH1[i][2],xH1[i][3]);
                System.out.printf("H %12.6f %12.6f %12.6f  %6.3f\n", xH2[i][0],xH2[i][1],xH2[i][2],xH2[i][3]);
                System.out.printf("H %12.6f %12.6f %12.6f  %6.3f\n", xS[i][0], xS[i][1], xS[i][2], xS[i][3]);
            }
        }

    }

    double TwoBodyInteraction(boolean isGrad){
		double V2body=0;
		double r11, r21, r12, r22, theta, cost, dt, dt2,fcutoff,temp;
		for(int i=0;i<xS.length;i++){
            for(int j=i+1;j<xS.length;j++){
				double rij2= (Math.pow( (xS[i][0]-xS[j][0]),2)
							+ Math.pow( (xS[i][1]-xS[j][1]),2)
							+ Math.pow( (xS[i][2]-xS[j][2]),2))/Math.pow(p_sigma, 2);

				double rij6 = Math.pow(rij2,3);
                
				V2body+=4.0*p_epsilon*(1.0- rij6)/(Math.pow(rij6,2));
			}            
		}
		return V2body;
	}

    double ChargeChargeInteraction(boolean isGrad){
		double Vq=0;
		for(int i=0; i<grandCoords.length; i++ ){
			for(int j=i+1; j< grandCoords.length;j++){
                double rij= Math.sqrt(
                              Math.pow( (grandCoords[i][0]-grandCoords[j][0]),2)
							+ Math.pow( (grandCoords[i][1]-grandCoords[j][1]),2)
							+ Math.pow( (grandCoords[i][2]-grandCoords[j][2]),2));

				double temp1=grandCoords[i][3]*grandCoords[j][3] / rij *rescale;
				Vq+= temp1;
//				if(isGrad){ //for gradient
//					for(int k=0;k<3;k++){
//						double temp2=(-q[i]*q[j]/r2[i][j])*dr[i][j][k]*rescale;
//						grad[i][k]+= temp2;
//						grad[j][k]+=-temp2;
//					}
//				}
			}
		}
		return Vq;
	}

    void CalcElectricalField(){
		// <part 4: matrix E>
		double[] tmp=new double[3];

		for(int i=0; i< grandCoords.length;i++){
			for(int k=0;k<3;k++) tmp[k] = 0;

			for(int j=0;j< grandCoords.length; j++) if(i!=j)
				for(int k=0;k<3;k++) {
					double rijk = grandCoords[j][k] - grandCoords[i][k];
					tmp[k] += grandCoords[j][3] * rijk;
				}

			for(int k=0;k<3;k++){
				E[i][k] = tmp[k];
			}
		}
	}

//	void CalcDipoleMoment(){
//		double[][] Dtmp=new double[3*nO][3*nO];
//		double[] Etmp=new double[3*nO];
//
//
//		for(int i=0; i< 3*nO;i++){
//			for(int j=i;j< 3*nO; j++){
//				if(i==j)
//					Dtmp[i][j] = Dtmp[j][i]=  ( 1 + p_alpha*D[i][j]);
//				else
//					Dtmp[i][j] = Dtmp[j][i]=  ( p_alpha*D[i][j]);
//			}
//			Etmp[i]=E[i];
//		}
//    }
//
//	double ChargeDipoleInteraction(boolean isGrad){
//		double Vcd = 0;
//		for(int i=0;i< 3*nO; i++){
//			Vcd -= mu[i] * E[i] / p_alpha *rescale;
//		}
//
//		if(isGrad){
//				for(int i=0;i<nAtom;i++){
//					for(int j=0;j<nO;j++){
//						if(i!=j){
//							for(int t=0;t<3;t++){
//								double tmp = 0;
//								for(int k=0;k<3;k++) tmp += mu[3*j+k] * rv[i][j][k];
//
//								tmp =q[i]*tmp*dScd[i][j]*dr[i][j][t];
//								tmp+=q[i]*mu[3*j+t]*Scd[i][j];
//
//								tmp*=rescale;
//
//								grad[i][t]+=tmp;
//								grad[j][t]-=tmp;
//							}
//						}
//				}
//			}
//		}
//		return Vcd;
//	}
//
//	double DipoleDipoleInteraction(boolean isGrad){
//		double Vdd = 0;
//		for(int i=0;i< 3*nO; i++){
//			for(int j=i+1;j< 3*nO; j++){
//				Vdd += mu[i]*D[i][j]*mu[j] *rescale;
//			}
//		}
//
//		if(isGrad){
//			for(int i=0;i<nO;i++){
//				for(int j=0;j<nO;j++){
//					if(i!=j){
//						double pt1 = 0;
//						for(int k=0;k<3;k++)
//							for(int l=0;l<3;l++)
//								pt1 += mu[3*i+k]*T[3*i+k][3*j+l]*mu[3*j+l];
//
//						double pt2 = 0;
//						double pt3 = 0;
//
//						for(int k=0;k<3;k++) {
//							pt2 += mu[i*3+k] * dr[i][j][k];
//							pt3 += mu[j*3+k] * dr[i][j][k];
//						}
//
//						for(int t=0;t<3;t++){
//							double temp = pt1*dSdd[i][j]*dr[i][j][t];
//							temp += 3.0*(2.0*dr[i][j][t]*pt2*pt3
//							- mu[j*3+t]*pt2
//							- mu[i*3+t]*pt3)/r[i][j]*Sdd[i][j];
//
//							grad[i][t]+=temp*rescale;
//						}
//					}
//				}
//			}
//		}
//		return Vdd;
//	}
//
//	double SelfDipoleInteraction(){
//		double Vsd = 0;
//		for(int i=0;i< 3*nO; i++) Vsd += mu[i]*mu[i];
//		Vsd /= 2*p_alpha;
//		Vsd*=rescale;
//		return Vsd;
//	}

	//================== for gradient
}
