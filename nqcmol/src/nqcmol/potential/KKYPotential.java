/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.potential;

import java.io.File;
import java.io.FileNotFoundException;
import java.util.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.Cluster;

/**
 *
 * @author nqc
 */
public class KKYPotential extends Potential{
	public KKYPotential(){
		nativeUnit="kJ/mol";
		unit=nativeUnit;

        //params.put(String("O_z"),Double(-0.96));
        params.clear();
        params.put("O_z",-0.92);
        params.put("O_a",1.728);
        params.put("O_b",0.1275);
        params.put("O_c",56.06);

        params.put("H_z",0.46);
        params.put("H_a",0.035);
        params.put("H_b",0.044);
        params.put("H_c",0.0);

        params.put("O-H_D1",57394.93);
        params.put("O-H_beta1",7.400);
        params.put("O-H_D2",-2189.30);
        params.put("O-H_beta2",3.130);
        params.put("O-H_D3",34.74);
        params.put("O-H_beta3",12.800);
        params.put("O-H_r3",1.283);

        params.put("H-O-H_f",69.253);
        params.put("H-O-H_theta",99.5);
        params.put("H-O-H_rm",1.43);
        params.put("H-O-H_gr",9.20);  
	};

	@Override
	public String getEquation(){
		String equation="KKY";
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
	public void setParam(String filename) {
		if(!filename.isEmpty()){
		try {
            params.clear();
			Scanner scanner = new Scanner(new File(filename));
            while(scanner.hasNext()){
                String key=scanner.next();
                String sValue=scanner.next();
                Double val=Double.valueOf(sValue);
                params.put(key, val);
            }

//            Set set = params.entrySet();
//            Iterator i = set.iterator();
//            System.out.println("Parameters:" );
//            while(i.hasNext()){
//              Map.Entry me = (Map.Entry)i.next();
//              System.out.println(me.getKey() + " = " + me.getValue() );
//            }
            
		} catch (FileNotFoundException ex) {
			Logger.getLogger(KKYPotential.class.getName()).log(Level.SEVERE, null, ex);
		}
		}
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

		double V_2body=TwoBodyInteraction(false);
		double V_3body=ThreeBodyInteraction(false);

		energy=V_2body + V_3body;

//	   System.out.println("Summary ");
//       System.out.println("     1. v_Coulomb         ="+v_Coulomb);
//       System.out.println("     2. v_ShortRangeRepuls="+v_ShortRangeRepuls);
//       System.out.println("     3. v_VanDerWaals     ="+v_VanDerWaals);
//       System.out.println("     4. v_CovalentBond    ="+v_CovalentBond);
//	   System.out.println(" 5. V_2body (1+2+3+4) = "+V_2body);
//       System.out.println(" 6. V_3body           = "+V_3body);
//	   System.out.println(" 7. Energy (5+6)      = "+energy);

		energy=ConvertUnit(energy,nativeUnit,unit);

		nEvals++;
		return energy;
	}

	@Override
	protected void Gradient_(double[] _p, double[] gradient) {
		//System.out.print(" Size of vec = "+ Integer.toString(_p.length));
		if(_p.length<=3) return ;
		
		AllocatePrivateVariables(_p.length/3);

		CalcDistanceAndCharge(_p,true);

		double V_2body=TwoBodyInteraction(true);
		double V_3body=ThreeBodyInteraction(true);


		for(int i=0;i<nAtom;i++)
			for(int k=0;k<3;k++)
				gradient[i*3+k]=ConvertUnit(this.grad[i][k],nativeUnit,unit);


		nEvals+=3;
	}

	

	//========================== implement OSS2
	int nAtom,nO;
	double[][] r,x,r2;


	//for gradient
	double[][][] rv,dr;
	double[][] grad;

    //parameters
    HashMap<String,Double> params = new HashMap<String,Double>();   
        
	

	double SQR(double x){
		return x*x;
	}

	private void AllocatePrivateVariables(int nAtoms_){
		if(nAtom!=nAtoms_){
			nAtom=nAtoms_;
			x=new double[nAtom][3];
			r=new double[nAtom][nAtom];
			r2=new double[nAtom][nAtom];

			//if(true){ ///for gradient
				grad=new double[nAtom][3];
				rv=new double[nAtom][nAtom][3];
				dr=new double[nAtom][nAtom][3];
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
        
        // distance vectors, euclidian distances and their derivations
        for (int i = 0; i < nAtom; i += 1){
            for (int j = i+1; j < nAtom; j += 1){
                r2[i][j] = 0;

                for (int dim = 0; dim < 3; dim += 1){
                    rv[i][j][dim] = x[i][dim] - x[j][dim];
                    rv[j][i][dim] = -rv[i][j][dim];
                    r2[i][j] += SQR(rv[i][j][dim]);
                }
                r2[j][i]=r2[i][j];
                r[i][j] = r[j][i] = Math.sqrt(r2[i][j]);
                              

                if(isGrad){
                    for (int dim = 0; dim < 3; dim += 1){
                        dr[i][j][dim] = rv[i][j][dim]/r[i][j];
                        dr[j][i][dim] =-dr[i][j][dim];
                    }
                }
            }
        }
		
	}

    double getParam1(int i,String key){
        double answer=0;        
        String query=cluster.getAtomicSymbol(i)+"_"+key;
        Double val= params.get(query);
        if(val!=null) answer=val.doubleValue();

        //System.out.println("Key="+query+ " Value="+answer);
        return answer;
    }

    double getParam2(int i,int j,String key){
        double answer=0;
        String query=cluster.getAtomicSymbol(i)+"-"+cluster.getAtomicSymbol(j)+"_"+key;
        Double val= params.get(query);
        if(val!=null) answer=val.doubleValue();
        else{
            query=cluster.getAtomicSymbol(j)+"-"+cluster.getAtomicSymbol(i)+"_"+key;
            val= params.get(query);
            if(val!=null) answer=val.doubleValue();
        }        
        //System.out.println("Key="+query+ " Value="+answer);
        return answer;
    }

    double getParam2(String symbol,String key){
        double answer=0;
        String query=symbol+"_"+key;
        Double val= params.get(query);
        if(val!=null) answer=val.doubleValue();
        else{
            query=symbol+"_"+key;
            val= params.get(query);
            if(val!=null) answer=val.doubleValue();
        }
        //System.out.println("Key="+query+ " Value="+answer);
        return answer;
    }

     double getParam3(int i,int j,int k,String key){
        double answer=0;
        String query=cluster.getAtomicSymbol(i)+"-"+cluster.getAtomicSymbol(j)+"-"+cluster.getAtomicSymbol(k)+"_"+key;
        Double val= params.get(query);
        if(val!=null){
            answer=val.doubleValue();           
        }else answer=Double.NaN;
        //System.out.println("Key="+query+ " Value="+answer);

        return answer;
    }

    double v_Coulomb=0;
    double v_ShortRangeRepuls=0;
    double v_VanDerWaals=0;
    double v_CovalentBond=0;

	double TwoBodyInteraction(boolean isGrad){
        double V_pair=0;
        v_Coulomb=0;
        v_ShortRangeRepuls=0;
        v_VanDerWaals=0;
        v_CovalentBond=0;
        
        final double p_f0=41.865*0.1; //constant parameters
        final double kCoulomb=1389.35444;//= 8.9875517873681764e9*(1.60217646e-19)^2/1e-10 * 0.001*6.0221415e23
        //final double kCoulomb=1389.354678;//=from OSS2 potential

        for (int i = 0; i < nAtom; i += 1){
            double p_zi,p_ai,p_bi,p_ci;
            p_zi=getParam1(i,"z");
            p_ai=getParam1(i,"a");
            p_bi=getParam1(i,"b");
            p_ci=getParam1(i,"c");


            for (int j = i+1; j < nAtom; j += 1){
                double p_zj,p_aj,p_bj,p_cj;
                p_zj=getParam1(j,"z");
                p_aj=getParam1(j,"a");
                p_bj=getParam1(j,"b");
                p_cj=getParam1(j,"c");

                double p_D1,p_D2,p_D3;
                p_D1=getParam2(i,j,"D1");
                p_D2=getParam2(i,j,"D2");
                p_D3=getParam2(i,j,"D3");

                double p_beta1,p_beta2,p_beta3;
                p_beta1=getParam2(i,j,"beta1");
                p_beta2=getParam2(i,j,"beta2");
                p_beta3=getParam2(i,j,"beta3");

                double p_r3=getParam2(i,j,"r3");


//                double expij   = Math.exp(( p_ai + p_aj - r[i][j])/(p_bi + p_bj));
//                double expij_1 = Math.exp( -p_beta1*r[i][j]);
//                double expij_2 = Math.exp( -p_beta2*r[i][j]);
//                double expij_3 = Math.exp( -p_beta3* SQR(r[i][j] - p_r3));

//                v_Coulomb+=kCoulomb*p_zi*p_zj / r[i][j];
//                v_ShortRangeRepuls+= p_f0 * (p_bi + p_bj) * expij;
//                v_VanDerWaals+= - p_ci*p_cj/Math.pow(r[i][j],6);
//                v_CovalentBond+= p_D1 * expij_1 + p_D2 * expij_2 + p_D3 * expij_3;
//
//                if(isGrad){
//                    double gd= (-kCoulomb*p_zi*p_zj / r2[i][j])
//                         +(-p_f0 * expij)
//                         +(6.0 * p_ci*p_cj /Math.pow(r[i][j],7))
//                         +(-(p_D1*p_beta1)*expij_1)
//                         +(-(p_D2*p_beta2)*expij_2)
//                         +(2.0*p_D3*p_beta3*(p_r3 - r[i][j]))*expij_3;
//
//                    for (int dim = 0; dim < 3; dim += 1){
//                        grad[i][dim] +=  gd*dr[i][j][dim];
//                        grad[j][dim] += -gd*dr[i][j][dim];
//                    }
//                }
            }
        }
      

       V_pair+=v_Coulomb +  v_ShortRangeRepuls +  v_VanDerWaals + v_CovalentBond;

        return V_pair;
	}	

	double ThreeBodyInteraction(boolean isGrad){
		double V_3b = 0;

        for(int i = 0; i < nAtom; i += 1){
            for(int j = i+1; j < nAtom; j += 1){
                for(int k = j+1; k < nAtom; k += 1){

                    //boolean isExist1 = false, isExist2=false, isExist3=false, isExist4=false;
                    
                    double p_f=getParam3(j,i,k,"f");
                    double p_theta=getParam3(j,i,k,"theta");
                    double p_rm=getParam3(j,i,k,"rm");
                    double p_gr=getParam3(j,i,k,"gr");                                     

                    if(!Double.isNaN(p_f) && !Double.isNaN(p_theta) && !Double.isNaN(p_rm) && !Double.isNaN(p_gr)){ //if parameters exist
                        double cost = (r2[i][j] + r2[i][k] - r2[j][k]) / (2.0*r[i][j]*r[i][k]);
                        double dth = Math.acos(cost) - Math.toRadians(p_theta);

                        double expj = Math.exp(p_gr * (r[i][j] - p_rm) );
                        double expk = Math.exp(p_gr * (r[i][k] - p_rm) );

                        double V_f = -p_f * (Math.cos( 2.0*dth) - 1.0 );
                        double V_l = Math.sqrt( ( 1.0 / (expj + 1.0) ) * ( 1.0 / (expk + 1.0 ) ) );

                        V_3b += V_f * V_l;

                        //System.out.println(" V_3b = "+V_3b);

                        if(isGrad){
                            double sindt = Math.sin(2.0*dth);

                            double dacosthHOH = -2.0 / Math.sqrt(1.0 - cost*cost);

                            double dV_lj = V_l * p_gr / (2.0/expj + 2.0);
                            double dV_lk = V_l * p_gr / (2.0/expk + 2.0);

                            for (int dim = 0; dim < 3; dim += 1){
                                 double dcosHOHj = -dr[i][k][dim]/r[i][j] + dr[i][j][dim]*cost/r[i][j];
                                 double dcosHOHk = -dr[i][j][dim]/r[i][k] + dr[i][k][dim]*cost/r[i][k];

                                 double dV_fj = p_f * sindt * dacosthHOH * dcosHOHj;
                                 double dV_fk = p_f * sindt * dacosthHOH * dcosHOHk;

                                 double gdHj = dV_fj * V_l + V_f * dV_lj * dr[i][j][dim];
                                 double gdHk = dV_fk * V_l + V_f * dV_lk * dr[i][k][dim];

                                 grad[i][dim] += -gdHj-gdHk;
                                 grad[j][dim] += gdHj;
                                 grad[k][dim] += gdHk;
                            }
                        }
                    }

                }
            }
        }
        return V_3b;
	}

	
	//================== for gradient
}
