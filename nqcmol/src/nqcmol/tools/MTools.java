/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.tools;

import org.apache.commons.math.linear.Array2DRowRealMatrix;
import org.apache.commons.math.linear.EigenDecomposition;
import org.apache.commons.math.linear.EigenDecompositionImpl;
import org.apache.commons.math.linear.RealMatrix;
import org.apache.commons.math.linear.RealVector;
import org.apache.commons.math.util.MathUtils;

/**
 *
 * @author nqc
 */
public  class  MTools{

	public static double DOTPRODUCT(double[] a,double[] b){
		double result=0;
		for(int i=0;i<a.length;i++) result+=a[i]*b[i];
		return result;
	}

	//vector product, only for vector [3]
	public static void CROSSPRODUCT(double[] c,double[] a,double[] b){
		c[0]= a[1]*b[2] - a[2]*b[1];
		c[1]=-a[0]*b[2] + a[2]*b[0];
		c[2]= a[0]*b[1] - a[1]*b[0];
	}

	public static double SQR(double x){
		return x*x;
	}

	//normalise for vector[3]
	public static void NORMALISE(double[] a){
		double sum_0000=Math.sqrt(DOTPRODUCT(a,a));
		for(int i=0;i<a.length;i++)
			if(sum_0000!=0) a[i]/=sum_0000;
	}

	//normalise for vector[3]
	public static void NORMALISE(double[][] A){
		for (int i = 0; i< A.length; i++){
			NORMALISE(A[i]);
		} 
	}

//	public static void VEC_MUL_VEC(a,B,C,alpha,s0) for(int l=0;l<(s0);l++) for(int m=0;m<(s0);m++) (a)[l][m]=(alpha)*(B)[l]*(C)[m];
//
//	public static void VEC_MUL_MAT2(a,B,C,alpha) for(l=0;l<DGR;l++) (a)[l]=(alpha)*DOTPRODUCT((B)[l],(C));
//
//	public static void VEC_PLUS_VEC(a,B,C,s0,alpha,beta) for(l=0;l<(s0);l++){\
//	 (a)[l]=(alpha)*(B)[l]+(beta)*(C)[l];}
//
//	public static void VEC_PLUS_NUM(a,B,s0,num) for(l=0;l<(s0);l++) (a)[l]=(B)[l]+(num);
//	public static void VEC_MUL_NUM(a,B,s0,num) for(l=0;l<(s0);l++) (a)[l]=(B)[l]*(num);
//	public static void VEC_EQU_NUM(a,s0,num) for(l=0;l<(s0);l++) (a)[l]=(num);
//
//	/*==================== for 2D matrix */
//	public static void MAT2_MUL_NUM(a,B,num,s0,s1) for(l=0;l<(s0);l++)	for(m=0;m<(s1);m++) (a)[l][m]=(B)[l][m]*(num);
//	public static void MAT2_EQU_MAT2(a,B,s0,s1) for(l=0;l<(s0);l++)	for(m=0;m<(s1);m++) (a)[l][m]=(B)[l][m];
//
//	public static void MAT2_PLUS_NUM(a,B,num,s0,s1) for(l=0;l<(s0);l++)	for(m=0;m<(s1);m++) (a)[l][m]=(B)[l][m]+(num);
	public static void MAT2_EQU_NUM(double[][] A,double num){
		for(int l=0;l<A.length;l++)
			for(int m=0;m<A[l].length;m++) A[l][m]=(num);
	}
//	public static void MUL_NUM(a,num,s0,s1) for(l=0;l<(s0);l++)	for(m=0;m<(s1);m++) (a)[l][m]*=(num);
//	public static void PLUS_NUM(a,num,s0,s1) for(l=0;l<(s0);l++)	for(m=0;m<(s1);m++) (a)[l][m]+=(num);
//	public static void MAT2_MUL_MAT2(a,B,C,s0,s1,s2) for(l=0;l<(s0);l++) for(n=0;n<(s2);n++){ (a)[l][n]=0;for(m=0;m<(s1);m++) (a)[l][n]+=(B)[l][m]*(C)[m][n];}
//
	public static void MAT2_MUL_VEC(double[] A,double B[][],double C[]){
		for(int l=0;l<A.length;l++){
			A[l]=0;
			for(int m=0;m<C.length;m++) A[l]+=B[l][m]*(C)[m];
		}
	}

//	public static void MAT2_PLUS_MAT2(a,B,C,s0,s1) for(l=0;l<(s0);l++)	for(m=0;m<(s1);m++) (a)[l][m]=(B)[l][m]+(C)[l][m];
//
//	public static void PLUS_MAT2(a,B,s0,s1) for(l=0;l<(s0);l++)	for(m=0;m<(s1);m++) (a)[l][m]+=(B)[l][m];
//
//	//for converting
//	public static void MAT2_TO_VEC(a,B,s0,s1) for(l=0;l<(s0);l++) for(m=0;m<(s1);m++) (B)[l*(s1)+m]=(a)[l][m];
//
//	public static void VEC_TO_MAT2(B,a,s0,s1) for(l=0;l<(s0);l++) for(m=0;m<(s1);m++) (a)[l][m]=(B)[l*(s1)+m];

	public static void SortVector(double[] A){
		
	}

	public static void Eigen(double[][] A,double[][] eigenVec,double[] eigenValue){
		RealMatrix X=new Array2DRowRealMatrix(A);
		EigenDecomposition eigen=new EigenDecompositionImpl(X,MathUtils.SAFE_MIN);
		double[] value=eigen.getRealEigenvalues();
//		System.err.printf(" Eigen values \n");
//		PrintArray(value);
		System.arraycopy(value,0,eigenValue,0,value.length);
		
		for(int i=0;i<value.length;i++){
			RealVector vec=eigen.getEigenvector(i);
			double[] tmp=vec.getData();
			System.arraycopy(tmp,0,eigenVec[i],0,tmp.length);
		}
	}

	public static void PrintArray(double[][] A){
		for(int i=0;i<A.length;i++){
			for(int j=0;j<A[i].length;j++)
				System.err.printf("%1.9f ",A[i][j]);
			System.err.printf("\n");
		}
	}

	public static void PrintArray(double[] a){
			for(int j=0;j<a.length;j++)
				System.err.printf("%1.9f ",a[j]);
			System.err.printf("\n");
	}
}
