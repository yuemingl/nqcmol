/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.tools;

import java.lang.reflect.Array;
import java.util.Random;
import javax.vecmath.AxisAngle4d;
import javax.vecmath.Matrix3d;
import javax.vecmath.Point3d;
import javax.vecmath.Vector3d;

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
	public static void NORMALIZE(double[] a){
		double sum=Math.sqrt(DOTPRODUCT(a,a));
		for(int i=0;i<a.length;i++)
			if(sum!=0) a[i]/=sum;
	}

	//normalise for vector[3]
	public static void NORMALIZE(double[][] A){
		for (int i = 0; i< A.length; i++){
			NORMALIZE(A[i]);
		} 
	}

	//normalise for vector[3]
	public static void NORMALIZE(double[][] A,int nRows,int nCols){
		for (int i = 0; i<nRows; i++){
			double sum = 0;
			for (int j = 0; j<nCols; j++){
				sum += A[i][j]*A[i][j];
			}
			sum = Math.sqrt(sum);
			for (int j = 0; j<nCols; j++){
				A[i][j]/=sum;
			}
		}
	}

//	public static void VEC_MUL_VEC(a,B,C,alpha,s0) for(int i=0;i<(s0);i++) for(int j=0;j<(s0);j++) (a)[i][j]=(alpha)*(B)[i]*(C)[j];
//
//	public static void VEC_MUL_MAT2(a,B,C,alpha) for(i=0;i<DGR;i++) (a)[i]=(alpha)*DOTPRODUCT((B)[i],(C));
//
	public static void VEC_PLUS_VEC(double[] a,double[] b,double[] c,double alpha,double beta){
		for(int i=0;i<a.length;i++)
			a[i]=(alpha)*b[i]+(beta)*c[i];
	}

	public static void VEC_PLUS_VEC(int[] a,int[] b,int[] c,int alpha,int beta){
		for(int i=0;i<a.length;i++)
			a[i]=(alpha)*b[i]+(beta)*c[i];
	}

	public static void VEC_PLUS_NUM(double[] a,double[] b,double num){
		for(int i=0;i<a.length;i++) a[i]=b[i]+(num);
	}

	public static void VEC_MUL_NUM(double[] a,double[] b,double num){
		for(int i=0;i<a.length;i++) a[i]=b[i]*(num);
	}

	public static void VEC_EQU_NUM(int[] a,int num){
		if(a!=null)	for(int i=0;i<a.length;i++) a[i]=num;
	}
//
	public static void VEC_EQU_NUM(double[] a,double num){
		if(a!=null)
			for(int i=0;i<a.length;i++) a[i]=num;
	}

	public static String VEC_TO_STRING(int[] a){
		String s="";
		if(a!=null)
			for(int i=0;i<a.length;i++) s+=Integer.toString(a[i])+" ";
		return s;
	}
//
//
//	/*==================== for 2D matrix */
//	public static void MAT2_MUL_NUM(a,B,num,s0,s1) for(i=0;i<(s0);i++)	for(j=0;j<(s1);j++) (a)[i][j]=(B)[i][j]*(num);
//	public static void MAT2_EQU_MAT2(a,B,s0,s1) for(i=0;i<(s0);i++)	for(j=0;j<(s1);j++) (a)[i][j]=(B)[i][j];
//
	public static void MAT2_PLUS_NUM(double[][] A,double[][] B,double num){
		if(A!=null)
		for(int i=0;i<A.length;i++)
			for(int j=0;j<A[i].length;j++) A[i][j]=B[i][j]+num;
	}

	public static void MAT2_EQU_NUM(double[][] A,double num){
		if(A!=null)
		for(int i=0;i<A.length;i++)
			for(int j=0;j<A[i].length;j++) A[i][j]=(num);
	}
//	public static void MUL_NUM(a,num,s0,s1) for(i=0;i<(s0);i++)	for(j=0;j<(s1);j++) (a)[i][j]*=(num);
//	public static void PLUS_NUM(a,num,s0,s1) for(i=0;i<(s0);i++)	for(j=0;j<(s1);j++) (a)[i][j]+=(num);
	public static void MAT2_MUL_MAT2(double[][] A,double[][] B,double[][] C){
        int s0=B.length;
        int s1=C.length;
        int s2=C[0].length;

        for(int i=0;i<(s0);i++)
            for(int n=0;n<(s2);n++){
                A[i][n]=0;
                for(int j=0;j<(s1);j++) A[i][n]+=B[i][j]*C[j][n];
            }

    }
//
	public static void MAT2_MUL_VEC(double[] A,double B[][],double C[]){
		for(int i=0;i<A.length;i++){
			A[i]=0;
			for(int m=0;m<C.length;m++) A[i]+=B[i][m]*(C)[m];
		}
	}



//	public static void MAT2_PLUS_MAT2(a,B,C,s0,s1) for(i=0;i<(s0);i++)	for(j=0;j<(s1);j++) (a)[i][j]=(B)[i][j]+(C)[i][j];
//
//	public static void PLUS_MAT2(a,B,s0,s1) for(i=0;i<(s0);i++)	for(j=0;j<(s1);j++) (a)[i][j]+=(B)[i][j];
//
//	//for converting
	public static void MAT2_TO_VEC(double[][] A,double[] b){
        int s0=A.length;

        for(int i=0;i<(s0);i++){
            int s1=A[i].length;
            for(int j=0;j<(s1);j++)
                b[i*(s1)+j]=A[i][j];
        }
    }
//
	public static void VEC_TO_MAT2(double[] a, double[][] B){
        for(int i=0;i<B.length;i++){
            int s1=B[i].length;
            for(int j=0;j<s1;j++) B[i][j]=a[i*s1+j];
        }
    }


	public static double CalculateRMS(double[] x){
		double rms=0;
		for(int i=0;i<x.length;i++){
			rms+=x[i]*x[i];
		}
		rms=Math.sqrt(rms/x.length);
		return rms;
	}

	// for eigen vector solver


    /**
     * Solve A*x=b, where A is a square matrix, x and b are vectors.
     * @param x the solution. x is needed to be generated first.
     */
	public static void LinearSolver_apache(double[][] A, double[] b, double[] x) {
//		if(A.length==1){
//			x[0]=b[0]/A[0][0];
//		}else{
		org.apache.commons.math.linear.RealMatrix At =new org.apache.commons.math.linear.Array2DRowRealMatrix(A,false);
		org.apache.commons.math.linear.DecompositionSolver solver = new org.apache.commons.math.linear.LUDecompositionImpl(At).getSolver();
        org.apache.commons.math.linear.RealVector bt = new org.apache.commons.math.linear.ArrayRealVector(b,false);
		org.apache.commons.math.linear.RealVector solution = solver.solve(bt);
		double[] sol=solution.getData();
		System.arraycopy(sol, 0, x, 0, sol.length);
//		}
	}

    /**
     * Inverse matrix A
     * @param X the solution. X is needed to be generated first.
     */
    public static void InverseMatrix_apache(double[][] A, double[][] X) {
//		if(A.length==1){
//			x[0]=b[0]/A[0][0];
//		}else{
		org.apache.commons.math.linear.RealMatrix At =new org.apache.commons.math.linear.Array2DRowRealMatrix(A,false);
		org.apache.commons.math.linear.DecompositionSolver solver = new org.apache.commons.math.linear.LUDecompositionImpl(At).getSolver();
        org.apache.commons.math.linear.RealMatrix solution = solver.getInverse();

		double[][] sol=solution.getData();
        MAT2_PLUS_NUM(X,sol,0);
//		}
	}

	/**
	 * Sort the eigen values and corresponding eigen vectors wiht Bubble sort algorithm
	 * @param IsAscending: ascending or not
	 */
	public static void SortEigenvaluesAndEigenVectors(double[] eigenValues, double[][] eigenVectors,boolean IsAscending){
		for (int i=0;i<eigenValues.length;i++)
			for (int j=eigenValues.length-1; j > i; j--)
				if(   (( IsAscending) && ( eigenValues[j-1] > eigenValues[j]))
				   || ((!IsAscending) && ( eigenValues[j-1] < eigenValues[j])) ) {
					double swap=eigenValues[j];
					eigenValues[j]=eigenValues[j-1];
					eigenValues[j-1]=swap;

					for(int k=0;k<eigenVectors[j].length;k++){
						swap=eigenVectors[j][k];
						eigenVectors[j][k]=eigenVectors[j-1][k];
						eigenVectors[j-1][k]=swap;
					}
				}
	}

	public static void Eigen_apache(double[][] A,double[][] eigenVec,double[] eigenValue){
		org.apache.commons.math.linear.RealMatrix X=new org.apache.commons.math.linear.Array2DRowRealMatrix(A);
		org.apache.commons.math.linear.EigenDecomposition eigen=new org.apache.commons.math.linear.EigenDecompositionImpl(X,org.apache.commons.math.util.MathUtils.SAFE_MIN);
		double[] value=eigen.getRealEigenvalues();
//		System.err.printf(" Eigen_apache values \n");
//		PrintArray(value);
		System.arraycopy(value,0,eigenValue,0,value.length);
		for(int i=0;i<value.length;i++){
			org.apache.commons.math.linear.RealVector vec=eigen.getEigenvector(i);
			double[] tmp=vec.getData();
			System.arraycopy(tmp,0,eigenVec[i],0,tmp.length);
		}
	}

	public static void Eigen_colt(double[][] A,double[][] eigenVec,double[] eigenValue){
		cern.colt.matrix.DoubleMatrix2D X=new cern.colt.matrix.impl.DenseDoubleMatrix2D(A);

		cern.colt.matrix.linalg.EigenvalueDecomposition eigen=new  cern.colt.matrix.linalg.EigenvalueDecomposition(X);
		double[] value=eigen.getRealEigenvalues().toArray();
//		System.err.printf(" Eigen_apache values \n");
//		PrintArray(value);
		System.arraycopy(value,0,eigenValue,0,value.length);

		double[][] vec=eigen.getV().toArray();
		for(int i=0;i<vec.length;i++)
			for(int j=0;j<vec[i].length;j++){
				//eigenVec[i][j]=vec[i][j];
				System.arraycopy(vec[i],0,eigenVec[i],0,vec[i].length);
			}
	}

	public static void Eigen_jama(double[][] A,double[][] eigenVec,double[] eigenValue){
		/*
		Matrix X=new cern.colt.matrix.impl.DenseDoubleMatrix2D(A);

		Jama.EigenvalueDecomposition eigen=new Jama.EigenvalueDecomposition(X);
		double[] value=eigen.getRealEigenvalues().toArray();
//		System.err.printf(" Eigen_apache values \n");
//		PrintArray(value);
		System.arraycopy(value,0,eigenValue,0,value.length);

		double[][] vec=eigen.getV().toArray();
		for(int i=0;i<vec.length;i++)
			for(int j=0;j<vec[i].length;j++){
				eigenVec[i][j]=vec[i][j];
				//System.arraycopy(vec[i],0,eigenVec[i],0,vec[i].length);
			}
		 * */
	}



	public static void PrintArray(int[][] A){
		for(int i=0;i<A.length;i++){
			for(int j=0;j<A[i].length;j++)
				System.err.printf("%d ",A[i][j]);
			System.err.printf("\n");
		}
	}

	public static void PrintArray(double[][] A){
		for(int i=0;i<A.length;i++){
			for(int j=0;j<A[i].length;j++)
				System.err.printf("%1.9f ",A[i][j]);
			System.err.printf("\n");
		}
	}

	public static void PrintArray(int[] a){
			for(int j=0;j<a.length;j++)
				System.err.printf("%d ",a[j]);
			System.err.printf("\n");
	}

	public static void PrintArray(double[] a){
			for(int j=0;j<a.length;j++)
				System.err.printf("%1.9f ",a[j]);
			System.err.printf("\n");
	}


	public static String getExtension(String fullPath) {
			int dot = fullPath.lastIndexOf(".");
			return fullPath.substring(dot + 1);
	  }

	public static String getFilename(String fullPath) { // gets filename without extension
		int dot = fullPath.lastIndexOf(".");
		int sep = fullPath.lastIndexOf("/");
		return fullPath.substring(sep + 1, dot);
	  }

	public static String getPath(String fullPath) {
		int sep = fullPath.lastIndexOf("/");
		return fullPath.substring(0, sep);
	  }

      /**
       * @param src source array
       * @param pos index of first element to be erased
       * @param num number of elements to be erased
       */
    public static void eraseElementsFromArray(Object src,Object dest,int pos,int num) {
          //result = (T[])new Object[src.length - num];
          //System.out.printf("d =%d",Array.getLength(src));
          System.arraycopy(src, 0, dest, 0, pos);
          System.arraycopy(src, pos + num, dest, pos, Array.getLength(src) - pos - num );
      }

    public static void generateRandomVector(double vec[],Random ran,double scale){
        for(int i=0;i<vec.length;i++){
            vec[i]=scale*ran.nextDouble();
        }
    }

    private static Vector3d rotate(Vector3d vector, Vector3d axis, double angle) {
        Matrix3d rotate = new Matrix3d();
        rotate.set(new AxisAngle4d(axis, angle));
        Vector3d result = new Vector3d();
        rotate.transform(vector, result);
        return result;
    }

    public static double[] zmatrixToCartesian(double distances, double[] first_atoms,
                                        double angles,   double[] second_atoms,
                                        double dihedrals, double[] third_atoms) {      

       Vector3d f1 = new Vector3d(first_atoms);
       Vector3d f2 = new Vector3d(second_atoms);
       Vector3d f3 = new Vector3d(third_atoms);

       Vector3d cd = new Vector3d();       cd.sub(f3, f2);

       Vector3d bc = new Vector3d();      bc.sub(f2, f1);

       Vector3d n1 = new Vector3d();
       n1.cross(cd, bc);

        Vector3d n2 = rotate(n1,bc,-dihedrals);
        Vector3d ba = rotate(bc,n2,-angles);

        ba.normalize();
        ba.scale(distances);

        Point3d result = new Point3d();
        result.add(f1, ba);
        double[] cartesianCoords = new double[3];
        cartesianCoords[0] = result.x;
        cartesianCoords[1] = result.y;
        cartesianCoords[2] = result.z;

        return cartesianCoords;
    }
}
