/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.tools;


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
//	public static void MAT2_PLUS_NUM(a,B,num,s0,s1) for(i=0;i<(s0);i++)	for(j=0;j<(s1);j++) (a)[i][j]=(B)[i][j]+(num);
	public static void MAT2_EQU_NUM(double[][] A,double num){
		if(A!=null)
		for(int i=0;i<A.length;i++)
			for(int j=0;j<A[i].length;j++) A[i][j]=(num);
	}
//	public static void MUL_NUM(a,num,s0,s1) for(i=0;i<(s0);i++)	for(j=0;j<(s1);j++) (a)[i][j]*=(num);
//	public static void PLUS_NUM(a,num,s0,s1) for(i=0;i<(s0);i++)	for(j=0;j<(s1);j++) (a)[i][j]+=(num);
//	public static void MAT2_MUL_MAT2(a,B,C,s0,s1,s2) for(i=0;i<(s0);i++) for(n=0;n<(s2);n++){ (a)[i][n]=0;for(j=0;j<(s1);j++) (a)[i][n]+=(B)[i][j]*(C)[j][n];}
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
//	public static void MAT2_TO_VEC(a,B,s0,s1) for(i=0;i<(s0);i++) for(j=0;j<(s1);j++) (B)[i*(s1)+j]=(a)[i][j];
//
//	public static void VEC_TO_MAT2(B,a,s0,s1) for(i=0;i<(s0);i++) for(j=0;j<(s1);j++) (a)[i][j]=(B)[i*(s1)+j];



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
}
