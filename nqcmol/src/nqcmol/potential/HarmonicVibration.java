/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.potential;

import nqcmol.Cluster;
import nqcmol.tools.MTools;

/**
 *
 * @author nqc
 */
public class HarmonicVibration {
	
	Cluster m;
	private int n_AtomOR3;
	private double[][] Hessian;
	private double[][] TR_D;
	private double[][] TR_Hessian_vect;
	private double[] TR_Hess_eig;
	
	private double[][] Mom_I=new double[3][3]; //!< inertia tensor
	private	double[] I_eig=new double[3]; //!< the eigenvalues of the inertia tensor
	private	double[][] I_vect=new double[3][3]; //!< the eigenvectors of the inertia tensor

	int n_TrD; //!< dimension of TR_D, usually 6 -- 5 in case of linear molecule
	private int nFreqs;

	/**
	 * Constructor
	 * @param m_ : reference of a cluster. Frequencies and other vibration info will be updated further.
	 */
	public HarmonicVibration(Cluster m_){
		setCluster(m_);
	}

	/**
	 * Constructor: require to call setCluster in the future if you want to perform calculation
	 */
	public HarmonicVibration(){
	}

	int setCluster(Cluster m_){
		m=m_;
		Hessian=m.getHessian().clone();
		if(n_AtomOR3!=m.getNcoords()){
			n_AtomOR3=m.getNcoords();
			TR_D=new double[n_AtomOR3][n_AtomOR3];

	//		mextra->Freqs.clear();
	//		mextra->RedMass.clear();
	//		mextra->kConsts.clear();
	//		mextra->lCART.clear();
		}
		return 0;
	}

	public int CalcFreq(){
		m.Center();
//			//#ifdef VIBANA_DEBUG
//			System.err.printf("Hessian \n");
//			MTools.PrintArray(Hessian);
//			//#endif
//			System.err.printf("Coordinates \n");
//			MTools.PrintArray(m.getCoords());

		m.getInertiaTensor(Mom_I);

//			System.err.printf("Inertia Tensor \n");
//			MTools.PrintArray(Mom_I);

		MTools.Eigen_apache(Mom_I,I_vect,I_eig); MTools.SortEigenvaluesAndEigenVectors(I_eig,I_vect,true);
		
		for(int j=0;j<3;j++){
			I_vect[1][j]=-I_vect[1][j];
			I_vect[2][j]=-I_vect[2][j];
		}
		//MTools.Eigen_colt(Mom_I,I_vect,I_eig);

//			System.err.printf("Inertial Tensor Eigen Values : ");
//			MTools.PrintArray(I_eig);
//			System.err.printf("Inertial Tensor Eigen Vectors : \n");
//			MTools.PrintArray(I_vect);

		MassweightedHessian();
//			//#ifdef VIBANA_DEBUG
//			System.err.printf("Mass weighted Hessian \n");
//			MTools.PrintArray(Hessian);
//			//#endif

		Calc_tr_D(); //calculate transform matrix D
//			System.err.printf("\n TrD \n");
//			MTools.PrintArray(TR_D);

		//#ifdef VIBANA_DEBUG
//			System.err.printf("Normalize..\n");
		//#endif
		//MTools.NORMALIZE(TR_D,n_TrD,n_AtomOR3);//normalize the lines of D
		MTools.NORMALIZE(TR_D);//normalize the lines of D
		//#ifdef VIBANA_DEBUG
//			System.err.printf("\n n_TrD=%d \n",n_TrD);
//			MTools.PrintArray(TR_D);
		//#endif
				
		//orthogonalization
		if(n_TrD<5) return -1;
		else Calc_Sch_Ort();

			//#ifdef VIBANA_DEBUG
//			System.err.printf("\n After n_TrD=%d \n",n_TrD);
//			System.err.printf("Tr_D \n");
//			MTools.PrintArray(TR_D);

//			System.err.printf("Apply the transformation to the Hessian	\n");
			//#endif

		Transform_Hessian();

			//#ifdef VIBANA_DEBUG
//			System.err.printf("\n TransformedHessian \n");
//			MTools.PrintArray(Hessian);
			//#endif
			
		//the submatrix
		
		nFreqs=n_AtomOR3-n_TrD;
		double[][] subHessian=new double[nFreqs][nFreqs];
		for (int i = n_TrD; i<n_AtomOR3; i++)
			for (int j = n_TrD; j<n_AtomOR3; j++)
				subHessian[i-n_TrD][j-n_TrD] = Hessian[i][j];

		//#ifdef VIBANA_DEBUG
//		System.err.printf("\n Diagonalize..\n");
		//#endif


		TR_Hessian_vect=new double[nFreqs][nFreqs];
		TR_Hess_eig=new double[nFreqs];
		//MTools.Eigen_apache(subHessian,TR_Hessian_vect, TR_Hess_eig);
		MTools.Eigen_colt(subHessian,TR_Hessian_vect, TR_Hess_eig);


		MTools.SortEigenvaluesAndEigenVectors(TR_Hess_eig,TR_Hessian_vect,true);// need to be changed in the future

		//#ifdef VIBANA_DEBUG
//		System.err.printf("\n\nTransformedHessian \n");
//		MTools.PrintArray(TR_Hessian_vect);
//		System.err.printf("\n\nEigen values \n");
//		MTools.PrintArray(TR_Hess_eig);
		//#endif

		//Sort_Freq();

		double[] freqs=new double[nFreqs];
		for(int i=0;i<nFreqs;i++){
			int sign=((TR_Hess_eig[i]>0)?1:-1);
			//freq[i]=sign*Math.sqrt(Math.abs(TR_Hess_eig[i]))*Math.sqrt(2.6255e+29)/(2.0*M_PI*29979245800);
			freqs[i]=sign*Math.sqrt(Math.abs(TR_Hess_eig[i]))*2.7202e+03;
		}
		m.setFreqs(freqs);

//		System.err.printf("\n\nFreq\n");
//		MTools.PrintArray(freqs);

		return 0;
	}
//
//void CalcReducedMass(){
//	mextra->lCART.clear();
//	//#ifdef VIBANA_DEBUG
//	System.out.printf("\n\nCalculate lCART \n";
//	//#endif
//	for(int k=0;k<nFreqs;k++){
//		vector<double> v;	v.clear();
//		for(int i=0;i<n_AtomOR3;i++){
//			double tmp=0;
//
//			for(int j=0;j<nFreqs;j++)	tmp+=TR_D[j+6][i]*TR_Hessian_vect[j][k];
//
//			OBAtom *a=m->GetAtom((int)(i/3)+1);
//			double m=Math.sqrt(a->GetExactMass());
//			tmp/=m;
//			v.push_back(tmp);
//		}
//		mextra->lCART.push_back(v);
//	}
//
//	mextra->RedMass.clear();
//	//#ifdef VIBANA_DEBUG
//	System.out.printf("\n\nCalculate RedMass \n";
//	//#endif
//	for(int i=0;i<nFreqs;i++){
//		double tmp=0;
//		for(int j=0;j<n_AtomOR3;j++){
//			tmp+=pow(mextra->lCART[i][j],2);
//		}
//		tmp=1.0/tmp;
//		mextra->RedMass.push_back(tmp);
//	}
//	//#ifdef VIBANA_DEBUG
//	System.out.printf("\n\nUpdate lCART \n";
//	//#endif
//
//	for(int i=0;i<nFreqs;i++){
//		for(int j=0;j<n_AtomOR3;j++)	mextra->lCART[i][j]*=Math.sqrt(mextra->RedMass[i]);
//	}
//}
//
//void CalcForceConstant(){
//	//#ifdef VIBANA_DEBUG
//	System.out.printf("\n\nCalculate kConsts \n";
//	//#endif
//	//double tmp=SQR(2.0*M_PI*29.978e9)*1.660538782e-27*1.0e-2;
//	double tmp=5.8913e-7;
//	mextra->kConsts.clear();
//	for(int i=0;i<nFreqs;i++)	mextra->kConsts.push_back(pow(mextra->Freqs[i],2)*mextra->RedMass[i]*tmp);
//
//}
//
//void CalcIRIntensity(){
//	//#ifdef VIBANA_DEBUG
//	System.out.printf("\n\nCalculate Intensity ";
//	//#endif
//	//double dipole_deriv[SIZE_MAX*3][SIZE_MAX*3];
//
//	//double delta=1e-5;
//	//EDipole_1stgrad(dipole_deriv,delta);
//	mextra->Intens.clear();
//	for(int k=0;k<nFreqs;k++){
//		double t[3];
//		for(int l=0;l<3;l++){
//			t[l]=0;
//			FOR_ATOMS_OF_MOL(a,m){
//				int j=a->GetIdx()-1;
//				int q=a->GetFormalCharge();
//				double mass=a->GetExactMass();
//				t[l]+=q/Math.sqrt(mass)*mextra->lCART[k][j*3+l]*4.79664/Math.sqrt(mextra->RedMass[k]);
//			}
//		}
//
//		double tmp=42.255*(t[0]*t[0]+t[1]*t[1]+t[2]*t[2]); //convert to KM/mol unit
//		mextra->Intens.push_back(tmp);
//	}
//	//#ifdef VIBANA_DEBUG
//	System.out.printf(" Done \n";
//	//#endif
//}
//
	private void MassweightedHessian(){
		for(int i=0;i<n_AtomOR3;i++){
			for(int l=i;l<n_AtomOR3;l++){				
				double mi=m.getMass(i/3);
				double ml=m.getMass(l/3);

				double fTmp=Hessian[l][i]/Math.sqrt(mi*ml);
				Hessian[i][l]=Hessian[l][i]=fTmp;
			}
		}
	}
///*
//double EDipole(double *mu){
//	double mu0=0;
//	for(int j=0;j<3;j++){
//		mu[j]=0;
//		for(int i=0;i<m->NumAtoms();i++){
//			int q=
//			switch(Nz[i]){
//				case HYDROGEN: q=+1;break;
//				case OXYGEN:   q=-2;break;
//			}
//			mu[j]+=q*x[i*3+j]*4.79664/Math.sqrt(getMass[Nz[i]]);
//		}
//		mu0+=SQR(mu[j]);
//	}
//	return Math.sqrt(mu0);
//}
//*/
///*
//void EDipole_1stgrad(double edipole_derive[][SIZE_MAX*3],double delta){
//	int k,l;
//	double mu1,mu0;
//	double mu0_v[3],mu1_v[3],delta_mu[3];
//	mu0=EDipole(mu0_v);
//	for(k=0;k<nFreqs;k++){
//		for(l=0;l<n_AtomOR3;l++)	x[l]+=lCART[k][l]*delta;
//		mu1=EDipole(mu1_v);
//		for(l=0;l<3;l++) edipole_derive[k][l]=(mu1_v[l]-mu0_v[l])/(delta);
//		for(l=0;l<n_AtomOR3;l++)	x[l]-=lCART[k][l]*delta;
//	}
//}
//*/
//
	private void Calc_tr_D(){
		double[] x=m.getCoords();

		//matrix P - the dot product of Coord_R_com and the corresponding eigenvector
		double[][] P=new double[3][m.getNAtoms()];

		for (int i = 0; i<3; i++){
			for (int j = 0; j< m.getNAtoms(); j++){
				P[i][j] = x[j*3+0]*I_vect[0][i] +
						x[j*3+1]*I_vect[1][i]+
						x[j*3+2]*I_vect[2][i];
			}
		}

		//change to the new base -- unit vector
		for (int i = 0; i<3; i++){
			for (int j = 0; j<3; j++){
				I_vect[i][j] = (i==j)?1:0;
			}
		}

		//the translation and rotation matrix
		int k = 0;
		for (int j = 0; j < n_AtomOR3; j = j + 3){

			double sqrtMassAtom=Math.sqrt(m.getMass(k));

			TR_D[0][j]   = sqrtMassAtom;
			TR_D[1][j+1] = sqrtMassAtom;
			TR_D[2][j+2] = sqrtMassAtom;

			TR_D[3][j]   = (P[1][k]*I_vect[0][2] - P[2][k]*I_vect[0][1])/sqrtMassAtom;
			TR_D[3][j+1] = (P[1][k]*I_vect[1][2] - P[2][k]*I_vect[1][1])/sqrtMassAtom;
			TR_D[3][j+2] = (P[1][k]*I_vect[2][2] - P[2][k]*I_vect[2][1])/sqrtMassAtom;

			TR_D[4][j]   = (P[2][k]*I_vect[0][0] - P[0][k]*I_vect[0][2])/sqrtMassAtom;
			TR_D[4][j+1] = (P[2][k]*I_vect[1][0] - P[0][k]*I_vect[1][2])/sqrtMassAtom;
			TR_D[4][j+2] = (P[2][k]*I_vect[2][0] - P[0][k]*I_vect[2][2])/sqrtMassAtom;

			TR_D[5][j]   = (P[0][k]*I_vect[0][1] - P[1][k]*I_vect[0][0])/sqrtMassAtom;
			TR_D[5][j+1] = (P[0][k]*I_vect[1][1] - P[1][k]*I_vect[1][0])/sqrtMassAtom;
			TR_D[5][j+2] = (P[0][k]*I_vect[2][1] - P[1][k]*I_vect[2][0])/sqrtMassAtom;

			k++;
		}

	double[] norm=new double[6];
	int jel_remove = 0;
	n_TrD = 6;
	for (int i = 0; i<6; i++){
		norm[i] = 0.0;
		for (int j = 0; j<n_AtomOR3; j++){
			norm[i] += TR_D[i][j]*TR_D[i][j];
		}
		norm[i] = Math.sqrt(norm[i]);
		if (norm[i] < 1e-9){
			jel_remove = 1;
			n_TrD--;
		}
	}

	//remove the line with zero elements
	k = 0;
	if (jel_remove == 1){
		for (int i = 0; i<6; i++){
			if (norm[i] > 1e-9){
				for (int j = 0; j<n_AtomOR3; j++){
					TR_D[k][j] = TR_D[i][j];
				}
				k++;
			}

		}
	}
}
	
	private int Calc_Sch_Ort(){
		//the orthogonalization procedure - Stabilized Gramm-Schmidt
		//first: 0-5 - they are linearly independent		
		for (int i = 0; i<n_TrD; i++){
			for (int j = 0; j<i; j++){
				double V_OR_D=MTools.DOTPRODUCT(TR_D[j],TR_D[i]);
				for (int k = 0; k<n_AtomOR3; k++){
					TR_D[i][k] -= V_OR_D*TR_D[j][k];
				}
			}

			double norm=MTools.DOTPRODUCT(TR_D[i],TR_D[i]);
			norm = Math.sqrt(norm);
			for (int k = 0; k<n_AtomOR3; k++){
				TR_D[i][k] = TR_D[i][k]/norm;
			}
		}

		//second: generate 6-3*n_Atom orthogonal vectors
		//complete the original set with a set of linerly
		//independent vectors and orthogonalize
		int cik = 0;
		for (int i = n_TrD; i<n_AtomOR3; i++){
			double norm = 0.0;
			while (Math.abs(norm) < 1e-9){
				for (int j = 0; j<n_AtomOR3; j++){
					//TR_D[i][j] = Unit_M[cik][j];
					TR_D[i][j] = (cik==j)?1:0;
				}

				for (int j = 0; j<n_AtomOR3; j++){
					for (int k = 0; k < i; k++)	{
						TR_D[i][j] -= TR_D[k][cik]*TR_D[k][j];
					}
				}

				norm=MTools.DOTPRODUCT(TR_D[i],TR_D[i]);
				norm = Math.sqrt(norm);

				//check for the linearly independence
				if ( Math.abs(norm) < 1e-6){
					cik++; //if not take another vector
					norm = 0.0;
				}

			}
			for (int k = 0; k<n_AtomOR3; k++){
				TR_D[i][k] = TR_D[i][k]/norm;
			}
			cik++;

			if(cik>=n_AtomOR3){
				n_TrD=-1;
				return -1;
			}
		}

		return 0;
	}


	
//
/////////////////////////////////////////////////////////////////////////////////////////
////Apply the transformation to the Hessian
/////////////////////////////////////////////////////////////////////////////////////////
private void Transform_Hessian(){
	//double X[SIZE_MAX*3][SIZE_MAX*3];
	double[][] X=new double[n_AtomOR3][n_AtomOR3];

	//multiply  TR_D with the Hessian
	for (int i = 0; i<n_AtomOR3; i++){
		for (int j = 0; j<n_AtomOR3; j++){
			X[i][j]=0;
			for (int k = 0; k<n_AtomOR3; k++){
				X[i][j] += TR_D[i][k]*Hessian[k][j];
			}
		}
	}

	//multiply the the result of the previous steep by the transpose of TR_D
	for (int i = 0; i<n_AtomOR3; i++){
		for (int j = 0; j<n_AtomOR3; j++){
			Hessian[i][j] = 0.0;
			for (int k = 0; k<n_AtomOR3; k++){
				Hessian[i][j] += X[i][k]*TR_D[j][k];
			}
		}
	}
}

	///////////////////////////////////////////////////////////////////////////////////////
	//Sort the frequencies -- Bubble sort algorithm
	///////////////////////////////////////////////////////////////////////////////////////
//	private void Sort_Freq(){
//		double dTmp;
//		for (int i=0;i<freqs.length;i++)
//			for (int j=nFreq-1; j > i; j--)
//				if(mextra->Freqs[j-1] > mextra->Freqs[j]){
//					swap(mextra->Freqs[j],mextra->Freqs[j-1]);
//					swap(mextra->RedMass[j],mextra->RedMass[j-1]);
//					swap(mextra->kConsts[j],mextra->kConsts[j-1]);
//					swap(mextra->Intens[j],mextra->Intens[j-1]);
//					for(int k=0;k<n_AtomOR3;k++){
//						swap(mextra->lCART[j][k],mextra->lCART[j-1][k]);
//					}
//				}
//	}

//
///*
//void WriteGauOut(ostream &os){
//	os<<"@IN "<<tag<<endl;
//	os<<" NAtoms=    "<<m->NumAtoms()<<endl;
//	os<<"                         Standard orientation:"<<endl;
//	os<<" ---------------------------------------------------------------------"<<endl;
//	os<<" Center     Atomic     Atomic              Coordinates (Angstroms)"<<endl;
//	os<<" Number     Number      Type              X           Y           Z"<<endl;
//	os<<" ---------------------------------------------------------------------"<<endl;
//	for(int i=0;i<m->NumAtoms();i++){
//		os<<format("%5d %10d             0     %11.6lf %11.6lf %11.6lf\n")%(i+1)%Nz[i]%x[i*3]%x[i*3+1]%x[i*3+2];
//	}
//	os<<endl;
//	if(bf)	os<<" Using "<<bf->equation()<<" empirical model !"<<endl;
//	os<<" Harmonic frequencies (cm**-1), IR intensities (KM/Mole), Raman scattering"<<endl;
//	os<<" activities (A**4/AMU), depolarization ratios for plane and unpolarized "<<endl;
//	os<<" incident light, reduced masses (AMU), force constants (mDyne/A),"<<endl;
//	os<<" and normal coordinates: "<<endl;
//
//	for(int i=0;i<nFreqs;i+=3){
//		os<<format("%22d%22d%22d")%i%(i+1)%(i+3)<<endl;
//		os<<"                     A                     A                     A"<<endl;
//		os<<format(" Frequencies --%11.4lf%23.4lf%23.4lf")%freq[i]%freq[i+1]%freq[i+2]<<endl;
//		os<<format(" Red. masses --%11.4lf%23.4lf%23.4lf")%RedMass[i]%RedMass[i+1]%RedMass[i+2]<<endl;
//		os<<format(" Frc consts  --%11.4lf%23.4lf%23.4lf")%kConsts[i]%kConsts[i+1]%kConsts[i+2]<<endl;
//		os<<format(" IR Inten    --%11.4lf%23.4lf%23.4lf")%IRIntensity[i]%IRIntensity[i+1]%IRIntensity[i+2]<<endl;
//		os<<" Atom AN      X      Y      Z        X      Y      Z        X      Y      Z"<<endl;
//		for(int j=0;j<m->NumAtoms();j++){
//			os<<format("%4d%4d")%j%Nz[j];
//			os<<format("  %7.2lf%7.2lf%7.2lf")
//					%lCART[i][j*3]%lCART[i][j*3+1]%lCART[i][j*3+2];
//			os<<format("  %7.2lf%7.2lf%7.2lf")
//					%lCART[i+1][j*3]%lCART[i+1][j*3+1]%lCART[i+1][j*3+2];
//			os<<format("  %7.2lf%7.2lf%7.2lf")
//					%lCART[i+2][j*3]%lCART[i+2][j*3+1]%lCART[i+2][j*3+2];
//			os<<endl;
//		}
//	}
//
//	os<<endl<<endl;
//	os<<" Principal axes and moments of inertia in atomic units:"<<endl;
//	os<<"                           1         2         3"<<endl;
//	os<<format("     EIGENVALUES --   %9.5lf %9.5lf %9.5lf")%I_eig[0]%I_eig[1]%I_eig[2]<<endl;
//	os<<format("			X         %9.5lf %9.5lf %9.5lf")%I_vect[0][0]%I_vect[0][1]%I_vect[0][2]<<endl;
//	os<<format("			Y         %9.5lf %9.5lf %9.5lf")%I_vect[1][0]%I_vect[1][1]%I_vect[1][2]<<endl;
//	os<<format("			Z         %9.5lf %9.5lf %9.5lf")%I_vect[2][0]%I_vect[2][1]%I_vect[2][2]<<endl;
//	os<<endl<<endl;
//}
//*/
//


}
