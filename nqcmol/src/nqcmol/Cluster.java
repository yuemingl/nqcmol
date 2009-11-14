/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.logging.*;
import nqcmol.tools.MTools;

/**
 *
 * @author nqc
 */
public class Cluster implements Cloneable{
	public Cluster(){
	};

	public Cluster(int nAtoms){
		setNAtoms(nAtoms);
	};

	public Cluster(Cluster cluster_){
		Get(cluster_);
	}
	
	@Override
	public Object clone() {
           return new Cluster(this);           
        }

	public void Get(Cluster src){
		//System.out.println(" Come here");
		setNAtoms(src.getNAtoms());
		//System.out.println(" Come here");
		setCoords(src.getCoords());
		//System.out.println(" Come here");
		setGradient(src.getGradient());
		//System.out.println(" Come here");
		setHessian(src.getHessian());
		//System.out.println(" Come here");
		setFreqs(src.getFreqs());
		//System.out.println(" Come here");
		setAtomicNumber(src.getAtomicNumber());
		//System.out.println(" Come here");
		tag=src.getTag();
		rmsGrad=src.getRmsGrad();		
	}

	
	protected static Logger logger=Logger.getLogger(Cluster.class.getName());

	/* Maximum number of atomic number
	 */
	static final int cTypeMax=15;

	/**
	 * Symbols of chemical elements
	 */
	static final String[] cElements={"$",
	"H","He","Li","Be","B","C","N","O","F","Ne", //1st and 2nd round
	"Na","Mg","Al","Si","P","S","Cl","Ar", 	//3rd round
	"K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr"};

	static final String[] cAtomicNo={"0",
	"1","2","3","4","5","6","7","8","9","10",
	"11","12","13","14","15","16","17","18",
	"19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36"};
	
	static final double[] cMass={0,
	1.00783,4,6,8,10,12,14,15.99491,18,20,
	22.989,24.305,27,28,31,32,35.453,40,
	39,49,44,47.867,51,52,55,55.845,59,58.693,63.546,65.38,69.723,72.64,75,79,80,83.798};


	//================== properties of Clusters
	/**
	 * Number of atoms of cluster
	 */
	protected int nAtoms = 0;

	/**
	 * Get the value of nAtoms
	 *
	 * @return the value of nAtoms
	 */
	public int getNAtoms() {
		return nAtoms;
	}

	/**
	 * Set the value of nAtoms
	 *
	 * @param nAtoms new value of nAtoms
	 */
	public void setNAtoms(int nAtoms) {
		if(this.nAtoms!=nAtoms){
			this.nAtoms = nAtoms;
			ncoords=nAtoms*3;
			coords=new double[ncoords];
			gradient=new double[ncoords];
			hessian=new double[ncoords][ncoords];
			Nz=new int [nAtoms];			
		}
	}

	/**
	 * Coordinates of clusters
	 */
	protected double[] coords;

	/**
	 * Get the value of coords
	 *
	 * @return the value of coords
	 */
	public double[] getCoords() {
		return coords;
	}

	/**
	 * Set the value of coords
	 *
	 * @param coords new value of coords
	 */
	public void setCoords(double[] coords) {
		assert ncoords<=coords.length;
		System.arraycopy(coords,0,this.coords,0,ncoords);
	}

	/**
	 * Get the value of coords at specified index
	 *
	 * @param index
	 * @return the value of coords at specified index
	 */
	public double getCoords(int index) {
		return this.coords[index];
	}

	/**
	 * Set the value of coords at specified index.
	 *
	 * @param index
	 * @param newCoords new value of coords at specified index
	 */
	public void setCoords(int index, double newCoords) {
		this.coords[index] = newCoords;
	}

	protected int ncoords = 0;

	public int getNcoords() {
		return ncoords;
	}

	protected double energy;

	/**
	 * Get the value of energy
	 *
	 * @return the value of energy
	 */
	public double getEnergy() {
		return energy;
	}

	/**
	 * Set the value of energy
	 *
	 * @param energy new value of energy
	 */
	public void setEnergy(double energy) {
		this.energy = energy;
	}

	protected double[] gradient;

	/**
	 * Get the value of gradient
	 *
	 * @return the value of gradient
	 */
	public double[] getGradient() {
		return gradient;
	}

	/**
	 * Set the value of gradient
	 *
	 * @param gradient new value of gradient
	 */
	public void setGradient(double[] gradient_) {
		assert ncoords<=gradient_.length;
		System.arraycopy(gradient_,0,this.gradient,0,ncoords);
	}

	/**
	 * Get the value of gradient at specified index
	 *
	 * @param index
	 * @return the value of gradient at specified index
	 */
	public double getGradient(int index) {
		return this.gradient[index];
	}

	/**
	 * Set the value of gradient at specified index.
	 *
	 * @param index
	 * @param newGrad new value of gradient at specified index
	 */
	public void setGradient(int index, double newGrad) {
		this.gradient[index] = newGrad;
	}

	protected double[][] hessian;

	public double[][] getHessian(){
		return hessian;
	}

	public void setHessian(double[][] hessian_){
		assert ncoords<=hessian_.length;
		for(int i=0;i<ncoords;i++){
			assert ncoords<=hessian_[i].length;
			System.arraycopy(hessian_[i],0,this.hessian[i],0,ncoords);
		}
	}

	protected double[] freqs = null;

	/**
	 * Get the value of freqs
	 *
	 * @return the value of freqs
	 */
	public double[] getFreqs() {
		return freqs;
	}

	/**
	 * Set the value of freqs
	 *
	 * @param freqs new value of freqs
	 */
	public void setFreqs(double[] freqs) {
		this.freqs = freqs.clone();
	}

	/**
	 * Get the value of freqs at specified index
	 *
	 * @param index
	 * @return the value of freqs at specified index
	 */
	public double getFreqs(int index) {
		return this.freqs[index];
	}

	protected double rmsGrad=0;

	/**
	 * Get the value of rmsGrad
	 *
	 * @return the value of rmsGrad
	 */
	public double getRmsGrad() {
		return rmsGrad;
	}

	/**
	 * Set the value of rmsGrad
	 *
	 * @param rmsGrad new value of rmsGrad
	 */
	public void setRmsGrad(double rmsGrad) {
		this.rmsGrad = rmsGrad;
	}


	protected int tag;

	/**
	 * Get the value of tag
	 *
	 * @return the value of tag
	 */
	public int getTag() {
		return tag;
	}

	/**
	 * Set the value of tag
	 *
	 * @param tag new value of tag
	 */
	public void setTag(int tag) {
		this.tag = tag;
	}

	/**
	 * Array holds  atomic numbers. Access by atom's index
	 */
	protected int[] Nz;

	/**
	 * Get the value of Nz
	 * @return the value of Nz
	 */
	public int[] getAtomicNumber() {
		return Nz;
	}

	/**
	 * Get the value of Nz
	 *
	 * @param index of atom
	 * @return the value of Nz
	 */
	public int getAtomicNumber(int index) {
		return Nz[index];
	}

	/**
	 * Set array of atomic number
	 * @param iType new value of Nz
	 */
	public void setAtomicNumber(int[] Nz_) {
		assert nAtoms<=Nz_.length;
		System.arraycopy(Nz_,0,this.Nz,0,nAtoms);
	}


	/**
	 * Set the value of Nz
	 * @param index index of atom
	 * @param iType new value of Nz
	 */
	public void setAtomicNumber(int index,int type) {
		Nz[index] = type;
	}


	//double[] a,alpha; //lattice parameters
	//double USRsig[12];	//for USR
	//int iCorType; //direct or fractional coordinates	
	
	int[] nType=new int[cTypeMax];

	/**
	 * @return Total mass of clusters
	 */
	public double getTotalMass() {
		double TotalMass=0;
		for(int i=0;i<nAtoms;i++)
			TotalMass+=getMass(i);
		return TotalMass;
	}

	/**
	 * @param i index of atom
	 * @return getMass of an atom
	 */
	public double getMass(int i){
		return cMass[Nz[i]];
	}

	//======================== Methods


	public void Clear(){
		nAtoms=0;	tag=0;
		energy=0; rmsGrad=0;
	}

	public void ClearGradients(){
		for(int i=0;i<ncoords;i++) gradient[i]=0;
	}

	public int getAtomicNumberFromSymbol(String atomsymbol){
		int k=0;
		for(int i=0;i< cElements.length ;i++){
			//System.out.print(cElements[i]+" \n");
			if(cElements[i].contentEquals(atomsymbol)){
				k=i; break;
			}
		}
		//System.out.printf(" %s %d \n",cElements[k],k);
		return k;
	}

	public void R(int i,double[] dr){
		for(int l=0;l<3;l++) dr[l]=coords[i*3+l];
	}


	public void dR(int i,int j,double[] dr){
			for(int l=0;l<3;l++)	dr[l]=coords[j*3+l]-coords[i*3+l];
	};

	public void dnR(int i,int j,double[] dr){ //!< normalised dR(i,j)
			dR(i,j,dr);	double d=distance(i,j);
			//cout<<" l= "<<d<<" "<<dr[0]<<" "<<dr[1]<<" "<<dr[2]<<endl;
			for(int l=0;l<3;l++)	dr[l]/=d;
		};

	public double distance(int i,int j){//!< distance between i and j atoms
			double[] vec=new double[3];	dR(i,j,vec);
			return  Math.sqrt(MTools.DOTPRODUCT(vec,vec));
	}

	public double distance(int i,double[] d){//!< distance between i and a given point
			double[] vec={coords[i*3],coords[i*3+1],coords[i*3+2]};
			if(d!=null){ vec[0]-=d[0];vec[1]-=d[1];vec[2]-=d[2];}
			return  Math.sqrt(MTools.DOTPRODUCT(vec,vec));
		}

	public double angle(int ow,int hw1,int hw2){ //!< angle between (hw1,ow,hw2)
			double[] l1=new double[3];
			double[] l2=new double[3];
			dR(hw1,ow,l1);
			dR(hw2,ow,l2);
			return Math.acos(MTools.DOTPRODUCT(l1,l2)/Math.sqrt(MTools.DOTPRODUCT(l1,l1)*MTools.DOTPRODUCT(l2,l2)));
		}

	public double angle_dgr(int ow,int hw1,int hw2){ //!< angle between (hw1,ow,hw2) in degree
			return angle(ow,hw1,hw2)*180.0/Math.PI;
	}

	public double torsion_dgr(int i1,int i2,int i3,int i4){ //!< torsion angle of four atoms i1-i2-i3-i4 in degree
			return torsion(i1,i2,i3,i4)*180.0/Math.PI;
	}


	public double torsion(int i1,int i2,int i3,int i4){
		double t=0;
		double[] b1=new double[3];
		double[] b2=new double[3];
		double[] b3=new double[3];
		double[] c1=new double[3];
		double[] c2=new double[3];
		double[] c3=new double[3];

		dR(i2,i1,b1);//b1 = a - b;
		dR(i3,i2,b2);//b2 = b - c;
		dR(i4,i3,b3);//b3 = c - d;

		MTools.CROSSPRODUCT(c1,b1,b2); //c1 = cross(b1,b2);
		MTools.CROSSPRODUCT(c2,b2,b3); //c2 = cross(b2,b3);
		MTools.CROSSPRODUCT(c3,c1,c2); //c3 = cross(c1,c2);

		if (MTools.DOTPRODUCT(c1,c1) * MTools.DOTPRODUCT(c2,c2) < 0.001) {
		  t = 0.0;
		  return t;
		}

		t = Math.acos(MTools.DOTPRODUCT(c1,c2)/Math.sqrt(MTools.DOTPRODUCT(c1,c1)*MTools.DOTPRODUCT(c2,c2)));
		if (MTools.DOTPRODUCT(b2,c3) > 0.0)    t = -t;

		return (t);
	}

	public void RotateX(double theta){
		double cost=Math.cos(theta),sint=Math.sin(theta);
		double xt,yt,zt;
		for(int i=0;i<nAtoms;i++){
			xt=coords[i*3]; yt=coords[i*3+1]; zt=coords[i*3+2];
			coords[i*3+0]=  xt;
			coords[i*3+1]=  yt*cost + zt*sint;
			coords[i*3+2]= -yt*sint + zt*cost;
		}
	}

	public void RotateY(double theta){
		double cost=Math.cos(theta),sint=Math.sin(theta);
		double xt,yt,zt;
		for(int i=0;i<nAtoms;i++){
			xt=coords[i*3]; yt=coords[i*3+1]; zt=coords[i*3+2];
			coords[i*3+0]=  xt*cost +zt*sint;
			coords[i*3+1]=  yt;
			coords[i*3+2]= -xt*sint + zt*cost;
		}
	}

	public void RotateZ(double theta){
		double cost=Math.cos(theta),sint=Math.sin(theta);
		double xt,yt,zt;
		for(int i=0;i<nAtoms;i++){
			xt=coords[i*3]; yt=coords[i*3+1]; zt=coords[i*3+2];
			coords[i*3+0]=  xt*cost + yt*sint;
			coords[i*3+1]= -xt*sint + yt*cost;
			coords[i*3+2]=  zt;
		}
	}

	public void RotateAxis(double theta,double[] axis){
		for(int i=0;i<nAtoms;i++) RotateAxis_atom(i,theta,axis);
	}

	public void Translate(double[] vec,double direction){
		for(int i=0;i<nAtoms;i++) Translate_atom(i,vec,direction);
	}

	public void Transform(double[][] D){
		for(int i=0;i<nAtoms;i++) Transform_atom(i,D);
	}


	public void RotateX_atom(int i,double theta){
		double cost=Math.cos(theta),sint=Math.sin(theta);
		double xt,yt,zt;

		xt=coords[i*3]; yt=coords[i*3+1]; zt=coords[i*3+2];
		coords[i*3+0]=  xt;
		coords[i*3+1]=  yt*cost + zt*sint;
		coords[i*3+2]= -yt*sint + zt*cost;
	}

	public void RotateY_atom(int i,double theta){
		double cost=Math.cos(theta),sint=Math.sin(theta);
		double xt,yt,zt;

		xt=coords[i*3]; yt=coords[i*3+1]; zt=coords[i*3+2];
		coords[i*3+0]=  xt*cost +zt*sint;
		coords[i*3+1]=  yt;
		coords[i*3+2]= -xt*sint + zt*cost;
	}

	public void RotateZ_atom(int i,double theta){
		double cost=Math.cos(theta),sint=Math.sin(theta);
		double xt,yt,zt;
		xt=coords[i*3]; yt=coords[i*3+1]; zt=coords[i*3+2];
		coords[i*3+0]=  xt*cost + yt*sint;
		coords[i*3+1]= -xt*sint + yt*cost;
		coords[i*3+2]=  zt;
	}

	public void RotateAxis_atom(int i,double theta,double axis[]){
		double X,Y,Z,sum;
		sum=Math.sqrt(MTools.SQR(axis[0])+MTools.SQR(axis[1])+MTools.SQR(axis[2]));
		X=axis[0]/sum;	Y=axis[1]/sum;	Z=axis[2]/sum;
		//construct rotation matrix
		double c,s,t;
		c=Math.cos(theta);	s=Math.sin(theta);	t=1-c;
		double[][] R={
			{t*X*X+c  , t*X*Y+s*Z, t*X*Z-s*Y,	0},
			{t*X*Y-s*Z, t*Y*Y+c,   t*Y*Z+s*X,	0},
			{t*X*Z+s*Y, t*Z*Y-s*X, t*Z*Z+c,		0},
			{0,			0,			0,			1}};

	  int l,m,n;
	  double[] xt=new double[4],xtp=new double[4];

	  for(l=0;l<3;l++)	xt[l]=coords[i*3+l];	xt[3]=1;
	  MTools.MAT2_MUL_VEC(xtp,R,xt);
	  for(l=0;l<3;l++){
		  coords[i*3+l]=xtp[l];
	  }
	}

	public void RotateAxisWithPoint_atom(int i,double theta,double axis[],double point[]){
		Translate_atom(i,point,-1);
		RotateAxis_atom(i,theta,axis);
		Translate_atom(i,point,1);
	}

	public void Translate_atom(int i,double[] vec,double direction){
		coords[i*3+0]+=direction*vec[0];
		coords[i*3+1]+=direction*vec[1];
		coords[i*3+2]+=direction*vec[2];
	}

	public void Transform_atom(int i, double[][] D){
		int l,m; double[] v=new double[3];
		R(i,v);
		for(l=0;l<3;l++){
			coords[i*3+l]=0;
			for(m=0;m<3;m++) coords[i*3+l]+=D[l][m]*v[m];
		}
	}

	public void RotateAxis_mol(int[] mk,double theta,double[] axis){
		double X,Y,Z,sum;
		sum=Math.sqrt(MTools.SQR(axis[0])+ MTools.SQR(axis[1]) + MTools.SQR(axis[2]));
		X=axis[0]/sum;	Y=axis[1]/sum;	Z=axis[2]/sum;
		//construct rotation matrix
		double c,s,t;
		c=Math.cos(theta);	s=Math.sin(theta);	t=1-Math.cos(theta);
		double[][] R={
			{t*X*X+c  , t*X*Y+s*Z, t*X*Z-s*Y,	0},
			{t*X*Y-s*Z, t*Y*Y+c,   t*Y*Z+s*X,	0},
			{t*X*Z+s*Y, t*Z*Y-s*X, t*Z*Z+c,		0},
			{0,			0,			0,			1}};

		  int l,m,n;
		  double[] xt=new double[4],xtp=new double[4];
		  for(int i=0;i<mk.length;i++){
			  for(l=0;l<3;l++)	xt[l]=coords[mk[i]*3+l];	xt[3]=1;
			  MTools.MAT2_MUL_VEC(xtp,R,xt);
			  for(l=0;l<3;l++)	coords[mk[i]*3+l]=xtp[l];
		  }
	}

	public void Translate_mol(int[] mk,double[] vec,double direction){
		for(int i=0;i<mk.length;i++){
			coords[mk[i]*3]+=direction*vec[0];
			coords[mk[i]*3+1]+=direction*vec[1];
			coords[mk[i]*3+2]+=direction*vec[2];
		}
	}

	public void RotateAxisWithPoint_mol(int[] mk,double theta,double[] axis,double[] point){
		for(int i=0;i<mk.length;i++){
			RotateAxisWithPoint_atom(mk[i],theta,axis,point);
		}
	}

	public void Center(){
		int i,k;
		//fix to origin
		double[] rc=new double[3];
		getMassCenter(rc);
		for(i=0;i<nAtoms;i++)
			for(k=0;k<3;k++)
				coords[i*3+k]-=rc[k];
	}


	public void getMassCenter(double[] rc){
		int i,k;
		for(k=0;k<3;k++) rc[k]=0;
		double mass=0;
		for(i=0;i<nAtoms;i++){
			mass+=cMass[Nz[i]];
			for(k=0;k<3;k++)
				rc[k]+=cMass[Nz[i]]*coords[i*3+k];
		}
		for(k=0;k<3;k++) rc[k]/=mass;
	}


	public void getInertiaTensor(double[][] I){
		MTools.MAT2_EQU_NUM(I,0);
		for(int i=0;i<nAtoms;i++){
			double mass_i=getMass(i);
			
			I[0][0]+= mass_i*(MTools.SQR(coords[i*3+1])+MTools.SQR(coords[i*3+2]));
			I[1][1]+= mass_i*(MTools.SQR(coords[i*3+0])+MTools.SQR(coords[i*3+2]));
			I[2][2]+= mass_i*(MTools.SQR(coords[i*3+0])+MTools.SQR(coords[i*3+1]));

			I[0][1]+=-mass_i*coords[i*3+0]*coords[i*3+1];
			I[0][2]+=-mass_i*coords[i*3+0]*coords[i*3+2];
			I[1][2]+=-mass_i*coords[i*3+1]*coords[i*3+2];

			System.err.printf(" i12 = %f \n",I[1][2]);
		}
		I[1][0]=I[0][1];
		I[2][1]=I[1][2];
		I[2][0]=I[0][2];
		//		System.err.printf(" I \n");
		//		MTools.PrintArray(I);
	}

	//======================== Read/Write Method

	public boolean Read(Scanner scanner,String format){
		boolean isReadable=false;
		if(format.contentEquals("xyz"))		isReadable=ReadXYZ(scanner);
		
		return isReadable;
	}
	protected boolean ReadXYZ(Scanner scanner){
			Clear();

			if(!scanner.hasNextInt()) return false;
			int nAtoms_=scanner.nextInt();
			//System.out.printf(" Here nAtoms = %d \n",nAtoms_);
			if (nAtoms_ <= 0) {	return false; }
			setNAtoms(nAtoms_);

			scanner.nextLine();

			String line=scanner.nextLine();

			StringTokenizer tokenizer;
			tokenizer = new StringTokenizer(line, "\t ,;"); //read no of atoms

			String info="";
				if (line.contains("@IN")) {
					tag=(int) energy;
					//for my format
					info = tokenizer.nextToken();		energy = Double.parseDouble(info);
					info = tokenizer.nextToken();		rmsGrad = Double.parseDouble(info);
				}else{
					info= tokenizer.nextToken();     	energy = Double.parseDouble(info);
				}
				//System.out.printf(" Here tag=%d Energy = %f rmsGrad=%f\n",tag,energy,rmsGrad);


				//System.out.printf(" Here tag=%d Energy = %f rmsGrad=%f\n",Nz[0],energy,rmsGrad);

				//gradient=new double[nAtoms*3];
				for (int i = 0; i < nAtoms; i++) {
					line = scanner.nextLine();
					//System.out.printf(" Here %s nAtoms = %d \n",line,nAtoms_);
					if (line == null) return false;

					if (line.startsWith("#") && line.length() > 1) {
						i--; // a comment line does not count as an atom
					} else {
						tokenizer = new StringTokenizer(line, "\t ,;");
						int fields = tokenizer.countTokens();
						if (fields < 4) {
							return false;
						} else {
							String atomtype = tokenizer.nextToken();
							Nz[i] =  getAtomicNumberFromSymbol(atomtype);
							//System.out.printf(" %s %d \n",atomtype,i);

							if(Nz[i]<=0) return false;

							//System.out.printf(" Here i=%d Nz=%d x = %f y= %f z=%f\n",i,Nz[i],coords[i*3+0],coords[i*3+1],coords[i*3+2]);
							nType[Nz[i]]++;
							//System.out.printf(" Here i=%d Nz=%d x = %f y= %f z=%f\n",i,Nz[i],coords[i*3+0],coords[i*3+1],coords[i*3+2]);
							coords[i*3+0] = Double.parseDouble(tokenizer.nextToken());
							coords[i*3+1] = Double.parseDouble(tokenizer.nextToken());
							coords[i*3+2] =  Double.parseDouble(tokenizer.nextToken());
							//System.out.printf(" Here i=%d Nz=%d x = %f y= %f z=%f\n",i,Nz[i],coords[i*3+0],coords[i*3+1],coords[i*3+2]);
							if (fields >= 7){
								gradient[i*3+0] = Double.parseDouble(tokenizer.nextToken());
								gradient[i*3+1] = Double.parseDouble(tokenizer.nextToken());
								gradient[i*3+2] =  Double.parseDouble(tokenizer.nextToken());
							}
						}
					}
				}
			
		//System.out.printf(" Here \n");
		return true;
	}

	public void Write(Writer writer,String format) {
		if(nAtoms>0){
			try {
				if(format.contentEquals("xyz")){
					WriteXYZ(writer);
				}
				writer.flush();
			} catch (IOException ex) {
				Logger.getLogger(Cluster.class.getName()).log(Level.SEVERE, null, ex);
			} 
		}
	}

	public void Write(OutputStream os,String format) {
		BufferedWriter writer=new BufferedWriter(new OutputStreamWriter(os));
		Write(writer,format);
	}

	protected void WriteXYZ(Writer writer) throws  IOException{
        String s1 = "" + nAtoms + "\n";		
        writer.append(s1);
		//System.out.printf(" Here %s \n",s1);

        s1=String.format("%d %1.10f %f ",tag,energy,rmsGrad);
		if(freqs!=null)
			if(freqs.length>0){
				s1+=String.format("@FREQ %d",freqs.length);
				for(int i=0;i<freqs.length;i++)
					s1+=String.format(" %1.3f",freqs[i]);
			}
		s1+="% @IN\n";
		writer.append(s1);
		//System.out.printf(" Here %s \n",s1);

        // Loop through the atoms and write them out:
        for(int i=0;i<nAtoms;i++){
			s1=cElements[Nz[i]]+String.format(" %1.6f %1.6f %1.6f\n",coords[i*3],coords[i*3+1],coords[i*3+2]);
			writer.append(s1);
			//System.out.printf(" Here %s \n",s1);
        }
	}
}
