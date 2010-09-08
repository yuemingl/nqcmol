/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.logging.*;
import nqcmol.tools.MTools;
import nqcmol.tools.XmlWriter;

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
        double[] p=src.getUSRsig();
        for (int i=0; i<12; i++){	USRsig[i] = p[i];}

        energy=src.getEnergy();
		tag=src.getTag();
		rmsGrad=src.getRmsGrad();	
	}

	
	protected static Logger logger=Logger.getLogger(Cluster.class.getName());

	/* Maximum number of atomic number
	 */
	

	/**
	 * Symbols of chemical elements
	 */
	protected static final String[] cElements={"$",
	"H","He","Li","Be","B","C","N","O","F","Ne", //1st and 2nd round
	"Na","Mg","Al","Si","P","S","Cl","Ar", 	//3rd round
	"K","Ca","Sc","Ti","V","Cr","Mn","Fe","Co","Ni","Cu","Zn","Ga","Ge","As","Se","Br","Kr",//4th round
    "Rb","Sr","Y","Zr","Nb","Mo","Tc","Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","Te","I","Xe",//5th round
    "Cs","Ba", //6th round
    "La","Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb", //Lanthaniod
    "Lu","Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At","Rn", //6th round
    "Fr","Ra",//7th round
    "Ac","Th","Pa","U","Np","Pu","Am","Cm","Bk","Cf","Es","Fm","Md","No",//Actonid
    "103","104","105","106","107","108","109","110","111","112","113","114","115","116","116","118", //7th round
    "OW ","HW1","HW2","LP1","LP2" //extra or dummy atom
};

	protected static final String[] cAtomicNo={"0",
	"1","2","3","4","5","6","7","8","9","10",//1st and 2nd round
	"11","12","13","14","15","16","17","18",//3rd round
	"19","20","21","22","23","24","25","26","27","28","29","30","31","32","33","34","35","36",//4th round
    "37","38","39","40","41","42","43","44","45","46","47","48","49","50","51","52","53","54",//5th round
    "55","56",//6th round
    "57","58","59","60","61","62","63","64","65","66","67","68","69","70",//Lanthaniod
    "71","72","73","74","75","76","77","78","79","80","81","82","83","84","85","86",//6th round
    "87","88",//7th round
    "89","90","91","92","93","94","95","96","97","98","99","100","101","102",//Actonid
    "103","104","105","106","107","108","109","110","111","112","113","114","115","116","116","118", //7th round
    "119","120","121","122","123" //extra or dummy atom
    };


	protected static final double[] cMass={0,
	1.00783,4,6,8,10,12,14,15.99491,18,20,
	22.989,24.305,27,28,31,32,35.453,40,
	39,49,44,47.867,51,52,55,55.845,59,58.693,63.546,65.38,69.723,72.64,75,79,80,83.798};

	protected static final int cTypeMax=cElements.length;

	public final static String[] format={"xyz","g03","g03c","car","xml","int","pdb"};

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
			pairwise=new PairwiseType[nAtoms][nAtoms];
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
		if(freqs!=null)
			this.freqs = freqs.clone();
        else this.freqs=null;
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

	protected double[] reducedMass = null;

	/**
	 * Get the value of reducedMass
	 *
	 * @return the value of reducedMass
	 */
	public double[] getReducedMass() {
		return reducedMass;
	}

	/**
	 * Set the value of reducedMass
	 *
	 * @param reducedMass new value of reducedMass
	 */
	public void setReducedMass(double[] reducedMass) {
		this.reducedMass = reducedMass;
	}

	/**
	 * Get the value of reducedMass at specified index
	 *
	 * @param index
	 * @return the value of reducedMass at specified index
	 */
	public double getReducedMass(int index) {
		return this.reducedMass[index];
	}

	/**
	 * Set the value of reducedMass at specified index.
	 *
	 * @param index
	 * @param newReducedMass new value of reducedMass at specified index
	 */
	public void setReducedMass(int index, double newReducedMass) {
		this.reducedMass[index] = newReducedMass;
	}


	protected double[][] normalModeVectors=null;

	/**
	 * Get the value of normalModeVectors
	 *
	 * @return the value of normalModeVectors
	 */
	public double[][] getNormalModeVectors() {
		return normalModeVectors;
	}

	/**
	 * Set the value of normalModeVectors
	 *
	 * @param normalModeVectors new value of normalModeVectors
	 */
	public void setNormalModeVectors(double[][] normalModeVector) {
		this.normalModeVectors = normalModeVector.clone();
	}

	/**
	 * Get the value of normalModeVectors at specified index
	 *
	 * @param index
	 * @return the value of normalModeVectors at specified index
	 */
	public double[] getNormalModeVectors(int index) {
		return this.normalModeVectors[index];
	}

	/**
	 * Set the value of normalModeVectors at specified index.
	 *
	 * @param index
	 * @param newNormalModeVector new value of normalModeVectors at specified index
	 */
	public void setNormalModeVectors(int index, double[] newNormalModeVector) {
		this.normalModeVectors[index] = newNormalModeVector.clone();
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


	String tag;

	/**
	 * Get the value of tag
	 *
	 * @return the value of tag
	 */
	public String getTag() {
		return tag;
	}

	/**
	 * Set the value of tag
	 *
	 * @param tag new value of tag
	 */
	public void setTag(String tag) {
		this.tag = tag;
	}

	/**
	 * Array holds  atomic numbers. Access by atom's index
	 */
	protected int[] Nz;

	/**
	 * @param index: index of ato,
	 * @return Atomic Symbol of atom
	 */
	public String getAtomicSymbol(int index) {
		return cElements[Nz[index]];
	}

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


	/**
	 * Set the value of Nz
	 * @param index index of atom
	 * @param x new coordinates of atom
	 */
	public void setAtomicCoords(int index,double x,double y,double z) {
		coords[index*3+0] = x;
		coords[index*3+1] = y;
		coords[index*3+2] = z;
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
	 * @param index index of atom
	 * @return getMass of an atom
	 */
	public double getMass(int i){
		return cMass[Nz[i]];
	}

    	/**
	 * @return Total charge of clusters
	 */
	public int getTotalCharge() {
		int TotalCharge=0;
		for(int i=0;i<nAtoms;i++)
			TotalCharge+=getCharge(i);
		return TotalCharge;
	}

	/**
	 * @param index index of atom
	 * @return getMass of an atom
	 */
	public int getCharge(int i){
        int charge = 0;
        switch(Nz[i]){
            case 1: charge= 1;break;
            case 8: charge=-2;break;
            case 7: charge=-1;break;
        }
		return charge;
	}


    public boolean isNAN(){
        boolean answer=false;

        answer=Double.isNaN(energy);
        for(int i=0;i<coords.length;i++) answer=Double.isNaN(coords[i]);
        return answer;
    }
	//======================== Methods


	public void Clear(){
		MTools.VEC_EQU_NUM(coords, 0);
		MTools.VEC_EQU_NUM(gradient, 0);
		MTools.VEC_EQU_NUM(nType, 0);
		MTools.VEC_EQU_NUM(Nz, 0);

		ncoords=0;
		nAtoms=0;	tag="0";
		energy=0; rmsGrad=0;
	}

	public void ClearGradients(){
		for(int i=0;i<ncoords;i++) gradient[i]=0;
	}

	public int getAtomicNumberFromSymbol(String atomsymbol){
		int k=0;
		for(int i=0;i< cElements.length ;i++){
			//System.out.print(cElements[index]+" \n");
			if(cElements[i].contentEquals(atomsymbol)){
				k=i; break;
			}
		}

		if(k==0){
			k=Math.min(Integer.parseInt(atomsymbol),cTypeMax-1);
		}
		//System.out.printf(" %s %d \n",cElements[k],k);
		return k;
	}

	public void R(int i,double[] dr){
        if(dr==null) dr=new double[3];
		for(int l=0;l<3;l++) dr[l]=coords[i*3+l];
	}


	public void dR(int i,int j,double[] dr){
			for(int l=0;l<3;l++)	dr[l]=coords[j*3+l]-coords[i*3+l];
	};

	public void dnR(int i,int j,double[] dr){ //!< normalised dR(index,j)
			dR(i,j,dr);	double d=distance(i,j);
			//cout<<" l= "<<d<<" "<<dr[0]<<" "<<dr[1]<<" "<<dr[2]<<endl;
			for(int l=0;l<3;l++)	dr[l]/=d;
		};

	public double distance(int i,int j){//!< distance between index and j atoms
			double[] vec=new double[3];	dR(i,j,vec);
			return  Math.sqrt(MTools.DOTPRODUCT(vec,vec));
	}

	public double distance(int i,double[] d){//!< distance between index and a given point
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
			return Math.toDegrees(angle(ow,hw1,hw2));
	}

	public double torsion_dgr(int i1,int i2,int i3,int i4){ //!< torsion angle of four atoms i1-i2-i3-i4 in degree
			return Math.toDegrees(torsion(i1,i2,i3,i4));
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

        double b2abs=Math.sqrt(MTools.DOTPRODUCT(b2, b2));
        t=-Math.atan2(b2abs*MTools.DOTPRODUCT(b1,c2), MTools.DOTPRODUCT(c1,c2));//see wikipedia

//		if (MTools.DOTPRODUCT(c1,c1) * MTools.DOTPRODUCT(c2,c2) < 0.001) {
//		  t = 0.0;
//		  return t;
//		}
//
//		t = Math.acos(MTools.DOTPRODUCT(c1,c2)/Math.sqrt(MTools.DOTPRODUCT(c1,c1)*MTools.DOTPRODUCT(c2,c2)));
//		if (MTools.DOTPRODUCT(b2,c3) > 0.0)    t = -t;

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

	public void CorrectOrder(){
		for(int i=0;i<nAtoms-1;i++)
			for(int j=nAtoms-1;j>i;j--)
				if(Nz[j]>Nz[j-1]){
					SwapAtom(j,j-1);
				}
	}

	public void SwapAtom(int i,int j){ //!< swapping index and j atom including force
			int tempI,l;
			double tempD;
			tempI=Nz[i];	Nz[i]=Nz[j];	Nz[j]=tempI;

			for(l=0;l<3;l++){
				tempD=coords[i*3+l];
				coords[i*3+l]=coords[j*3+l];
				coords[j*3+l]=tempD;

				tempD=gradient[i*3+l];
				gradient[i*3+l]=gradient[j*3+l];
				gradient[j*3+l]=tempD;
			}
	}

	//========================= extra get method
	public void getMassCenter(double[] rc){
		int i,k;
		for(k=0;k<3;k++) rc[k]=0;
		double totalmass=0;
		for(i=0;i<nAtoms;i++){
            double iM=1.0;
            if(Nz[i]<cMass.length) iM=cMass[Nz[i]];
			totalmass+=iM;
			for(k=0;k<3;k++)
				rc[k]+=iM*coords[i*3+k];
		}
		for(k=0;k<3;k++) rc[k]/=totalmass;
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

			//System.err.printf(" i12 = %f \n",I[1][2]);
		}
		I[1][0]=I[0][1];
		I[2][1]=I[1][2];
		I[2][0]=I[0][2];
		//		System.err.printf(" I \n");
		//		MTools.PrintArray(I);
	}

	/**
	 * @return the number of hydrogen atoms
	 */
	public int getHydrogenNum(){
		return nType[1];
	}

	public int getNonHydrogenNum(){
		return nAtoms-nType[1];
	}

	/**
	 * @param i index of atom
	 * @return true if atom is hydrogen
	 */
	public boolean IsHydrogen(int i){
        boolean answer=false;
        if((Nz[i]==1)||(Nz[i]==120)||(Nz[i]==121)||(Nz[i]==122)||(Nz[i]==123)) answer=true;
		return answer;
	}

	/**
	 * @param i index of atom
	 * @return true if atom is not hydrogen
	 */
	public boolean IsNonHydrogen(int i){
		return !IsHydrogen(i);
	}

	public String getFormula(){
		String answer="";
		for(int i=nType.length-1;i>=0;i--){
			if(nType[i]>0){
				answer+=cElements[i]+Integer.toString(nType[i]);
			}
		}
		return answer;
	}

    /**
     * replace a number of atoms (from index) as the whole atoms from mol
     * @param index starting position to replace
     * @param mol
     */
	public void replaceMolecule(int index,Cluster mol){
		for(int j=0;j<mol.getNAtoms();j++){
			coords[(index+j)*3+0]=mol.getCoords(j*3+0);
			coords[(index+j)*3+1]=mol.getCoords(j*3+1);
			coords[(index+j)*3+2]=mol.getCoords(j*3+2);
			Nz[index+j]=mol.getAtomicNumber(j);
		}
		for(int i=0;i<cTypeMax;i++) nType[i]=0;
		for(int i=0;i<nAtoms;i++) nType[Nz[i]]++;
	}

    /**
     * erase an atom, the number of atoms will be deduced by 1
     * @param index of atom to be erased
     */
	public void eraseAtom(int index){
        double[] newCoords=new double[coords.length-3];
        MTools.eraseElementsFromArray(coords,newCoords, index*3,3);
        coords=newCoords;

        int[] newNz = new int[Nz.length-1];
        MTools.eraseElementsFromArray(Nz,newNz, index,1);
        Nz=newNz;
        
//        System.arraycopy(coords, 0, newCoords, 0, index*3);
//        if (coords.length != index*3) {
//            System.arraycopy(coords, (index + 1)*3, newCoords, index*3, coords.length - 3*index );
//        }
//        coords=newCoords;
//
//        int[] newNz = new double[(nAtoms-1)*3];
//
//        System.arraycopy(Nz, 0, newNz, 0, index*3);
//        if (coords.length != index*3) {
//            System.arraycopy(coords, (index + 1)*3, newCoords, index*3, coords.length - 3*index );
//        }
//        coords=newCoords;
//
//		for(int j=0;j<mol.getNAtoms();j++){
//			coords[(index+j)*3+0]=mol.getCoords(j*3+0);
//			coords[(index+j)*3+1]=mol.getCoords(j*3+1);
//			coords[(index+j)*3+2]=mol.getCoords(j*3+2);
//			Nz[index+j]=mol.getAtomicNumber(j);
//		}

        nAtoms--;
		for(int i=0;i<cTypeMax;i++) nType[i]=0;
		for(int i=0;i<nAtoms;i++) nType[Nz[i]]++;
	}

	//======================== Bond 2
	/**
	 * Pair-wise bond type
	 * NONE: no connection
	 * HYD_NEAR: covalent hydrogen bond
	 * HYD_FAR: hydrogen bond
	 * SINGLEBOND: single bond
	 */
	public enum PairwiseType { NONE, HYD_NEAR, HYD_FAR, SINGLEBOND, DOUBLEBOND };

	protected PairwiseType[][] pairwise;

	/**
	 * Calculate Bond2 matrix
	 * @return true if connected, false if disconnected
	 */
	public boolean getPairwiseBond(){
		for(int i=0;i<nAtoms;i++){
			for(int j=0;j<nAtoms;j++)
				pairwise[i][j]=pairwise[j][i]=getPairwiseBond(i,j);
			pairwise[i][i]=PairwiseType.NONE;
		}


//		System.err.println(" We are here");
//		for(int index=0;index<nAtoms;index++){
//			for(int j=0;j<nAtoms;j++)
//				System.err.printf("%d",pairwise[index][j].ordinal());
//			//pairwise[index][index]=PairwiseType.NONE;
//			System.err.println("");
//		}

		return IsConnected();
	}

	/**
     * @param i,j index of atoms
     * @return Bond type defined in PairwiseType enum
	 */
	public PairwiseType getPairwiseBond(int i,int j){
		int ni=Nz[i];
		int nj=Nz[j];
		String s=cElements[ni]+cElements[nj];
		double dmin=1000;	double dmax=0;
		PairwiseType type=PairwiseType.NONE;

		double d=distance(i,j);

		if(s.contentEquals("HO")){
				dmin=1.20; dmax=2.25; type=PairwiseType.HYD_FAR;	if((d>=dmin)&&(d<=dmax))	return type;
				dmin=0.30; dmax=1.20; type=PairwiseType.HYD_NEAR;	if((d>=dmin)&&(d<=dmax))	return type;
		}
		if(s.contentEquals("HF")){
				dmin=1.40; dmax=2.25; type=PairwiseType.HYD_FAR;	if((d>=dmin)&&(d<=dmax))	return type;
				dmin=0.30; dmax=1.40; type=PairwiseType.HYD_NEAR;	if((d>=dmin)&&(d<=dmax))	return type;
		}
		else if(s.contentEquals("OH")){ dmin=0.30; dmax=1.20; type=PairwiseType.HYD_NEAR;}
		else if(s.contentEquals("OO")){ dmin=2.00; dmax=3.20; type=PairwiseType.SINGLEBOND;}
		else if(s.contentEquals("FH")){ dmin=0.20; dmax=1.4; type=PairwiseType.HYD_NEAR;}
		else if(s.contentEquals("FF")){ dmin=1.00; dmax=3.2; type=PairwiseType.SINGLEBOND; }
		else if(s.contentEquals("FO")||s.contentEquals("OF")){ dmin=1.00; dmax=3.2; type=PairwiseType.SINGLEBOND; }
		else if(s.contentEquals("CO")||s.contentEquals("OC")){ dmin=1.0; dmax=1.5; type=PairwiseType.DOUBLEBOND;}
		else if(s.contentEquals("CH")){ dmin=0.2; dmax=1.2; type=PairwiseType.HYD_NEAR;}

		if((d>=dmin)&&(d<=dmax)){
			return type;
		}else return PairwiseType.NONE;
	}

    public int[] getMolecule(int i){
        Vector<Integer> tmp=new Vector<Integer>();
        tmp.add(i);
        for(int k=0;k<nAtoms;k++) if(k!=i){
           if(pairwise[i][k]==PairwiseType.HYD_NEAR){
                tmp.add(k);
           }
        }
        int[] answer=new int[tmp.size()];
        for(int k=0;k<tmp.size();k++) answer[k]=tmp.elementAt(k);
        return answer;
    }

	//======================= Structure query
	/**
	 * Check whether two atoms index and j are connected
	 * @return true if connected, false if not
	 */
	public boolean IsConnected(int i,int j){
		PairwiseType type=getPairwiseBond(i,j);
		return (type!=PairwiseType.NONE);
	}

	/**
	 * Check whether any atom of cluster is connected to at least one atom 
	 * @return true if connected, false if not
	 */
	public boolean IsConnected(){
		boolean[][] d=new boolean[nAtoms][nAtoms];
		for(int i=0;i<nAtoms;i++){
			for(int j=0;j<nAtoms;j++)
				d[j][i]=IsConnected(i,j);
			
			d[i][i]=true;
		}

		int n=nAtoms;
		for(int i=0;i<n;i++){
			for(int j=0;j<n;j++)
				if((i!=j)&&(d[i][j])){//if connected,jump to jth row and flip
					for(int l=0;l<n;l++) if(d[i][l]) d[j][l]=true;
				}
		}
		for(int i=0;i<n;i++)
			for(int j=0;j<n;j++) if(!d[i][j]) return false;

		return true;
	}


	public static final int cNumCountMax=7;
	protected int[] coordNumCount=new int[cNumCountMax];

	/**
	 * Get the coordination number of Nonhydrogen atoms based on the pairwise table
	 * @param bUpdate if true, pairwise table will be updated
	 * @return the value of coordNumCount
	 */
	public int[] getCoordNumCount(boolean bUpdate) {
		if(bUpdate) getPairwiseBond();

		for(int i=0;i<coordNumCount.length;i++) coordNumCount[i]=0;

		for(int i=0;i<nAtoms;i++)	if(!IsHydrogen(i)){
			int count=getCoordNum(i);
			
			if(count> coordNumCount.length){
				//logger.log(Level.SEVERE,String.format(" A strange structure appears ! Coordination number %d is too large! ",count));
			}else coordNumCount[count]++;
		}

		return coordNumCount;
	}

	/**
	 * Return the number of rings based on Cauchy formula: n(edges) - n(vertices) + 1
	 * @param bUpdate if true, pairwise table will be updated
	 * @return the number of rings
	 */
	public int getNumberOfSmallestRingByCauchyFormula(boolean bUpdate){
		getCoordNumCount(bUpdate);
		int d=0;
		for(int i=0;i<coordNumCount.length;i++){
			d+=coordNumCount[i]*i;
		}
		//return d/2 - getNonHydrogenNum() +1;
		return d/2 - nType[8]-nType[9] +1;
	}

	/**
	 * Get the value of coordNumCount at specified index
	 *
	 * @param index
	 * @return the value of coordNumCount at specified index
	 */
	public int getCoordNumCount(int index) {
		return this.coordNumCount[index];
	}


	/**
	 * Get the value of coordNumCount
	 * @param index
	 * @return String of coordNumCount, starting from 0 coord until nonzero coord.
	 */
	public int[] getCoordNumCount() {
		return coordNumCount;
	}

	/**
	 * Get coordination number of nonhydgen atoms
	 * @param index of atom
	 * @return the value of coordNumCount at specified index
	 */
	public int getCoordNum(int index) {
		int count=0;
		for(int j=0;j<nAtoms;j++)	
			//if(!IsHydrogen(j) && (pairwise[index][j]!=PairwiseType.NONE) )	count++;
			if( ( (Nz[j]==8)||(Nz[j]==9)) && (pairwise[index][j]==PairwiseType.SINGLEBOND) )	count++;
		return count;
	}

	/**
	 * @param bUpdate if true, pairwise table will be updated
	 * @return morphology (in string)
	 */
	public String getMorphology(boolean bUpdate){
		if(bUpdate) getPairwiseBond();
		if(!IsConnected()) return "Disconnected";
		else{
			int nRings=getNumberOfSmallestRingByCauchyFormula(bUpdate);

			//System.err.printf(" nRings = %d , nH = %d\n",nRings,getHydrogenNum() );
			//MTools.PrintArray(coordNumCount);

			switch(nRings){
				case 0:
					int sum=0;
					for(int i=3;i<coordNumCount.length;i++) sum+=coordNumCount[i];
					if(sum>0) return "Treelike";
					else return "Linear";
				case 1: return "SingleRing";
				case 2:	return "DoubleRing";
			}

			if(nRings>=3){
					if(coordNumCount[3]+coordNumCount[4]==getNonHydrogenNum()) return "Cage";
					else return "MultiRing";
			}
		}
		return "Undefined";
	}

	/**
	 * @return  Compactness measured by deviation of distance
	 */
	public double getCompactness(){
		int k=0;
		double[] dev={0,0,0};
		double[] center=new double[3];
		getMassCenter(center);

		for(int i=0;i<nAtoms;i++) if(IsNonHydrogen(i)){
			double[] dr=new double[3];
			dr[0]=Math.pow(coords[i*3+0]-center[0],2);
			dr[1]=Math.pow(coords[i*3+1]-center[1],2);
			dr[2]=Math.pow(coords[i*3+2]-center[2],2);
			
			MTools.VEC_PLUS_VEC(dev,dev,dr,1.0,1.0);
			k++;
		}

		if(k!=0)
			MTools.VEC_MUL_NUM(dev, dev, 1.0/k);
		double P=Math.sqrt(dev[0]+dev[1]+dev[2]);
		return P;
	}
	
	//======================== Similarity index
	protected double[] USRsig=new double[12];

    public double[] getUSRsig(){
        return USRsig;
    }

	/**
	 *
	 * @param p: Cluster will be compared. Supposed that USR signature of p and (this) is both computed already
	 * @return the similarity between (this) and p.
	 */
	public double CalcUSRSimilarity(Cluster p){
		double S =0;
			for (int i=0; i<12; i++) {
				S += Math.abs(USRsig[i] - p.USRsig[i]);
			}
			S /= 12.0;
			S = 1.0/(1.0+S);
			return S;
	}

	public void CalcUSRsignature(){	
	//clear USRsig
	for (int i=0; i<12; i++){	USRsig[i] = 0;}

	if(nAtoms==2){ USRsig[0]=distance(0,1); return ;}

	double[] c={0,0,0};
	//compute centroid (C)

	for (int i=0; i<nAtoms; i++) {
		c[0] += coords[i*3+0]/nAtoms;
		c[1] += coords[i*3+1]/nAtoms;
		c[2] += coords[i*3+2]/nAtoms;
	}
	getMassCenter(c);
	
	//find distances to the centroid
	double[] distToCentroid=new double[nAtoms];
	double maxDist = -1;
	double minDist = 1e14;
	int minItem = 0;
	int maxItem = 0;
	for (int i=0; i<nAtoms; i++) {
		double d = distance(i,c);
		if (d < minDist) {
			minDist = d;		minItem = i;
		}
		if (d > maxDist) {
			maxDist = d;		maxItem = i;
		}
		distToCentroid[i]=d;
	}

	//find closest to centroid (A)
	//A is minItem	
	c[0] = coords[minItem*3+0];
	c[1] = coords[minItem*3+1];
	c[2] = coords[minItem*3+2];

	double[] distToA=new double[nAtoms];
	for (int i=0; i<nAtoms; i++) {
		double d = distance(i,c);
		distToA[i]=d;
	}

	//find furthest from the centroid (B)
	//B is maxItem and update maxitem for C
	c[0] = coords[maxItem*3+0];
	c[1] = coords[maxItem*3+1];
	c[2] = coords[maxItem*3+2];

	double[] distToB=new double[nAtoms];
	maxDist = -1;
	for (int i=0; i<nAtoms; i++) {
		double d = distance(i,c);
		if (d > maxDist) {
			maxDist = d;	maxItem = i;
		}
		distToB[i]=d;
	}

	//find oxy furthest from B (C)
	//C is in maxItem now
	c[0] = coords[maxItem*3+0];
	c[1] = coords[maxItem*3+1];
	c[2] = coords[maxItem*3+2];

	double[] distToC=new double[nAtoms];
	for (int i=0; i<nAtoms; i++) {
		double d = distance(i,c);
		distToC[i]=d;
	}

	//compute USRsigs from each distribution


	//compute means of the distributions
	for (int i=0; i<nAtoms; i++) {
		USRsig[0] += distToCentroid[i];
		USRsig[3] += distToA[i];
		USRsig[6] += distToB[i];
		USRsig[9] += distToC[i];
	}
	USRsig[0] /= (double) nAtoms;
	USRsig[3] /= (double) nAtoms;
	USRsig[6] /= (double) nAtoms;
	USRsig[9] /= (double) nAtoms;

	//compute variances of the distributions
	for (int i=0; i<nAtoms; i++) {
		USRsig[1] += Math.pow(distToCentroid[i] - USRsig[0], 2.0);
		USRsig[4] += Math.pow(distToA[i] - USRsig[3], 2.0);
		USRsig[7] += Math.pow(distToB[i] - USRsig[6], 2.0);
		USRsig[10] += Math.pow(distToC[i] - USRsig[9], 2.0);
	}
	USRsig[1] /= (double) (nAtoms-1);
	USRsig[4] /= (double) (nAtoms-1);
	USRsig[7] /= (double) (nAtoms-1);
	USRsig[10] /= (double) (nAtoms-1);

	//compute skewness of the distributions
	for (int i=0; i<nAtoms; i++) {
		USRsig[2] += Math.pow(distToCentroid[i] - USRsig[0], 3.0);
		USRsig[5] += Math.pow(distToA[i] - USRsig[3], 3.0);
		USRsig[8] += Math.pow(distToB[i] - USRsig[6], 3.0);
		USRsig[11] += Math.pow(distToC[i] - USRsig[9], 3.0);
	}
	USRsig[2] /= (double) (nAtoms-1);
	USRsig[5] /= (double) (nAtoms-1);
	USRsig[8] /= (double) (nAtoms-1);
	USRsig[11] /= (double) (nAtoms-1);
	USRsig[2] /= (double) Math.pow(Math.sqrt(USRsig[1]),3.0);
	USRsig[5] /= (double) Math.pow(Math.sqrt(USRsig[4]),3.0);
	USRsig[8] /= (double) Math.pow(Math.sqrt(USRsig[7]),3.0);
	USRsig[11] /= (double) Math.pow(Math.sqrt(USRsig[10]),3.0);

	//MTools.PrintArray(USRsig);
	}
	//======================== Read/Write Method
	
	public boolean Read(Scanner scanner,String format){
		boolean isReadable=false;
		//if(format.contentEquals("xyz"))		isReadable=ReadXYZ(scanner);
		if(format.contentEquals("int"))	isReadable=ReadInt(scanner);
		else if(format.contentEquals("g03"))	isReadable=ReadGauOut(scanner);
		else if(format.contentEquals("g03c"))	isReadable=ReadGauOut_C(scanner);
        else if(format.contentEquals("car"))	isReadable=ReadCAR(scanner);
		else
			isReadable=ReadXYZ(scanner);
		
		return isReadable;
	}

	protected boolean ReadXYZ(Scanner scanner){
		//System.out.printf(" Here nAtoms = %d \n",nAtoms);
		Clear();

		if(!scanner.hasNextInt()) return false;
		int nAtoms_=scanner.nextInt();
		//System.err.printf(" Here nAtoms = %d \n",nAtoms_);
		if (nAtoms_ <= 0) {	return false; }
		setNAtoms(nAtoms_);

		scanner.nextLine();

		String line=scanner.nextLine(); line=line.trim();

		StringTokenizer token;
		token = new StringTokenizer(line, "\t ,;"); //read no of atoms

		String info="";
			if (line.contains("@IN")) {
				//for my format
				info = token.nextToken();		tag = info;
				info = token.nextToken();		energy = Double.parseDouble(info);
				info = token.nextToken();		rmsGrad = Double.parseDouble(info);
				while(token.hasMoreTokens()){
					if(token.nextToken().contentEquals("@FREQ")){
						info = token.nextToken();
						freqs= new double[Integer.parseInt(info)];
						for(int i=0;i<freqs.length;i++){
							info = token.nextToken();		freqs[i] = Double.parseDouble(info);
						}
					}
				}
			}else if(!line.isEmpty()){
				info= token.nextToken();     	tag=info;//energy = Double.parseDouble(info);
			}



			//System.err.printf(" Here tag=%d Energy = %f rmsGrad=%f\n",tag,energy,rmsGrad);


			//System.out.printf(" Here tag=%d Energy = %f rmsGrad=%f\n",Nz[0],energy,rmsGrad);

			//gradient=new double[nAtoms*3];
			for (int i = 0; i < nAtoms; i++) {
				line = scanner.nextLine();
				//System.out.printf(" Here %s nAtoms = %d \n",line,nAtoms_);
				if (line == null) return false;

				if (line.startsWith("#") && line.length() > 1) {
					i--; // a comment line does not count as an atom
				} else {
					token = new StringTokenizer(line, "\t ,;");
					int fields = token.countTokens();
					if (fields < 4) {
						return false;
					} else {
						String atomtype = token.nextToken();
						Nz[i] =  getAtomicNumberFromSymbol(atomtype);
						//System.err.printf(" AtomType=%s Nz=%d index=%d \n",atomtype,Nz[index],index);

						if(Nz[i]<=0) return false;

						//System.out.printf(" Here index=%d Nz=%d x = %f y= %f z=%f\n",index,Nz[index],coords[index*3+0],coords[index*3+1],coords[index*3+2]);
						nType[Nz[i]]++;
						//System.out.printf(" Here index=%d Nz=%d x = %f y= %f z=%f\n",index,Nz[index],coords[index*3+0],coords[index*3+1],coords[index*3+2]);
						coords[i*3+0] = Double.parseDouble(token.nextToken());
						coords[i*3+1] = Double.parseDouble(token.nextToken());
						coords[i*3+2] =  Double.parseDouble(token.nextToken());
						//System.err.printf(" Here index=%d Nz=%d x = %f y= %f z=%f\n",index,Nz[index],coords[index*3+0],coords[index*3+1],coords[index*3+2]);
						if (fields >= 7){
							gradient[i*3+0] = Double.parseDouble(token.nextToken());
							gradient[i*3+1] = Double.parseDouble(token.nextToken());
							gradient[i*3+2] =  Double.parseDouble(token.nextToken());
						}
					}
				}
			}
			
		//System.out.printf(" Here \n");
		return true;
	}

	protected boolean ReadInt(Scanner scanner){
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
				info = tokenizer.nextToken();		tag = info;
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
                    //System.out.printf(" Number of fieldds = %d \n",fields);
					if (fields < 1) {
						return false;
					} else {
						String atomtype = tokenizer.nextToken();
						Nz[i] =  getAtomicNumberFromSymbol(atomtype);
						//System.out.printf(" %s %d \n",atomtype,Nz[i]);

						if(Nz[i]<=0) { nAtoms=0;	return false;}

						//System.out.printf(" Here index=%d Nz=%d x = %f y= %f z=%f\n",index,Nz[index],coords[index*3+0],coords[index*3+1],coords[index*3+2]);
						nType[Nz[i]]++;
                        double dst,ang,tor;
                        int i1,i2,i3;
                        double[] v1=new double[3];
                        double[] v2=new double[3];
                        double[] v3=new double[3];
                        double[] n=new double[3];
                        double[] nn=new double[3];

                        switch(i){
                          case	0: coords[0]=coords[1]=coords[2]=0;
                          break;
                          case	1:
                                i1=Integer.parseInt(tokenizer.nextToken());
                                dst=Double.parseDouble(tokenizer.nextToken());
                                coords[i*3+0]=dst;coords[i*3+1]=coords[i*3+2]=0;
                          break;
                          case 2:
                            i1=Integer.parseInt(tokenizer.nextToken());
                            dst=Double.parseDouble(tokenizer.nextToken());
                            i2=Integer.parseInt(tokenizer.nextToken());
                            ang=Double.parseDouble(tokenizer.nextToken());
                            i1--;i2--;
                            //cout<<i1<<" "<<dst<<" "<<i2<<" "<<ang<<" "<<i3<<" "<<tor<<endl;
                            ang=Math.toRadians(ang);
                            coords[i*3+0] = -Math.cos(ang)*dst + coords[1*3+0];
                            coords[i*3+1] =  Math.sin(ang)*dst;
                            coords[i*3+2] =  0;
                            if (i1== 0)
                                 coords[i*3+0] = Math.cos(ang)*dst;

                            //cout<<dst*cos(ang)<<" "<<dst*sin(ang)<<endl;
//                            coords[i*3+0]=coords[i1*3+0]+dst*Math.cos(ang);
//                            coords[i*3+1]=0;
//                            coords[i*3+2]=coords[i1*3+2]+dst*Math.sin(ang);
                          break;
                          default:
                            i1=Integer.parseInt(tokenizer.nextToken());
                            dst=Double.parseDouble(tokenizer.nextToken());
                            i2=Integer.parseInt(tokenizer.nextToken());
                            ang=Double.parseDouble(tokenizer.nextToken());
                            i3=Integer.parseInt(tokenizer.nextToken());
                            tor=Double.parseDouble(tokenizer.nextToken());
                            i1--;i2--;i3--;
                            ang=Math.toRadians(ang);
                            tor=Math.toRadians(tor);
                            //System.out.printf("%d %f %d %f %d %f\n",i1,dst,i2,ang,i3,tor);

                            double[] r1=new double[3]; R(i1,r1);
                            double[] r2=new double[3]; R(i2,r2);
                            double[] r3=new double[3]; R(i3,r3);
                            double[] r4=MTools.zmatrixToCartesian(dst, r1, ang, r2, tor, r3);
                            for(int j=0;j<3;j++) coords[i*3+j]=r4[j];

//                            dR(i2,i1,v1);
//                            dR(i3,i1,v2);
//                            MTools.CROSSPRODUCT(n,v1,v2);
//                            MTools.CROSSPRODUCT(nn,v1,n);
//                            MTools.NORMALIZE(n);
//                            MTools.NORMALIZE(nn);
//                            //VEC_MUL_NUM(n,n,3,-sin(tor));
//                            //VEC_MUL_NUM(nn,nn,3,cos(tor));
//                            MTools.VEC_PLUS_VEC(v3,n,nn,-Math.sin(tor),Math.cos(tor));
//                            MTools.NORMALIZE(v3);
//                            MTools.VEC_MUL_NUM(v3,v3,dst * Math.sin(ang));
//                            MTools.NORMALIZE(v1);
//                            MTools.VEC_MUL_NUM(v1,v1,dst * Math.cos(ang));
//
//                            for(int j=0;j<3;j++) coords[i*3+j]=coords[i1*3+j]+v3[j]-v1[j];
                          break;
                        }
					}
				}
			}

		//System.out.printf(" Here \n");
		return true;
	}

     protected boolean ReadCAR(Scanner scanner){
        Clear();
		int flag=0;
		tag="0";
        return false;
     }

	protected boolean ReadGauOut(Scanner scanner){
		Clear();
		int flag=0;
		tag="0";

		Vector<Double> freqs_tmp=new Vector<Double>();
		Vector<Double> reducedMass_tmp=new Vector<Double>();
		Vector<double[]> normalmode_tmp=new Vector<double[]>();

		while(scanner.hasNext()){
			String s=scanner.nextLine();
			//System.out.println(s);
			//read atomic number and coordinates

			//if( (s.contains("Input orientation")|| s.contains("Standard orientation") && (flag==0) )){
			if( ( s.contains("Input orientation") && (flag==0) ) || s.contains("Standard orientation") ){
				scanner.nextLine();scanner.nextLine();scanner.nextLine();scanner.nextLine();
				int i=0;
				Vector coords_tmp=new Vector();
				Vector Nz_tmp=new Vector();
				while(scanner.hasNext()){
					s=scanner.nextLine();
					//System.err.println(s);
					if(s.charAt(1)=='-') break;

					StringTokenizer tokenizer = new StringTokenizer(s, "\t ,;"); //read no of atoms
					tokenizer.nextToken();
					Nz_tmp.add(Integer.parseInt(tokenizer.nextToken()));
					tokenizer.nextToken();
					coords_tmp.add(Double.parseDouble(tokenizer.nextToken()));
					coords_tmp.add(Double.parseDouble(tokenizer.nextToken()));
					coords_tmp.add(Double.parseDouble(tokenizer.nextToken()));
					i++;
				}
				setNAtoms(i);
				for(i=0;i<nAtoms;i++){
					Nz[i]=(Integer)Nz_tmp.get(i);
					nType[Nz[i]]++;
					coords[i*3+0]=(Double)coords_tmp.get(i*3+0);
					coords[i*3+1]=(Double)coords_tmp.get(i*3+1);
					coords[i*3+2]=(Double)coords_tmp.get(i*3+2);
				}
				//Write(System.err,"xyz");

				flag|=1;//indicate coordinates have been collected
			}

			//read the energy
			if(s.contains("E(")&&s.contains("SCF Done")){
				int pos1=s.indexOf("=")+1;
				int pos2=s.indexOf("A.U.");
				String substr=s.substring(pos1,pos2);
				energy=Double.parseDouble(substr);
				flag|=2; //indicate energy has been collected
			}


			if ( s.contains("EUMP2 =")  ){	// wow, energy
				String substr=s.substring(s.indexOf("EUMP2 =")+7);
				substr=substr.replace('D','E');
                //System.out.println(substr);
				energy=Double.parseDouble(substr);
				flag|=2; //indicate energy has been collected
			}

			//read the RMS
			if(s.contains("RMS     Force")){
				String substr=s.substring(25,37);
				rmsGrad=Double.parseDouble(substr);
			}

			if(flag==3){ 	flag=0;	}

			if(s.contains("Frequencies --")){
				String substr=s.substring(16);

				StringTokenizer token = new StringTokenizer(substr, "\t ,;"); 
                int nToken=token.countTokens();
                for(int m=0;m<nToken;m++)
                    freqs_tmp.add(Double.parseDouble(token.nextToken()));

				substr=scanner.nextLine().substring(16);//getline(is,s);	sscanf(s.c_str(), "Red. masses  --  %lf ",&(dtmp[0]));
				token = new StringTokenizer(substr, "\t ,;");
                for(int m=0;m<nToken;m++)
                    reducedMass_tmp.add(Double.parseDouble(token.nextToken()));
				

				while(scanner.hasNext()){
					substr=scanner.nextLine();//getline(is,s);	sscanf(s.c_str(), "Frc consts  --  %lf ",&(dtmp[0]));
					if( substr.contains("Atom AN") || substr.contains("Atom  AN") ) break;
                    //System.err.printf("%s\n",substr);
				}

				//while(s.find("Atom")!=0) getline(fin,s);
				double[][] v1=new double[nToken][ncoords];
				for(int k=0;k<nAtoms;k++){
					substr=scanner.nextLine();
					token = new StringTokenizer(substr, "\t ,;");
					token.nextToken();
					token.nextToken();
                    for(int m=0;m<nToken;m++)
                        for(int j=0;j<3;j++)
                            v1[m][k*3+j]=Double.parseDouble(token.nextToken());
					
					
				}

                for(int m=0;m<nToken;m++)
                    normalmode_tmp.add(v1[m]);
			}
		}

		if(normalmode_tmp.size()>0){
			normalModeVectors=new double[normalmode_tmp.size()][];
			for(int i=0;i<normalmode_tmp.size();i++)
				normalModeVectors[i]=normalmode_tmp.get(i).clone();
		}

		freqs=MolExtra.ConvertToDoubleArray(freqs_tmp);
		reducedMass=MolExtra.ConvertToDoubleArray(reducedMass_tmp);
		return (nAtoms!=0);
	}

	protected boolean ReadGauOut_C(Scanner scanner){
		Clear();
		int flag=0;

		while(scanner.hasNext()){
			String s=scanner.nextLine();
			//cout<<s<<endl;
			//read atomic number and coordinates

			if( (s.contains("Input orientation")|| s.contains("Standard orientation") && (flag==0) )
					   ){
				scanner.nextLine();scanner.nextLine();scanner.nextLine();scanner.nextLine();
				int i=0;
				Vector coords_tmp=new Vector();
				Vector Nz_tmp=new Vector();
				while(scanner.hasNext()){
					s=scanner.nextLine();
					//System.err.println(s);
					if(s.charAt(1)=='-') break;

					StringTokenizer tokenizer = new StringTokenizer(s, "\t ,;"); //read no of atoms
					tokenizer.nextToken();
					Nz_tmp.add(Integer.parseInt(tokenizer.nextToken()));
					tokenizer.nextToken();
					coords_tmp.add(Double.parseDouble(tokenizer.nextToken()));
					coords_tmp.add(Double.parseDouble(tokenizer.nextToken()));
					coords_tmp.add(Double.parseDouble(tokenizer.nextToken()));
					i++;
				}
				setNAtoms(i);
				for(i=0;i<nAtoms;i++){
					Nz[i]=(Integer)Nz_tmp.get(i);
					nType[Nz[i]]++;
					coords[i*3+0]=(Double)coords_tmp.get(i*3+0);
					coords[i*3+1]=(Double)coords_tmp.get(i*3+1);
					coords[i*3+2]=(Double)coords_tmp.get(i*3+2);
				}
				//Write(System.err,"xyz");

				flag|=1;//indicate coordinates have been collected
			}

			//read the energy
			if( ( (flag & 1) ==1) && s.contains("E(")&& s.contains("SCF Done") ){
				int pos1=s.indexOf("=")+1;
				int pos2=s.indexOf("A.U.");
               // System.err.printf("pos 1 = %d pos 2 = %d\n",pos1,pos2);
				String substr=s.substring(pos1,pos2);
				energy=Double.parseDouble(substr);
				
				flag|=2; //indicate energy has been collected

                //System.err.printf("ERHF=%f - %d\n",energy,flag);
			}


			if ( ((flag&1)==1)&& s.contains("EUMP2")  ){	// wow, energy
                String substr=s.substring(s.indexOf("EUMP2")+8);
				substr=substr.replace('D','E');
               //System.out.println(substr);
				energy=Double.parseDouble(substr);
               
				flag|=2; //indicate energy has been collected
                // System.err.printf("EUMP2=%f- %d\n",energy,flag);
			}

			//read the RMS
			if( (flag==3) && s.contains("RMS     Force")){
				String substr=s.substring(25,37);
				rmsGrad=Double.parseDouble(substr);
                flag=0; return true;
			}
            
			if((flag==3) && s.contains("Population")){ 	flag=0; return true;}
		}
		return false;
	}

	public void Write(Writer writer,String format) {
		if(nAtoms>0){
			try {				
				if(format.contentEquals("cml")){
					WriteCML(writer);
				}else	if(format.contentEquals("int")){
					WriteInt(writer);
				}else	if(format.contentEquals("eps")){
					WriteEPS(writer);
                }else	if(format.contentEquals("car")){
					WriteCAR(writer);
                }else	if(format.contentEquals("pdb")){
					WritePDB(writer);
				}else{
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
		if(nAtoms>0){
			String s1 = "" + nAtoms + "\n";
			writer.append(s1);
			//System.out.printf(" Here %s \n",s1);

			s1=String.format("%s %1.10f %f ",tag,energy,rmsGrad)+"% ";
			if(freqs!=null)
				if(freqs.length>0){
					s1+=String.format("@FREQ %d",freqs.length);
					for(int i=0;i<freqs.length;i++)
						s1+=String.format(" %1.3f",freqs[i]);
				}
			s1+=" @IN\n";
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

	protected void WriteCML(Writer writer) throws IOException{
		BufferedWriter w=new BufferedWriter(writer);
		XmlWriter xmlwriter=new XmlWriter(w);
		xmlwriter.writeEntity("molecule");
		xmlwriter.writeAttribute("id","m"+tag);
		xmlwriter.writeEntity("name");
		xmlwriter.writeText(getFormula());
		xmlwriter.endEntity();

		xmlwriter.writeEntity("scalar");
		xmlwriter.writeAttribute("title", "energy");
		xmlwriter.writeNormalText(Double.toString(energy));
		xmlwriter.endEntity();
		
		xmlwriter.writeEntity("atomArray");
		for(int i=0;i<nAtoms;i++){
			xmlwriter.writeEntity("atom");
			xmlwriter.writeAttribute("id","a"+Integer.toString(i));
			xmlwriter.writeAttribute("elementType",cElements[Nz[i]]);
			xmlwriter.writeAttribute("x3",Double.toString(coords[i*3]));
			xmlwriter.writeAttribute("y3",Double.toString(coords[i*3+1]));
			xmlwriter.writeAttribute("z3",Double.toString(coords[i*3+2]));			
			xmlwriter.endEntity();
			
        }
		xmlwriter.endEntity();

		xmlwriter.endEntity();
		xmlwriter.flush();
	}

	protected void WriteInt(Writer writer) throws IOException{
		if(nAtoms>=1){
			String s1 = "" + nAtoms + "\n";
			writer.append(s1);
			//System.out.printf(" Here %s \n",s1);

			s1=String.format("%s %1.10f %f ",tag,energy,rmsGrad);
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
						//if(nAtom>=2) os<<Elements[Nz[1]]<<" 1 "<<distance(1,0)<<endl;
			//if(nAtom>=3) os<<Elements[Nz[2]]<<" "<<distance(2,1)<<" "<<angle_dgr(2,1,0)<<endl;
			for(int i=0;i<nAtoms;i++){
				double dist,distmin;

				int j1min=-1; distmin=1e9;
				for(int j=i-1;j>=0;j--){
					dist=distance(j,i);		if(dist<distmin){ distmin=dist; j1min=j;}
				}

				int j2min=-1;	distmin=1e9;
				for(int j=i-1;j>=0;j--) if(j!=j1min){
					dist=distance(j,j1min);		if(dist<distmin){ distmin=dist; j2min=j;}
				}

				int j3min=-1;	distmin=1e9;
				for(int j=i-1;j>=0;j--) if((j!=j1min)&&(j!=j2min)){
					dist=distance(j,j2min);		if(dist<distmin){ distmin=dist; j3min=j;}
				}

				s1=cElements[Nz[i]];

				if(j1min!=-1) s1+=String.format(" %3d %8.3f ",(j1min+1),distance(i,j1min));
				if(j2min!=-1) s1+=String.format(" %3d %8.3f ",(j2min+1),angle_dgr(j1min,i,j2min));
				if(j3min!=-1) s1+=String.format(" %3d %8.3f ",(j3min+1),torsion_dgr(i,j1min,j2min,j3min));
				s1+="\n";
				writer.append(s1);
			}
		}
	}

    protected void WriteCAR(Writer writer) throws  IOException{
        writer.append("Not implemented yet.\n");
    }

	protected void WriteEPS(Writer writer) throws IOException{
		if(nAtoms==0) return;

 		double zoom=100; //determine how many point size 1 Angtrom should be converted
		double ro=0.3; //radius of atom measured in Angstrom, will be scaled according to atomic mass, this value is for oxygen.
		double xo=100,yo=100; //coordinates of origin
		double width=300,height=300; //size of screen
		double angle=Math.toRadians(90.0); //perspective angle

		String s="%!PS-Adobe-2.0 EPSF-2.0 \n"
				+"%%BoundingBox:  " + String.format("%d %d %d %d \n",(int)xo,(int)yo,(int)(xo+width),(int)(yo+height))
				+"%%DocumentFonts: Helvetica-Bold \n"
				+"%%Title: Snapshot postscript file from Moviemol \n"
				+"%%Pages: 1 \n"
				+"%%EndComments \n%%EndProlog \n%%Page: 1 1 \n";
		//define macro
		s+="/mm { 2.83465 mul } def  \n";
		s+="/doAAtom \n{ newpath \n"
				+"0 360 arc \n"
				+"gsave \nsethsbcolor \nfill \ngrestore \n.0 setgray \nsetlinewidth \n"
				+"stroke } def \n";
		s+="/doALamp { newpath    0 360 arc    1.0 setgray    fill } def\n";
		s+="/bondArcn { newpath  arc } def\n";
		s+="/bondFill  { lineto      lineto      closepath      sethsbcolor    fill } def \n";
		s+="/bondFillSimple  { newpath moveto lineto lineto lineto closepath sethsbcolor    fill } def\n";
		s+="/bondOtl { newpath    moveto    lineto   moveto   lineto } def \n";
		s+="/bondWidt  {  .0 setgray     setlinewidth    stroke } def\n";
		s+="/showSymbol  { /Helvetica-Bold findfont   exch scalefont setfont   moveto  .0 setgray  show } def\n";

		s+="595 0 translate\n";
		 s+="90 rotate\n";
		 s+="1 setlinejoin\n";
		 s+="420 297 translate\n";
		 s+="-90 rotate\n";
		 s+="0.71 0.71 scale\n";
		 s+="-420 -297 translate\n";

		 s+="/black       {  0.00 0.00 0.00 } def\n";
		 s+="/darkblue    {  0.66 1.00 0.50 } def\n";
		 s+="/green       {  0.33 1.00 1.00 } def\n";
		 s+="/lightblue   {  0.66 0.35 1.00 } def\n";
		 s+="/red         {  0.00 1.00 1.00 } def\n";
		 s+="/purple      {  0.78 1.00 0.90 } def\n";
		 s+="/toffee      {  0.14 0.80 0.90 } def\n";
		 s+="/lightgray   {  0.00 0.00 0.90 } def\n";
		 s+="/gray        {  0.00 0.00 0.70 } def\n";
		 s+="/blue        {  0.66 1.00 1.00 } def\n";
		 s+="/greenyellow {  0.20 1.00 1.00 } def\n";
		 s+="/turquoise   {  0.50 1.00 1.00 } def\n";
		 s+="/salmonpink  {  0.00 0.18 1.00 } def\n";
		 s+="/magenta     {  0.83 1.00 1.00 } def\n";
		 s+="/yellow      {  0.16 1.00 1.00 } def\n";
		 s+="/white       {  0.00 0.00 1.00 } def\n";
		 s+="/lightyellow {  0.16 0.30 1.00 } def\n";
		 s+="/pink        {  0.00 0.36 1.00 } def\n";
		 s+="/steelgray   {  0.66 0.10 0.90 } def\n";
		 s+="/darkgray    {  0.00 0.00 0.50 } def\n";
		 s+="/bluegreen   {  0.50 0.00 0.70 } def\n";
	//                         H       He       Li         B                 C         F      O
	String[] AtomColor={"red","yellow","yellow","yellow","yellow","yellow","yellow","yellow","red"};

		Center();

	double D=width*0.5/Math.tan(angle/2);
	double maxX=0,maxY=0;

		for(int i=0;i<nAtoms;i++){
			maxX=Math.max(Math.abs(coords[i*3]),maxX);
			maxY=Math.max(Math.abs(coords[i*3+1]),maxY);
		}
		maxX+=0.5;
		maxY+=0.5;
		zoom=Math.min(width/(2.0*maxX),height/(2.0*maxX));

		for(int i=0;i<nAtoms;i++){
			double[] v={coords[i*3],coords[i*3+1],coords[i*3+2]};

			double xpos = xo + v[0] * zoom*(v[2]+D)/D + width * 0.5;
			double ypos = yo + v[1] * zoom*(v[2]+D)/D + height * 0.5;
			double rpos	=  	ro *Math.pow(getMass(i)/16.0,0.15)* zoom*(v[2]+D)/D;

			s+=String.format("1.05 %s  %1.1f  %1.1f %1.1f  doAAtom \n",AtomColor[Nz[i]],xpos,ypos,rpos);
		}

		s+=" showpage \n";

		//writing extra information
		//output
		writer.append(s);
	}

    protected void WritePDB(Writer writer) throws IOException{
        if(nAtoms>=1){
            String s1 = "TITLE     "+tag+ "\n";
			//writer.append(s1);
           // s1="AUTHOR nqcmol\n";
            //s1="  P 1           1\n";   writer.append(s1);
           // s1+="MODEL        1\n"; writer.append(s1);
            int nMol=0;
            for(int i=0;i<nAtoms;i++){
                s1=String.format("ATOM  %5d",i+1);
                if(IsNonHydrogen(i)) nMol++;
                s1+=String.format(" %4s SOL  %4d    %8.3f%8.3f%8.3f  1.00  0.00\n",cElements[Nz[i]],nMol,coords[i*3],coords[i*3+1],coords[i*3+2]);
                writer.append(s1);
            }
            writer.append("END\n");
        }
    }

}
