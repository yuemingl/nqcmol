/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.logging.*;

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
			TotalMass+=Mass(i);
		return TotalMass;
	}

	/**
	 * @param i index of atom
	 * @return Mass of an atom
	 */
	public double Mass(int i){
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

        s1=String.format("%d %1.10f %f ",tag,energy,rmsGrad) +"% @IN\n";
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
