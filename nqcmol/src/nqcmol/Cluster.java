/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.StringTokenizer;
import java.util.logging.*;

/**
 *
 * @author nqc
 */
public class Cluster {
	public Cluster(){
	};
	
	protected static Logger logger=Logger.getLogger(Cluster.class.getName());

	/* Maximum number of atomic number
	 */
	static final int cTypeMax=15;

	/* Name of chemica elements
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
		this.nAtoms = nAtoms;
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
		this.coords = coords;
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


	protected double rmsGrad;

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
	 *
	 * @param index of atom
	 * @return the value of Nz
	 */
	public int getAtomicNumber(int index) {
		return Nz[index];
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

	public boolean Read(InputStream is,String format){
		boolean isReadable=false;
		if(format.contentEquals("xyz"))		isReadable=ReadXYZ(is);
		//if(format.contentEquals("xyz"))		ReadXYZ(is);
//			case F_GRO: ReadGro(is);	break; //for GROMACS format
//			case F_PDB: ReadPDB(is);	break; //for protein data bank format
//			case F_GAUOUT: ReadGauOut(is);	break; //for Gaussian output format
//			case F_VASPOUT: ReadVASPOut(is); break; //for VASP output format
//		}
		return isReadable;
	}

	protected boolean ReadXYZ(InputStream is){
		try {
			Clear();
			
			BufferedReader input = new BufferedReader(new InputStreamReader(is));
			
			StringTokenizer tokenizer;
			String line = input.readLine();
			
			if (input.ready() && line != null) {
				// parse frame by frame
				tokenizer = new StringTokenizer(line, "\t ,;"); //read no of atoms
				String info = tokenizer.nextToken();	nAtoms = Integer.parseInt(info);
				
				//System.out.printf(" Here %s nAtoms = %d \n",info,nAtoms);

				if (nAtoms <= 0) {	return false; }
				line = input.readLine();		tokenizer = new StringTokenizer(line, "\t ,;"); //read tag, energy and so on

				info = tokenizer.nextToken();	tag = Integer.parseInt(info);

				if (line.contains("@IN")) {
					//for my format
					info = tokenizer.nextToken();		energy = Double.parseDouble(info);
					info = tokenizer.nextToken();		rmsGrad = Double.parseDouble(info);
				}
				//System.out.printf(" Here tag=%d Energy = %f rmsGrad=%f\n",tag,energy,rmsGrad);
				
				coords = new double[nAtoms * 3]; //start reading coordinates
				Nz=new int [nAtoms];

				//System.out.printf(" Here tag=%d Energy = %f rmsGrad=%f\n",Nz[0],energy,rmsGrad);

				//grad=new double[nAtoms*3];
				for (int i = 0; i < nAtoms; i++) {
					line = input.readLine();
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

						}
					}
				}	
				
			}
		} catch (IOException ex) {
			Logger.getLogger(Cluster.class.getName()).log(Level.SEVERE, null, ex);
		}

		//System.out.printf(" Here \n");
		return true;
	}

	public void Write(Writer os,String format) throws IOException {
		if(nAtoms>0){
			if(format.contentEquals("xyz")) WriteXYZ(os);
		}
	}

	public void Write(OutputStream os,String format) throws IOException {
		Write( new BufferedWriter(new OutputStreamWriter(os)),format);
		if(nAtoms>0){
			if(format.contentEquals("xyz")) WriteXYZ(os);
		}
	}

	protected void WriteXYZ(Writer os) throws  IOException{		
        String s1 = "" + nAtoms + "\n";		
        os.append(s1);
		//System.out.printf(" Here %s \n",s1);

        s1=String.format("%d %1.10f %f ",tag,energy,rmsGrad) +"% @IN\n";
		os.append(s1);
		//System.out.printf(" Here %s \n",s1);

        // Loop through the atoms and write them out:
        for(int i=0;i<nAtoms;i++){
			s1=cElements[Nz[i]]+String.format(" %1.6f %1.6f %1.6f\n",coords[i*3],coords[i*3+1],coords[i*3+2]);
			os.append(s1);
			System.out.printf(" Here %s \n",s1);
        }
	}
}