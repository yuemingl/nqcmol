/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.cluster;

import java.util.Vector;

/**
 *
 * @author chinh
 */
public class VibrationData implements Cloneable{

    public VibrationData(int nAtoms) {
        this.setNAtoms(nAtoms);
    }

    public VibrationData(VibrationData aThis) {
        set(aThis);
    }  
           
    @Override
	public Object clone() {
           return new VibrationData(this);
        }

    public void initialize(int nAtoms){
        this.setNAtoms(nAtoms);
    }

     public void set(VibrationData src){
		setNAtoms(src.getNAtoms());
		setFreqs(src.getFreqs());
		setForceConst(src.getForceConst());
        setReducedMass(src.getReducedMass());
        setNormalModeVectors(src.getNormalModeVectors());
	}

    protected int nModes = 0;

    /**
     * Get the value of nModes
     *
     * @return the value of nModes
     */
    public int getNModes() {
        return nModes;
    }

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
        if(nAtoms<=1) nModes=0;
        else if(nAtoms==2) nModes=1;
        else nModes=nAtoms*3-6;

        if(nModes>0){
            freqs=new double[nModes];
            kConst=new double[nModes];
            reducedMass=new double[nModes];
            normalModeVectors=new double[nModes][nAtoms*3];
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
		this.freqs = (freqs!=null)?freqs.clone():null;
	}

    /**
	 * Set the value of freqs
	 *
	 * @param freqs new value of freqs
	 */
	public void setFreqs(Vector<Double> freqs) {
        this.freqs=MolExtra.ConvertToDoubleArray(freqs);
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
       this.reducedMass = (reducedMass!=null)?reducedMass.clone():null;
	}

    /**
	 * Set the value of reducedMass
	 *
	 * @param reducedMass new value of reducedMass
	 */
	public void setReducedMass(Vector<Double> reducedMass) {
       this.reducedMass = MolExtra.ConvertToDoubleArray(reducedMass);
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
	 * Set the value of normalModeVectors
	 *
	 * @param normalModeVectors new value of normalModeVectors
	 */
	public void setNormalModeVectors(Vector<double[]> normalModeVector) {
		if(normalModeVector.size()>0){
			this.normalModeVectors=new double[normalModeVector.size()][];
			for(int i=0;i<normalModeVector.size();i++)
				this.normalModeVectors[i]=normalModeVector.get(i).clone();
		}
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
	public void setNormalModeVectors(int index, double[] normalModeVectors) {
		this.normalModeVectors[index] = (normalModeVectors!=null)?normalModeVectors.clone():null;
	}

    protected double[] kConst;

    /**
     * Get the value of kConst
     *
     * @return the value of kConst
     */
    public double[] getForceConst() {
        return kConst;
    }

    /**
     * Set the value of kConst
     *
     * @param kConst new value of kConst
     */
    public void setForceConst(double[] kConst) {
        this.kConst = kConst;
    }

    /**
     * Get the value of kConst at specified index
     *
     * @param index
     * @return the value of kConst at specified index
     */
    public double getForceConst(int index) {
        return this.kConst[index];
    }

    /**
     * Set the value of kConst at specified index.
     *
     * @param index
     * @param newKConst new value of kConst at specified index
     */
    public void setForceConst(int index, double newKConst) {
        this.kConst[index] = newKConst;
    }

    protected double[] IRIntensity = null;

    /**
     * Get the value of IRIntensity
     *
     * @return the value of IRIntensity
     */
    public double[] getIRIntensity() {
        return IRIntensity;
    }

    /**
     * Set the value of IRIntensity
     *
     * @param IRIntensity new value of IRIntensity
     */
    public void setIRIntensity(double[] IRIntensity) {
            this.IRIntensity = (IRIntensity!=null)?IRIntensity.clone():null;
    }

    /**
     * Set the value of IRIntensity
     *
     * @param IRIntensity new value of IRIntensity
     */
    public void setIRIntensity(Vector<Double> IRIntensity) {
        this.IRIntensity = MolExtra.ConvertToDoubleArray(IRIntensity);
    }

    /**
     * Get the value of IRIntensity at specified index
     *
     * @param index
     * @return the value of IRIntensity at specified index
     */
    public double getIRIntensity(int index) {
        return this.IRIntensity[index];
    }

    /**
     * Set the value of IRIntensity at specified index.
     *
     * @param index
     * @param newIRIntensity new value of IRIntensity at specified index
     */
    public void setIRIntensity(int index, double newIRIntensity) {
        this.IRIntensity[index] = newIRIntensity;
    }

}
