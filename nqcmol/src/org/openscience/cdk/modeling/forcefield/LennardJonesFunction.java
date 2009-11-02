/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package org.openscience.cdk.modeling.forcefield;

import javax.vecmath.GMatrix;
import javax.vecmath.GVector;
import org.openscience.cdk.interfaces.IAtomContainer;


/**
 *
 * @author nqc
 */
public class LennardJonesFunction implements IPotentialFunction {
	String energyFunctionShape = " Lennard-Jones energy ";
	double energy = 0;
	GVector energyGradient = null;
	GVector order2ndErrorApproximateGradient = null;
	GMatrix energyHessian = null;
	double[] forHessian = null;
	GMatrix order2ndErrorApproximateHessian = null;
	double[] forOrder2ndErrorApproximateHessian = null;
	int functionEvaluationNumber = 0;

	//private LoggingTool logger;

	/**
	 *  Constructor 
	 *
	 */
	public LennardJonesFunction()  {
		//logger.debug("LJEnergyFunction Constructor");
		//logger = new LoggingTool(this);
	}


	@Override
	public String toString(){
		return "Lennard-Jones";
	}


	/**
	 *  Evaluate the LJEnergyFunction for a given molecule
	 *
	 *@param  molecule  Current molecule.
	 *@return        LennardJones energy function value.
	 */
	public double energyFunctionOfAMolecule(IAtomContainer molecule) {

		GVector coords3d = ForceFieldTools.getCoordinates3xNVector(molecule);

		energy = energyFunction(coords3d);

		return energy;
	}


	/**
	 *  Evaluate the LennardJones energy function for a given 3xN point
	 *
	 *@param  coords3d  Current molecule 3xN coordinates.
	 *@return        LennardJones energy function value.
	 */
	public double energyFunction(GVector _p) {

		//double rc[3]={0,0,0};

		energy=0;

		//System.out.print(" Size of vec = "+ Integer.toString(coords3d.getSize()));

		int nAtom=_p.getSize()/3;
		for(int i=0;i< nAtom;i++){
			//		rc[0]+=_p[i*3+0];
			//		rc[1]+=_p[i*3+1];
			//		rc[2]+=_p[i*3+2];


			for(int j=i+1;j<nAtom;j++){
				double rij2=  Math.pow( (_p.getElement(i*3+0)-_p.getElement(j*3+0) ),2)
							+ Math.pow( (_p.getElement(i*3+1)-_p.getElement(j*3+1) ),2)
							+ Math.pow( (_p.getElement(i*3+2)-_p.getElement(j*3+2) ),2);
				double rij6 = Math.pow(rij2,3);
				energy+=4.0*(1.0- rij6)/(Math.pow(rij6,2));
			}
		}
			/*
			rc[0]/=nAtom;rc[1]/=nAtom;rc[2]/=nAtom;
			for(i=0;i<nAtom;i++){
		   rij2=   (pow(_p[i*3+0]-rc[0],2)
		   +pow(_p[i*3+1]-rc[1],2)
		   +pow(_p[i*3+2]-rc[2],2));
		   e+=pow(rij2/4,10);
		}
		*/

		//logger.debug("energy = " + bs.LennardJonesSumEB + "  " + ab.LennardJonesSumEA + "  " + sbi.getFunctionLennardJonesSumEBA() + "  " + t.LennardJonesSumET + "  " + vdwi.getFunctionLennardJonesSumEvdW() + "  " + ei.LennardJonesSumEQ);

		return energy;
	}


	/**
	 *  Evaluate the gradient for the LennardJones energy function in a given 3xN point
	 *
	 *@param  coords3d  Current molecule coordinates.
	 */
	public void setEnergyGradient(GVector _p) {
		//setOrder2ndErrorApproximateEnergyGradient(coords3d);

		//logger.debug("coords3d : " + coords3d);
		energyGradient = new GVector(_p.getSize());

		
		//logger.debug("bs.getGradientLennardJonesSumEB() = " + bs.getGradientLennardJonesSumEB());
		//logger.debug("ab.getGradientLennardJonesSumEA() = " + ab.getGradientLennardJonesSumEA());
		int nAtom=_p.getSize()/3;
		for (int i=0; i < nAtom ; i++) {
//			energyGradient.setElement(i,
			double[] a={0,0,0};
			for(int j=0;j<nAtom;j++)	if(i!=j){
				double rij2=  Math.pow( (_p.getElement(i*3+0)-_p.getElement(j*3+0) ),2)
							+ Math.pow( (_p.getElement(i*3+1)-_p.getElement(j*3+1) ),2)
							+ Math.pow( (_p.getElement(i*3+2)-_p.getElement(j*3+2) ),2);
				double rij6 = Math.pow(rij2,3);

				double rij14= rij6*rij6*rij2;
				double f=24.0*(rij6-2.0)/(rij14);
				for(int k=0;k<3;k++)	a[k]+=f*(_p.getElement(i*3+k)-_p.getElement(j*3+k));
			}
			for(int k=0;k<3;k++)	energyGradient.setElement(i*3+k,a[k]);
		}		
	}


	/**
	 *  Get the gradient for the potential energy function in a given point
	 *
	 *@return        Gradient value
	 */
	public GVector getEnergyGradient() {
		return energyGradient;
		//return order2ndErrorApproximateGradient;
	}


	/**
	 *  Evaluate a 2nd order approximation of the gradient
	 *
	 *
	 *@param  coord3d  Current molecule coordinates.
	 */
	public void set2ndOrderErrorApproximateGradient(GVector point) {
		order2ndErrorApproximateGradient = new GVector(point.getSize());
		double sigma = Math.pow(0.000000000000001,0.33);
		GVector xplusSigma = new GVector(point.getSize());
		GVector xminusSigma = new GVector(point.getSize());

		for (int m = 0; m < order2ndErrorApproximateGradient.getSize(); m++) {
			xplusSigma.set(point);
			xplusSigma.setElement(m,point.getElement(m) + sigma);
			xminusSigma.set(point);
			xminusSigma.setElement(m,point.getElement(m) - sigma);
			order2ndErrorApproximateGradient.setElement(m,(energyFunction(xplusSigma) - energyFunction(xminusSigma)) / (2 * sigma));
		}

		//logger.debug("order2ndErrorApproximateGradient : " + order2ndErrorApproximateGradient);
	}


	/**
	 *  Get the 2nd order error approximate gradient of the angle bending term.
	 *
	 *
	 *@return           Angle bending 2nd order error approximate gradient value.
	 */
	public GVector get2ndOrderErrorApproximateGradient() {
		return order2ndErrorApproximateGradient;
	}


	/**
	 *  Evaluate the hessian for the potential energy function in a given point
	 *
	 *@param  point  Current molecule coordinates.
	 */
	public void setEnergyHessian(GVector point) {
		/*double[] forHessian = {2,0,0,0,4,0,0,0,1};
		energyHessian = new GMatrix(3, 3, forHessian);
		*/
		set2ndOrderErrorApproximateHessian(point);
	}


	/**
	 *  Get the hessian for the potential energy function in a given point
	 *
	 *@return        Hessian value
	 */
	public GMatrix getEnergyHessian() {
		//return energyHessian;
		return order2ndErrorApproximateHessian;
	}


	/**
	 *  Get the hessian of the potential energy function in a given point.
	 *
	 *@return        Hessian energy value in the wished point.
	 */
	 public double[] getForEnergyHessian() {
		 //return forHessian;
		 return forOrder2ndErrorApproximateHessian;
	 }


	/**
	 *  Evaluate a 2nd order approximation of the Hessian
	 *  given the atoms coordinates.
	 *
	 *@param  coord3d  Current molecule coordinates.
	 */
	public void set2ndOrderErrorApproximateHessian(GVector coord3d) {
		forOrder2ndErrorApproximateHessian = new double[coord3d.getSize() * coord3d.getSize()];
		double sigma = Math.pow(0.000000000000001,0.33);
		GVector xplusSigma = new GVector(coord3d.getSize());
		GVector xminusSigma = new GVector(coord3d.getSize());
		GVector gradientAtXplusSigma = new GVector(coord3d.getSize());
		GVector gradientAtXminusSigma = new GVector(coord3d.getSize());

		int forHessianIndex;
		for (int i = 0; i < coord3d.getSize(); i++) {
			xplusSigma.set(coord3d);
			xplusSigma.setElement(i,coord3d.getElement(i) + sigma);
			setEnergyGradient(xplusSigma);
			gradientAtXplusSigma.set(this.getEnergyGradient());
			xminusSigma.set(coord3d);
			xminusSigma.setElement(i,coord3d.getElement(i) - sigma);
			setEnergyGradient(xminusSigma);
			gradientAtXminusSigma.set(this.getEnergyGradient());
			for (int j = 0; j < coord3d.getSize(); j++) {
				forHessianIndex = i*coord3d.getSize()+j;
				forOrder2ndErrorApproximateHessian[forHessianIndex] = (gradientAtXplusSigma.getElement(j) - gradientAtXminusSigma.getElement(j)) / (2 * sigma);
				//(energyFunction(xplusSigma) - 2 * fx + energyFunction(xminusSigma)) / Math.pow(sigma,2);
			}
		}
		forOrder2ndErrorApproximateHessian[8] = 1;
		order2ndErrorApproximateHessian = new GMatrix(coord3d.getSize(), coord3d.getSize());
		order2ndErrorApproximateHessian.set(forOrder2ndErrorApproximateHessian);
		//logger.debug("order2ndErrorApproximateHessian : " + order2ndErrorApproximateHessian);
	}


	/**
	 *  Get the 2nd order error approximate Hessian
	 *
	 *
	 *@return           2nd order error approximate Hessian value.
	 */
	public GMatrix get2ndOrderErrorApproximateHessian() {
		return order2ndErrorApproximateHessian;
	}
}
