/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.potential;

import org.apache.commons.math.analysis.MultivariateRealFunction;

/**
 *
 * @author nqc
 */
public class PotentialFitting implements MultivariateRealFunction {

	protected Potential f;

	/**
	 * Get the value of f
	 *
	 * @return the value of f
	 */
	public Potential getF() {
		return f;
	}

	/**
	 * Set the value of f
	 *
	 * @param f new value of f
	 */
	public void setF(Potential f) {
		this.f = f;
	}

	@Override
	public double value(double[] p){
		return f.getEnergy(p);
	}
}
