/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.potential;

/**
 *
 * @author nqc
 */
public class LennardJonesPotential extends Potential{
	public LennardJonesPotential(){};

	@Override
	public String getEquation(){
		String equation="LennardJones";
		return equation;
	}

	


	@Override
	public boolean isValidSetup() {
		return true;
	}

	@Override
	public boolean HasAnalyticalGradients() {
		return true;
	}


	@Override
	protected double Energy_(double[] _p){
		double energy=0;

		//System.out.print(" Size of vec = "+ Integer.toString(coords3d.getSize()));

		int nAtom=_p.length/3;
		for(int i=0;i< nAtom;i++){
			//		rc[0]+=_p[i*3+0];
			//		rc[1]+=_p[i*3+1];
			//		rc[2]+=_p[i*3+2];


			for(int j=i+1;j<nAtom;j++){
				double rij2=  Math.pow( (_p[i*3+0]-_p[j*3+0] ),2)
							+ Math.pow( (_p[i*3+1]-_p[j*3+1] ),2)
							+ Math.pow( (_p[i*3+2]-_p[j*3+2] ),2);
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

	@Override
	protected void Gradient_(double[] _p, double[] grad) {
		int nAtom=_p.length/3;
		for (int i=0; i < nAtom ; i++) {
//			energyGradient.setElement(i,
			for(int k=0;k<3;k++)	grad[i*3+k]=0;
			for(int j=0;j<nAtom;j++)	if(i!=j){
				double rij2=  Math.pow( (_p[i*3+0]-_p[j*3+0] ),2)
							+ Math.pow( (_p[i*3+1]-_p[j*3+1] ),2)
							+ Math.pow( (_p[i*3+2]-_p[j*3+2] ),2);
				double rij6 = Math.pow(rij2,3);

				double rij14= rij6*rij6*rij2;
				double f=24.0*(rij6-2.0)/(rij14);
				for(int k=0;k<3;k++)	grad[i*3+k]+=f*(_p[i*3+k]-_p[j*3+k]);
			}

		}
	}
}
