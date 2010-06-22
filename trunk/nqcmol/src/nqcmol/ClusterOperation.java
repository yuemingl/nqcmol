/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;
import java.util.Random;
import nqcmol.tools.MTools;

/**
 *
 * @author chinh
 */
public class ClusterOperation {
    Random ranMov = new Random();

    double pAtomMove=1.0, AtomdR=0.2;
    double pMolMove=1.0, MoldR, MolMovePhi=180;
    double pMolFlip=1.0, MoldFlipPhi=180;
    double pMolRotate=1.0, MolRotatePhi=180;

    /**
     * a random number of atoms are translated
     * @param m
     */
    public void RandomAtomicMove(Cluster m){
        double[] vec=new double[3];
        for(int i=0;i<m.getNAtoms();i++)
				if(ranMov.nextDouble()< pAtomMove){
					MTools.generateRandomVector(vec, ranMov, AtomdR);
					m.Translate_atom(i,vec,1);
				}
    }

    public void RandomMolMove(Cluster m){
        double theta;
        double[] axis=new double[3];
        double[] vec=new double[3];
        double[] point=new double[3];


        m.getPairwiseBond();
        for(int i=0;i<m.getNAtoms();i++) if(m.IsNonHydrogen(i))
            if(ranMov.nextDouble()< pMolMove){
                theta=Math.toRadians(MolMovePhi*(2.0*ranMov.nextDouble()-1.0));
                MTools.generateRandomVector(axis, ranMov, 1.0);
                //m.Write(System.out, "xyz");
                m.R(i,point);
                
                int[] mol=m.getMolecule(i);
                //MTools.PrintArray(axis);

                m.RotateAxisWithPoint_mol(mol, theta, axis, point);
                //m.Write(System.out, "xyz");

                MTools.generateRandomVector(vec, ranMov, MoldR);
                m.Translate_mol(mol, vec,1);
                //m.Write(System.out, "xyz");
                break;
        }
    }

    public void RandomMolFlip(Cluster m){
        double theta;
        double[] axis=new double[3];
        //double[] vec=new double[3];
        double[] point=new double[3];

        m.getPairwiseBond();
        for(int i=0;i<m.getNAtoms();i++) if(m.IsNonHydrogen(i))
            if(ranMov.nextDouble()< pMolFlip){
                theta=Math.toRadians(MoldFlipPhi*(2.0*ranMov.nextDouble()-1.0));
                MTools.generateRandomVector(axis, ranMov, 1.0);
                m.R(i,point);
                int[] mol=m.getMolecule(i);
                m.RotateAxisWithPoint_mol(mol, theta, axis, point);
                //MTools.generateRandomVector(vec, ranMov, MoldR);
               // m.Translate_mol(mol, vec,1);
        }
    }

    public void RandomMolRotate(Cluster m){
        double theta;
        double[] axis=new double[3];
        //double[] vec=new double[3];
        //double[] point=new double[3];

        m.getPairwiseBond();
        for(int i=0;i<m.getNAtoms();i++) if(m.IsNonHydrogen(i))
            if(ranMov.nextDouble()< pMolRotate){
                theta=Math.toRadians(MolRotatePhi*(2.0*ranMov.nextDouble()-1.0));
                MTools.generateRandomVector(axis, ranMov, 0);
                //m.R(i,point);
                int[] mol=m.getMolecule(i);
                m.RotateAxis_mol(mol, theta, axis);
                //MTools.generateRandomVector(vec, ranMov, MoldR);
               // m.Translate_mol(mol, vec,1);
        }
    }

}
