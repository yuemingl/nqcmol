/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;



import org.kohsuke.args4j.*;

import nqcmol.symmetry.PointGroup;
import nqcmol.symmetry.Symmetry;

/**
 *
 * @author nqc
 */
public class TaskCalculateSymmetryPointGroup extends Task{

    @Option(name="-final",usage="Final criterion for atom equivalence",metaVar="DOUBLE")
    double tolFinal= 5e-2;

    @Option(name="-primary",usage="Initial criterion for atom equivalence",metaVar="DOUBLE")
    double tolInitial= 5e-2;

    @Option(name="-same",usage="Atoms are colliding if distance falls below this value",metaVar="DOUBLE")
    double tolSame= 1e-3;

	@Override
	public String getName(){
		return "CalculateSymmetryPointGroup";
	}

    static final public String Option="sym";

    static final public String Descriptions="\t "+Option+" \t - "+ "Calculate symmetry point group based on brute-force algorithm\n";

	@Override
	protected void Process() {				
        //System.out.println("What the fuck!"+sFileIn);
        int i = 0;
        Cluster mol=new Cluster();
        while (mol.Read(fileIn, sFormatIn)) {
            Symmetry sym=new Symmetry();
            sym.setToleranceFinal(tolFinal);
            sym.setTolerancePrimary(tolInitial);
            sym.setToleranceSame(tolSame);
            Symmetry.CalculateSymmetry(sym, mol,verbose);

            xmllog.writeEntity("Cluster");
            xmllog.writeAttribute("id", Integer.toString(i));
            xmllog.writeAttribute("nAtoms", Integer.toString(mol.getNAtoms()));
            xmllog.writeAttribute("SymmetryCode", sym.getSymmetryCode());
            PointGroup pg=sym.identifyPointGroup();
            if(pg!=null)
                xmllog.writeAttribute("PointGroup", pg.getGroupName());
            else
                xmllog.writeAttribute("PointGroup", "Unknown");

            xmllog.writeAttribute("PointGroupOrder", Integer.toString(sym.getPointGroupOrder()));

            xmllog.endEntity();
            xmllog.flush();
            i++;
            //break;
        }
        fileIn.close();
	}	
}
