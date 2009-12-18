/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;



import nqcmol.symmetry.PointGroup;
import nqcmol.symmetry.Symmetry;

/**
 *
 * @author nqc
 */
public class TaskCalculateSymmetryPointGroup extends Task{

	@Override
	public String getName(){
		return "CalculateSymmetryPointGroup";
	}

	@Override
	protected void Process() {
		try {			
			//System.out.println("What the fuck!"+sFileIn);
			int i = 0;
			Cluster mol=new Cluster();
			while (mol.Read(fileIn, sFormatIn)) {
				Symmetry sym=new Symmetry();
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
				break;
			}
			fileIn.close();
		} catch (IOException ex) {
			Logger.getLogger(TaskCalculateSymmetryPointGroup.class.getName()).log(Level.SEVERE, null, ex);
		} 
	}

	
}
