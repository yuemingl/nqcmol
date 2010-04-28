/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import org.kohsuke.args4j.*;



/**
 *
 * @author nqc
 */
public class TaskCutFromLattice extends Task {

	@Option(name = "-formula", usage = "Formula of clusters. It MUST be in quote. For example: (H3O)2(H2O)4(NH3)4", metaVar = "STRING")
	String sFormula = "(H2O)4";

	@Option(name = "-num", usage = "Number of clusters", metaVar = "INTEGER")
	int num=1;

	@Option(name = "-rad", usage = "Confined radius", metaVar = "DOUBLE")
	double radius=5;

	@Option(name = "-DMin", usage = "Minimum distance of any neighbour molecules", metaVar = "DOUBLE")
	double DMin=2.5; 

	@Option(name = "-DMax", usage = "Maximum distance of any neighbour molecules", metaVar = "DOUBLE")
	double DMax=4.0;

    //@Option(name = "-cubic", usage = "Confine molecules in a cubic box instead of a sphere as default", metaVar = "BOOLEAN")
	//boolean bBox=false;

    @Option(name = "-nocenter", usage = "Move mass center to origin", metaVar = "BOOLEAN")
	boolean bNoCenter=false;

	@Override
	public String getName(){
		return "CutFromLattice";
	}

	
	@Override
	protected void Process() {			
		try {
			xmllog.writeAttribute("Formula", sFormula).writeAttribute("NumOfClusters", Integer.toString(num));
			xmllog.writeAttribute("ConfinedRadius", Double.toString(radius));
			xmllog.writeAttribute("DMin", Double.toString(DMin));
			xmllog.writeAttribute("DMax", Double.toString(DMax));

            Lattice lat=new Lattice();
            lat.Read(fileIn, sFormatIn);
            lat.Write(System.out, "car");
		
			for(int i=0;i<num;i++){
				xmllog.writeEntity("Cluster").writeAttribute("id", Integer.toString(i));


				
				xmllog.endEntity().flush();
			}
		} catch (IOException ex) {
			Logger.getLogger(TaskCutFromLattice.class.getName()).log(Level.SEVERE, null, ex);
		}
			
	}
	
}

