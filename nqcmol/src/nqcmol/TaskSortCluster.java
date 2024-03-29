/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
import java.io.*;
import java.util.ArrayList;
import java.util.Collections;
import java.util.Comparator;
import java.util.Scanner;
import nqcmol.tools.MTools;
import org.kohsuke.args4j.*;

/**
 *
 * @author nqc
 */
public class TaskSortCluster extends Task{

	@Option(name = "-r", usage = "sort with decreasing of energies. Default is increasing.")
	boolean bReverse=false;

	@Override
	public String getName(){
		return "SortCluster";
	}

    static final public String Option="sort";

    static final public String Descriptions="\t "+Option+" \t - "+ "sort clusters according to energies\n";

    public static void main(String[] args) throws IOException  {
        new TaskSortCluster().Execute(args);
	}

	@Override
	protected void Initialize() {
		super.Initialize();		
		//xmllog.writeAttribute("USRthreshold", Double.toString(lim));
		//xmllog.writeAttribute("deltaE", Double.toString(deltaE));
	}

	@Override
	protected void Process() {
		try {
			xmllog.writeEntity("Note");
			xmllog.writeText("Sort clusters according to their energies\n");
			xmllog.endEntity();

			Scanner fileIn = new Scanner(new File(sFileIn));

			int i = 0;
			ArrayList<Cluster> pop=new ArrayList<Cluster>();
			while (fileIn.hasNext()) {
				Cluster mol = new Cluster();
				mol.Read(fileIn, sFormatIn);
				
				if(mol.getNAtoms()<=0) break;
				pop.add(mol);
				//if(i>2)pop.get(i-1).Write(System.err,MTools.getExtension(sFileOut));
				i++;
				
			}
            fileIn.close();
			xmllog.writeEntity("NumOfInputClusters").writeText(Integer.toString(pop.size())).endEntity();

			Collections.sort( pop, new Comparator()
                     {
					 @Override
                     public int compare( Object a, Object b ){
							int result=(int)(Math.signum( ((Cluster)a).getEnergy() - ((Cluster)b).getEnergy() ));
							return result;
                        }
                     } );
					 
			if(bReverse) Collections.reverse(pop);



			if(!sFileOut.isEmpty()){
				FileWriter fileOut=new FileWriter(new File(sFileOut));
			
				for(i=0;i<pop.size();i++){
					pop.get(i).Write(fileOut,MTools.getExtension(sFileOut));
				}
			}
			

		} catch (IOException ex) {
			//Logger.getLogger(TaskSingleCluster.class.getName()).log(Level.SEVERE, null, ex);
		}
	}


}
