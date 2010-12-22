/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashSet;
import java.util.Scanner;
import java.util.logging.*;
import org.kohsuke.args4j.*;

/**
 *
 * @author nqc
 */
public class TaskDetectEquivalentAtoms extends Task{
	@Option(name="-tol",usage="distance tolerance in detection",metaVar="DOUBLE")
    double tol=0.01;
	int decimalToRoundOff=2;

	@Override
	public String getName(){
		return "DetectEquivalentAtoms";
	}

    static final public String Option="equiAtom";

    static final public String Descriptions="\t "+Option+" \t - "+ "Detect topologically equivalent atoms in clusters\n";

    public static void main(String[] args) throws IOException  {
        new TaskDetectEquivalentAtoms().Execute(args);
	}

	@Override
	protected void Initialize() {
		super.Initialize();		
        decimalToRoundOff=(int) Math.round(Math.abs(Math.log10(tol)));
        xmllog.writeAttribute("Tolerance",Double.toString(tol));
        xmllog.writeAttribute("DecimalToRoundOff",Integer.toString(decimalToRoundOff));
        xmllog.writeEntity("Note");
        xmllog.writeText("Detect topologically equivalent atoms in clusters using the algorithm described in J. Chem. Inf. Comput. Sci., Vol 37, No 1, 1997.\n");
        xmllog.endEntity();
	}

	@Override
	protected void Process() {
		try {
            FileWriter fileOut=null;
			if (!sFileOut.isEmpty()) {
				fileOut = new FileWriter(new File(sFileOut));
			}
			
			Cluster mol = new Cluster();

			Scanner fileIn = new Scanner(new File(sFileIn));
			int i = 0;
			while (mol.Read(fileIn, sFormatIn)) {	
				//System.out.print(MTools.getExtension(sFileOut));
				xmllog.writeEntity("Cluster");
					xmllog.writeAttribute("id", Integer.toString(i));
					xmllog.writeAttribute("Tag", mol.getTag());
					DetectEquivalentAtoms(mol);
				xmllog.endEntity().flush();
				i++;
				//break;
			}
            fileIn.close();
            fileOut.close();
		} catch (IOException ex) {
			Logger.getLogger(TaskDetectEquivalentAtoms.class.getName()).log(Level.SEVERE, null, ex);
		}		
	}

	void DetectEquivalentAtoms(Cluster mol) throws IOException{
		int n=mol.getNAtoms();

		//first, calculate distance matrix
		double[][] D=new double[n][n];
		 ArrayList listDistance = new ArrayList();
//		for(int i=0;i<n;i++){
//			D[i][i]=0;
//			for(int j=i+1;j<n;j++){
//				D[i][j]=D[j][i]=mol.distance(i,j);
//			}
//		}
//		MTools.PrintArray(D);
		for(int i=0;i<n;i++){
			D[i][i]=0;
			for(int j=i+1;j<n;j++){
				String arg="%1."+decimalToRoundOff+"f";
				//System.out.println(arg);
				double d=mol.distance(i,j)*Math.sqrt(mol.getMass(i)*mol.getMass(j));
				D[i][j]=D[j][i]=Double.parseDouble(String.format(arg, d));
				listDistance.add(D[i][j]);
			}
		}
		HashSet hashSet = new HashSet(listDistance); //Create a HashSet which allows no duplicates
		listDistance = new ArrayList(hashSet); //Assign the HashSet to a new ArrayList
		Collections.sort(listDistance); //Ensure correct order, since HashSet doesn't

		//System.out.println("DistanceMatrix");		MTools.PrintArray(D);

		//System.out.println("ListOfDistances");	for (Object item : listDistance)	System.out.println(item);

		//second, calculate distance class number
		int[][] DCN=new int[n][n];
		for(int i=0;i<n;i++){
			DCN[i][i]=0;
			for(int j=i+1;j<n;j++){
				int pos=Collections.binarySearch(listDistance,D[i][j]);
				assert(pos!=-1);
				//System.out.println(arg);
				DCN[i][j]=DCN[j][i]=pos+1;
			}
		}
		//System.out.println("DistanceClassNumber");		MTools.PrintArray(DCN);

		//third, sort DCN
		boolean[] bCheck=new boolean[n];
		for(int i=0;i<n;i++){
			Arrays.sort(D[i]);
			bCheck[i]=false;
		}

		int id=0;
		for(int i=0;i<n;i++)
			if(bCheck[i]==false){
				bCheck[i]=true;
				xmllog.writeEntity("ClassOfAtoms");
				xmllog.writeAttribute("id",Integer.toString(id));
				String members=Integer.toString(i);
				for(int j=i+1;j<n;j++)
					if(Arrays.equals(D[i], D[j])){
						bCheck[j]=true;
						members+=" "+Integer.toString(j);
					}
				xmllog.writeAttribute("Members",members);
				xmllog.endEntity();
				id++;
			}
	}

}
