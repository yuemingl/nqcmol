/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.Collections;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.logging.*;
import nqcmol.tools.MTools;
import org.apache.commons.math.complex.Complex;
import org.apache.commons.math.transform.FastFourierTransformer;
import org.kohsuke.args4j.*;

/**
 *
 * @author nqc
 */
public class TaskAutoCorrelationAndTransform extends Task{
	final int cMaxSteps=50000;
	@Option(name="-from",usage="Index of step where to start. Count from 0. Defaul is 0.",metaVar="INTEGER")
    int iFrom=0;

	@Option(name="-to",usage="Index of step where to stop. Default is running up to maximum possible step",metaVar="DOUBLE")
	int iTo=-1;


	double timestep=0.5e-15;

	@Override
	public String getName(){
		return "AutoCorrelationAndTransform";
	}

	@Override
	protected void Initialize() {
		super.Initialize();
		try {			
			xmllog.writeEntity("Note");
			xmllog.writeText("Calculate auto-correlation and Fourier Transform");
			xmllog.endEntity();
			
		} catch (IOException ex) {
			Logger.getLogger(TaskRemoveDuplicateCluster.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	@Override
	protected void Process() {
		try {			
			if (!sFileOut.isEmpty()) {
				fileOut = new FileWriter(new File(sFileOut));
			}

			xmllog.writeEntity("ReadingInputFiles");
			ReadVelocityFromMDPWCF();
			xmllog.writeAttribute("iFrom", Integer.toString(iFrom));
			xmllog.writeAttribute("iTo", Integer.toString(iTo));
			xmllog.writeAttribute("MaxSteps", Integer.toString(nT));
			xmllog.writeAttribute("NumSteps", Integer.toString(T));
			xmllog.writeAttribute("TimeStep", Double.toString(timestep));
			xmllog.endEntity().flush();

			//Vector<double> s=new Vector<double>();

					
				//System.out.print(MTools.getExtension(sFileOut));
//				xmllog.writeEntity("Cluster");
//					xmllog.writeAttribute("id", Integer.toString(i));
//					xmllog.writeAttribute("Tag", mol.getTag());
//
//				xmllog.endEntity().flush();
					
		} catch (IOException ex) {
			Logger.getLogger(TaskAutoCorrelationAndTransform.class.getName()).log(Level.SEVERE, null, ex);
		}		
	}

	void ReadVelocityFromMDPWCF() throws FileNotFoundException{
			Scanner scanner = new Scanner(new File(sFileIn));
			String s=scanner.nextLine();

			StringTokenizer token = new StringTokenizer(s, "\t ,;"); //read no of atoms
			token.nextToken(); nT=Integer.parseInt(token.nextToken());


			int c=0;
			a0=0;
			boolean bFirstRead=true;
			double[][][] velRead=null; //!< velocity
    		while(scanner.hasNext()){
				s=scanner.nextLine();
				if(s.contains("PRIMVEC 1")&&(bFirstRead)){
					s=scanner.nextLine();	token = new StringTokenizer(s, "\t ,;"); //read no of atoms
					a0=Double.parseDouble(token.nextToken());
					bFirstRead=false;
				}

				if(s.contains("PRIMCOORD")){
					nAtom=0;
					s=scanner.nextLine();	token = new StringTokenizer(s, "\t ,;"); //read no of atoms
				    nAtom=Integer.parseInt(token.nextToken());
					if(velRead==null){
						velRead=new double[cMaxSteps][nAtom][3];
						R0=new double[nAtom][3];
					}

					  for(int i=0;i<nAtom;i++){
					  	s=scanner.nextLine();	token = new StringTokenizer(s, "\t ,;"); //read no of atoms
						token.nextToken();

						if(c==0){
							R0[i][0]=Double.parseDouble(token.nextToken());
							R0[i][1]=Double.parseDouble(token.nextToken());
							R0[i][2]=Double.parseDouble(token.nextToken());
						}else{
							token.nextToken();
							token.nextToken();
							token.nextToken();
						}
						velRead[c][i][0]=Double.parseDouble(token.nextToken());
						velRead[c][i][1]=Double.parseDouble(token.nextToken());
						velRead[c][i][2]=Double.parseDouble(token.nextToken());
					  }
					  c++;
					  if(c>=velRead.length) break;
				  }
    		}

			if(iTo==-1) iTo=velRead.length-1;
			T=iTo-iFrom+1;
			c=0;
			vel=new double[T][nAtom][3];
			for(int k=iFrom;k<=iTo;k++)
				for(int i=0;i<nAtom;i++){
						vel[c][i][0]=velRead[k][i][0];
						vel[c][i][1]=velRead[k][i][1];
						vel[c][i][2]=velRead[k][i][2];
				}
	  }

	void CalculateAutoCorrelation(){
		double[] kVector=new double[nAtom];
		for(int k=0;k<nAtom;k++){
			kVector[k]=Math.cos((K[0]*R0[k][0] + K[1]*R0[k][1] + K[2]*R0[k][2])*2.0*Math.PI/a0);
		}
		
		double s0=0;
	  	for(int t=0;t<T;t++){
			double s=0;
			int count=0;
			// 			for(int k=0;k<nAtom;k++){
			// 				double Ckt=CorrelAtT(k,t);
			// 				corr[k].push_back(Ckt);
			// 				s+=Ckt*kVector[k]/corr[k][0];
			// 				count++;
			// 			}

				for(int k=nAtom-16-1+0;k<nAtom-16-1+2;k+=2){ //for A
				//for(int k=nAtom-16-1+2;k<nAtom-16-1+8;k+=2){ //for B
				//for(int k=nAtom-16-1+8;k<nAtom-16-1+16;k+=2){ //for C
					double Ckt=0;

						for(int i=0;i<T-t;i++)
							Ckt+=Math.pow(vel[i][k][0]-vel[i+t][k+1][0],2)
								+Math.pow(vel[i][k][1]-vel[i+t][k+1][1],2)
								+Math.pow(vel[i][k][2]-vel[i+t][k+1][2],2);
						Ckt/=(T-t);

							//corr[k].push_back(Ckt);
							//s+=Ckt*kVector[k]/corr[k][0];
							count++;
					}

			signal[t]=s/count;


		}
	}

	double[][][] vel=null; //!< actual velocity taken account in
	double[][] R0=null;
	


	double[] K={1,1,1};//k vector
	int nT,T;//number of steps
	int nAtom;
	double a0;
	double[] signal;
	double Fs;//!< sampling frequency


	Complex[] FourierTransform(double[] signal){
		FastFourierTransformer fft=new FastFourierTransformer();
		Complex[] intens=fft.transform(signal);
		return intens;
	}
}
