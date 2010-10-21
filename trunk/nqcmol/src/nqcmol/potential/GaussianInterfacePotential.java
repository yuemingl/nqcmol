/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol.potential;

import java.io.*;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.cluster.Cluster;
import nqcmol.cluster.MolExtra;
import nqcmol.tools.XmlWriter;

/**
 *
 * @author nqc
 */
public class GaussianInterfacePotential extends Potential{
	public GaussianInterfacePotential(){
		nativeUnit="Hartree";
		unit=nativeUnit;
	}

	@Override
	public String XMLInfo(int verbose) {
		String info = "";
		try {
			//!< print setting parameters
			Writer writer = new java.io.StringWriter();
			XmlWriter xmlwriter = new XmlWriter(writer);
			xmlwriter.writeEntity("Potential");
			if ((verbose & 1) !=0) {
				xmlwriter.writeAttribute("Equation", getEquation()).writeAttribute("Unit",getUnit());
				xmlwriter.writeAttribute("BasisSet", basisSet);
				xmlwriter.writeAttribute("ChargeAndMultiplicity", chargeAndMultiplicity);
			}
			xmlwriter.endEntity();
			xmlwriter.close();
			info = writer.toString();

		} catch (IOException ex) {
			Logger.getLogger(Potential.class.getName()).log(Level.SEVERE, null, ex);
		}
		return info;
	}



	@Override
	public String getEquation(){
		String equation="Gaussian";
		return equation;
	}
	
	@Override
	public boolean HasAnalyticalGradients() {
		return false;
	}

	@Override
	public boolean isValidSetup() {
		return true;
	}

	
	@Override
	public void setCluster(Cluster cluster_) {
		super.setCluster(cluster_);
		//AllocatePrivateVariables(cluster_.getNAtoms());
		cluster.CorrectOrder();
		//CalcOH();
		int nO=cluster.getNonHydrogenNum();
		for(int i=0;i<nO;i++){
			int count=0;
			for(int j=nO;j<cluster.getNAtoms();j++)
				if((cluster.distance(i,j)<1.3)&&(nO+2*i+count<cluster.getNAtoms()-1)){
					cluster.SwapAtom(nO+2*i+count,j);
					count++;
				}
		}
	}

	protected String chargeAndMultiplicity = "0 1";

    /**
     * Get the value of chargeAndMultiplicity
     *
     * @return the value of chargeAndMultiplicity
     */
    public String getChargeAndMultiplicity() {
        return chargeAndMultiplicity;
    }

    /**
     * Set the value of chargeAndMultiplicity
     *
     * @param chargeAndMultiplicity new value of chargeAndMultiplicity
     */
    public void setChargeAndMultiplicity(String chargeAndMultiplicity) {
        this.chargeAndMultiplicity = chargeAndMultiplicity;
    }

	protected String basisSet="b3lyp/6-31+G";

	/**
	 * Get the value of basisSet
	 *
	 * @return the value of basisSet
	 */
	public String getBasisSet() {
		return basisSet;
	}

	/**
	 * Set the value of basisSet
	 *
	 * @param basisSet new value of basisSet
	 */
	public void setBasisSet(String basisSet) {
		this.basisSet = basisSet;
	}

	@Override
	protected double Energy_(double[] _p){
		double energy = 0;
		if (_p.length <= 3) return 0;
		//System.out.print(" Size of vec = "+	Integer.toString(_p.length));
       
       int nAtoms=_p.length/3;
       String input="";

       input+="%mem=1GB \n%nproc=1 \n";
        input+="#p "+basisSet+"\n";
        input+="\nSingle energy calculation \n\n";

        input+=chargeAndMultiplicity+"\n";//Charge and multiplicity

        for(int i=0;i<cluster.getNAtoms();i++){
            input+=cluster.getAtomicSymbol(i)+" ";
            if(i!=nAtoms-1)
                input+=String.format("%12.8f %12.8f %12.8f\n",_p[i*3+0],_p[i*3+1],_p[i*3+2]);
            else
                input+=String.format("%12.8f %12.8f %12.8f\n\n ",_p[i*3],_p[i*3+1],_p[i*3+2]);
        }

        //System.err.println(input);
        double[] energies=getEnergies(input);
        if(energies!=null) energy=energies[energies.length-1];
        else energy=0;


		energy=ConvertUnit(energy,nativeUnit,unit);

		nEvals++;
		return energy;
	}

	@Override
	protected void Gradient_(double[] _p, double[] gradient) {
		double energy = 0;
		if (_p.length <= 3) return ;
		//System.out.print(" Size of vec = "+	Integer.toString(_p.length));
       try {
		   int nAtoms=_p.length/3;
		   int nO=nAtoms/3;
		   String s=String.format("echo -e \"%d\n\n",nAtoms);
		   for(int i=0;i<nAtoms;i++){
			   if(i<nO) s+="O ";
			   else s+="H ";
			   if(i!=nAtoms-1)
					s+=String.format("%1.15f %1.15f %1.15f\n",_p[i*3],_p[i*3+1],_p[i*3+2]);
			   else
				   s+=String.format("%1.15f %1.15f %1.15f\" | $ttm21f",_p[i*3],_p[i*3+1],_p[i*3+2]);
		   }

			//input="./ttm21f < input.xyz";
			//System.err.println(input);
			Runtime r = Runtime.getRuntime();
			Process p = r.exec(new String[] {"/bin/bash","-c", s});
			InputStream in = p.getInputStream();
			BufferedInputStream buf = new BufferedInputStream(in);
			InputStreamReader inread = new InputStreamReader(buf);
			BufferedReader bufferedreader = new BufferedReader(inread);
			// Read the ls output
			String line= bufferedreader.readLine();
			if(!line.contains("Error")){
				StringTokenizer tokenizer;
				tokenizer = new StringTokenizer(line, "\t ,;"); //read no of atoms
				//System.err.println(line);
				energy=Double.parseDouble(tokenizer.nextToken());
			}
			// Check for ls failure
			try {
				if (p.waitFor() != 0) {
					System.err.println("exit value = " + p.exitValue());
				}
			} catch (InterruptedException e) {
				System.err.println(e);
			} finally {
				// Close the InputStream
				bufferedreader.close();
				inread.close();
				buf.close();
				in.close();
			}
          } catch (IOException ex) {
                  Logger.getLogger(GaussianInterfacePotential.class.getName()).log(Level.SEVERE, null,ex);
          }
		nEvals+=3;
	}


	/**
	 * @param GaussianInput String of Gaussian input
	 * @return array of energy collected
	 */
	public double[] getEnergies(String sGaussianInput){
		double[] energies= null;
		Vector<Double> vecEnergyMP2=new Vector<Double>();
		Vector<Double> vecEnergyDFT=new Vector<Double>();
		try {
			String command = "echo -e \"" + sGaussianInput + "\" | g03 ";
			//System.err.print(command);
			Runtime r = Runtime.getRuntime();
			Process p = r.exec(new String[]{"/bin/bash", "-c", command});
			// Query the results of the process
			BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
			String s;
			while ((s = in.readLine()) != null) {
				//System.err.println("Test = "+s);
				double energy=0;
                if(s.contains("E(")&&s.contains("SCF Done")){
                    int pos1=s.indexOf("=")+1;
                    int pos2=s.indexOf("A.U.");
                    String substr=s.substring(pos1,pos2);
                    energy=Double.parseDouble(substr);
                    vecEnergyDFT.add(energy);
                }

                if ( s.contains("EUMP2 =")  ){	// wow, energy
                    String substr=s.substring(s.indexOf("EUMP2 =")+7);
                    substr=substr.replace('D','E');
                    //System.out.println(substr);
                    energy=Double.parseDouble(substr);
                    vecEnergyMP2.add(energy);
                }
			}

			try {
				if (p.waitFor() != 0) {
					System.err.println("Error termination with Gaussian. Exit value = " + p.exitValue());
				}
			} catch (InterruptedException e) {
				System.err.println(e);
			} finally {
				// Close the InputStream
				in.close();
			}

		} catch (IOException ex) {
			Logger.getLogger(GaussianInterfacePotential.class.getName()).log(Level.SEVERE, null, ex);
		}

		if(vecEnergyMP2.size()>0){
			energies=MolExtra.ConvertToDoubleArray(vecEnergyMP2);
		}else{
			if(vecEnergyDFT.size()>0){
				energies=MolExtra.ConvertToDoubleArray(vecEnergyDFT);
			}
		}

		return energies;
	}

    @Override
    public boolean Optimize() {
        try {
            if (cluster.getNAtoms() <= 1) {
                return true;
                //System.out.print(" Size of vec = "+	Integer.toString(_p.length));
            }
            String sGaussianInput = "";
            sGaussianInput += "%mem=1GB \n%nproc=1 \n";
            sGaussianInput += "#opt " + basisSet + "\n";
            sGaussianInput += "\nSingle energy calculation \n\n";
            sGaussianInput += chargeAndMultiplicity + "\n"; //Charge and multiplicity
            double[] _p = cluster.getCoords();
            for (int i = 0; i < cluster.getNAtoms(); i++) {
                sGaussianInput += cluster.getAtomicSymbol(i) + " ";
                sGaussianInput += String.format("%12.8f %12.8f %12.8f\n", _p[i * 3 + 0], _p[i * 3 + 1], _p[i * 3 + 2]);
                if (i == cluster.getNAtoms() - 1) {
                    sGaussianInput += "\n";
                }
            }
            String command = "echo -e \"" + sGaussianInput + "\" | g03 ";
            //System.err.print(command);
            Runtime r = Runtime.getRuntime();
            Process p = r.exec(new String[]{"/bin/bash", "-c", command});
            // Query the results of the process
            //BufferedReader in = new BufferedReader(new InputStreamReader(p.getInputStream()));
            cluster.Read(new Scanner(p.getInputStream()),"g03");
        } catch (IOException ex) {
            Logger.getLogger(GaussianInterfacePotential.class.getName()).log(Level.SEVERE, null, ex);
        }
        return true;
    }



}
