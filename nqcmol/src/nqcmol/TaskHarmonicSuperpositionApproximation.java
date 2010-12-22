/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
import nqcmol.cluster.MolExtra;
import java.io.*;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.Scanner;
import java.util.logging.*;
import nqcmol.potential.Potential;
import nqcmol.symmetry.PointGroup;
import nqcmol.symmetry.Symmetry;
import nqcmol.tools.MTools;
import org.kohsuke.args4j.*;


/**
 *
 * @author nqc
 */
public class TaskHarmonicSuperpositionApproximation extends Task{	
	@Override
	public String getName(){
		return "HSA";
	}

    static final public String Option="hsa";

    static final public String Descriptions="\t "+Option+" \t - "+ "Harmonic Superposition Approximation\n";

    double Emin=-2;
	double Emax;

	@Option(name="-u",usage="unit of energy (Hartree,kJ/mol,kcal/mol or eV). [Hartree]",metaVar="STRING")
    String sUnit="Hartree";

	//@Option(name="-T",usage="Range of Temperature: Tmin Tmax interval, in Kelvin. For example \"-T 10 100 1\". [10 300 1]",metaVar="DOUBLE DOUBLE DOUBLE")
	double Tmin=10,Tmax=300,dT=1;

	@Option(name = "-regime", usage = "Regime calculation: CLASSICAL, QUANTUM or NOFREQ", metaVar = "STRING")
	REGIME regime=REGIME.CLASSICAL;

	@Option(name = "-Ecutoff", usage = "Energy cutoff. Cluster with the energies larger than Ecutoff will be ignored in calculation.", metaVar = "INTEGER")
	double Ecutoff=0;

    @Option(name="-auto",usage="if specified, automatically classify clusters according to topologies. [false]")
    boolean bAutoClassify=false;

    @Option(name = "-index", usage = "if specified, clusters are classified according to index file. []", metaVar = "STRING")
	String sIndexFile="";

	double beta=315774.659661423;
	double planck=4.55633500389233e-06;

	final double cBeta[]={315774.659661423,503.219018217726,120.272231887602};//for hartree,kcal/mol, kJ/mol
	final double cPlanck[]={4.55633500389233e-06,0.00285914300348446,0.0119626543265790};//for hartree,kcal/mol, kJ/mol
	final double cBetaPlanck=1.43877513515753; //cm*K

    public static void main(String[] args) throws IOException  {
        new TaskHarmonicSuperpositionApproximation().Execute(args);
	}

	@Override
	protected void Initialize() {		
        super.Initialize();
        xmllog.writeAttribute("Unit", sUnit);
        xmllog.writeAttribute("Regime", regime.toString());
        xmllog.writeAttribute("Ecutoff", Double.toString(Ecutoff));
        xmllog.writeEntity("Note");
        xmllog.writeText("Harmonic superposition approximation calculation\n");
        xmllog.endEntity();
        ReadInputs();
	}
		
	public void ReadInputs() {
		try {
            //determine the range of temperature
			int d = MolExtra.getIndexOfArguments("-T", args);
			if ((d != -1) && (d <= args.length - 1 - 3)) {
				Tmin = Double.parseDouble(args[d + 1]);
				Tmax = Double.parseDouble(args[d + 2]);
				dT = Double.parseDouble(args[d + 3]);
			}

			xmllog.writeEntity("Initialize");

            ///reading input clusters and assign groupname for each cluster

             if(bAutoClassify){//read from 1 file, group name of each cluster is determined automatically
                    Scanner scanner = new Scanner(new File(sFileIn));
                    while (scanner.hasNext()) {
                        PartFunc mol = new PartFunc();
                        if (!mol.Read(scanner, sFormatIn)) 	break;
                        String groupname = mol.getMorphology(true);
                        mol.setGroup(groupname);
                        pop.add(mol);
                    }
                    scanner.close();
                    
            }else  if(!sIndexFile.isEmpty()){//read from 1 file, group name of clusters is determined from index file.
                Scanner scanner = new Scanner(new File(sFileIn));
                while (scanner.hasNext()) {
                    PartFunc mol = new PartFunc();
                    if (!mol.Read(scanner, sFormatIn)) 	break;
                    pop.add(mol);
                }
                scanner.close();

                Scanner scannerIndex=new Scanner(new File(sIndexFile));
                String groupname=scannerIndex.next();
                while(scannerIndex.hasNext()){                    
                    String tag=scannerIndex.next();
                    if(tag.contentEquals("-1") && (scannerIndex.hasNext())) groupname=scannerIndex.next();
                    else{
                        for(int i=0;i<pop.size();i++) if(pop.get(i).getTag().contentEquals(tag)){
                            pop.get(i).setGroup(groupname);
                            break;
                        }
                    }

                }
                scannerIndex.close();

            }else { //read from multiple files, each one corresponds to one group
				for (int i = 0; i < 50; i++) {
					String opt = String.format("-i%d", i + 1);
					int j = MolExtra.getIndexOfArguments(opt, args);
					if ((j != -1) && (j <= args.length - 1 - 2)) {
						sFileIn = args[j + 1];
						String groupname = args[j + 2];
                        
						Scanner scanner = new Scanner(new File(sFileIn));
						while (scanner.hasNext()) {
							PartFunc mol = new PartFunc();
							if (!mol.Read(scanner, sFormatIn))	break;
							mol.setGroup(groupname);							
							pop.add(mol);
						}
						scanner.close();
					}
				}
			}

            //report if verbose >1
            if(verbose>=1){
                for(int i=0;i<pop.size();i++){
                    PartFunc mol=pop.get(i);
                    xmllog.writeEntity("Cluster");
                    xmllog.writeAttribute("id", Integer.toString(i));
                    xmllog.writeAttribute("Tag", mol.getTag());
                    xmllog.writeAttribute("Energy", Double.toString(mol.getEnergy()));
                    xmllog.writeAttribute("SymmetryPointGroup", mol.calcPointGroup());
                    xmllog.writeAttribute("Group", mol.getGroup());
                    xmllog.endEntity().flush();
                }
            }
            
			xmllog.endEntity().flush();
		} catch (IOException ex) {
			Logger.getLogger(TaskHarmonicSuperpositionApproximation.class.getName()).log(Level.SEVERE, null, ex);
		}
	}


	@Override
	protected void Process() {
		try {			

			//				if (mol.getEnergy()>Ecutoff) {// if energy is higher than Ecutoff, ignore
			//					xmllog.writeAttribute("Note", "Removed");
			//					xmllog.endEntity().flush();
			//					continue;
			//				}

							//find energy range
			for(int i=0;i<pop.size();i++){
				PartFunc mol=pop.get(i);
				if(i==0){ Emin=Emax=mol.getEnergy();}
				else{
						Emin=Math.min(Emin, mol.getEnergy());
						Emax=Math.max(Emax, mol.getEnergy());
				}
			}

			beta=cBeta[Potential.getIndexOfUnit(sUnit)];
			planck=cPlanck[Potential.getIndexOfUnit(sUnit)];

			xmllog.writeEntity("Calculate");
			xmllog.writeAttribute("Tmin", Double.toString(Tmin));
			xmllog.writeAttribute("Tmax", Double.toString(Tmax));
			xmllog.writeAttribute("deltaT", Double.toString(dT));
			xmllog.writeAttribute("Emin", Double.toString(Emin));
			xmllog.writeAttribute("Emax", Double.toString(Emax));
			xmllog.writeAttribute("Erange", Double.toString(Emax-Emin));
			xmllog.writeAttribute("ErangeInKelvin", Double.toString((Emax-Emin)*beta));

            FileWriter fileOut=null;
             if (!sFileOut.isEmpty()) {
                 fileOut = new FileWriter(new File(sFileOut));
            }

			for(double T=Tmin;T<=Tmax;T+=dT){
					ScanLnZ(T);
					calcZp_all(2);

					calcPopularity_group();

					if(fileOut!=null){
						fileOut.write(String.format("T=%fK Percent\n",T));
						for(PartFunc func:pop){
						  fileOut.write(String.format("%s %f\n",func.getTag(),Math.exp(func.lnz-(Zp[0].lnZ+refZp))));
						 }
					}

					xmllog.writeEntity("Cluster");
					xmllog.writeAttribute("T", Double.toString(T));
					xmllog.writeAttribute("Cv", Double.toString(calcCv(T)));
					for (String s : group.keySet()) {
						xmllog.writeAttribute(s, group.get(s).toString());
						group.put(s,(Double)0.0);
					}
					xmllog.endEntity().flush();
			}		

			xmllog.endEntity().flush();


		} catch (IOException ex) {
			Logger.getLogger(TaskHarmonicSuperpositionApproximation.class.getName()).log(Level.SEVERE, null, ex);
		}
		
	}

	ArrayList<PartFunc> pop=new ArrayList<PartFunc>();
	HashMap<String,Double> group=new HashMap<String,Double>();
	Zp_str[] Zp={new Zp_str(),new Zp_str(),new Zp_str(),new Zp_str(),new Zp_str()};

	double minLnZ,maxLnZ,midLnZ;
	double refZp;

	void ScanLnZ(double T){
		for(int i=0;i<pop.size();i++){
			pop.get(i).calclnZ(T);
			if(i==0){ minLnZ=maxLnZ=pop.get(0).lnz; }
			else{
				minLnZ=Math.min( minLnZ,pop.get(i).lnz);
				maxLnZ=Math.max( maxLnZ,pop.get(i).lnz);
			}
		}
		midLnZ=(minLnZ+maxLnZ)/2.0;
		refZp=maxLnZ-80.0;
		//clog<<" min = "<<minLnZ<<" max = "<<maxLnZ<<" refZp= "<<refZp<<endl;
	}

	void calcZp_all(int maxorder){
			for(int i=0;i<=maxorder;i++){
				if(Zp[i]==null) Zp[i]=new Zp_str();
				calcZp(Zp[i],i);
				//System.err.printf("Z[%d]= %f %f ; \n",i,Zp[i].lnZ,Zp[i].Z);
			}
	}

	void calcZp(Zp_str Zp1,int order_){
		int i;
		//double hs=cBeta[unit];
		Zp1.order=order_;
		if(pop.isEmpty()) return;

		Zp1.Z=0;
		for(PartFunc func : pop ){
			if(order_!=0){
				Zp1.Z+=Math.exp(func.lnz - refZp)*Math.pow(func.Ea,order_);
				//System.err.printf("		Ea= %f Z=%f ; \n",func.Ea,Zp1.Z);
			}else Zp1.Z+=Math.exp(func.lnz - refZp);
			//if(order_==1) if(indv[i].Ea<0) clog<<" Ea("<<i<<")= "<<indv[i].Ea<<" E= "<<indv[i].energy<<" Zp "<<indv[i].ZeroE<<endl;
		}
		Zp1.lnZ=Math.log(Zp1.Z);
	}

	void calcPopularity_group(){
		double TotalZ=0;
		for(PartFunc func: pop){
			double z=Math.exp(func.lnz-(Zp[0].lnZ+refZp));
			TotalZ+=z;
			if (group.containsKey(func.group)) {
					group.put(func.group, group.get(func.group) + z );
			} else {
					group.put(func.group, z );
			}
		}
		for (String s : group.keySet()) {
			group.put(s, group.get(s) / TotalZ);
		}
	}

	double calcCv(double T){
		double Cv=0;
		switch(regime){
			case CLASSICAL:
					Cv=calcCv_classical(T);
					//System.err.printf(" I am here %d %s",regime.ordinal(),regime.name());
				break;
			case QUANTUM:
					Cv=calcCv_quantum(T);
				break;
			case NOFREQ:	Cv=calcCv_classical(T);
			break;
		}
		return Cv;

	}
	
	double calcCv_classical(double T){
		double Cv=0,term1=0,term2=0;
		double hs=beta;
		term1=2.0*(Zp[1].lnZ-Zp[0].lnZ);
		term2=(Zp[2].lnZ-Zp[0].lnZ);
        int nModes=pop.get(0).getVibrationData().getNModes();
		//Cv=(indv[0].nFreq-exp(term1)+exp(term2))/indv[0].nAtom;
		Cv=(nModes/2+(-Math.exp(term1)+Math.exp(term2))*MTools.SQR(hs/T))/pop.get(0).getNAtoms();
		//clog<<" Cvgterm1= "<<term1<<" Cvgterm2= "<<term2<<" Cv= "<<Cv<<endl;
		//System.err.printf(" Cvgterm1= %f  Cvgterm2= %f  Cv=%f \n",term1,term2,Cv);
		return Cv;
	}

	double calcCv_quantum(double T){
		double Cv=0;

		Zp_str Z01=new Zp_str();
		Zp_str Z11=new Zp_str();
		Zp_str Z02=new Zp_str();

		if(pop.isEmpty()) return 0;

		Z01.Z=Z11.Z=Z02.Z=0;

		
		for(PartFunc func : pop){
			double tmp1=0,tmp2=0;
			for(double freq : func.getVibrationData().getFreqs()){
				double exprFreq=planck*freq/(Math.exp(cBetaPlanck*freq/T)-1.0);
				tmp1+=exprFreq;
				tmp2+=MTools.SQR(exprFreq)*Math.exp(cBetaPlanck*freq/T);
			}

			Z01.Z+=Math.exp(func.lnz - refZp)*tmp1;
			Z11.Z+=Math.exp(func.lnz - refZp)*(func.Ea)*tmp1;
			Z02.Z+=Math.exp(func.lnz - refZp)*(tmp2+tmp1*tmp1);

			//cerr<<T<<" "<<tmp1<<" "<<tmp2<<" "<<Z01.Z<<" "<<Z11.Z<<" "<<Z02.Z<<endl;
		}

		//Cv=-(Zp[1].Z/(Zp[0].Z))*(Zp[1].Z+Z01.Z)/(Zp[0].Z) + (Zp[2].Z+Z11.Z)/Zp[0].Z;
		//cerr<<Cv<<endl;
		Cv=(Zp[2].Z+2.0*Z11.Z+Z02.Z)/Zp[0].Z-Math.pow((Zp[1].Z+Z01.Z)/Zp[0].Z,2);
		Cv*=MTools.SQR(beta/T)/pop.get(0).getNAtoms();

		return Cv;
	}

	enum REGIME{CLASSICAL, QUANTUM, NOFREQ};

	class Zp_str{
		public double Z;
		public double lnZ;
		public int order;
	}

	class PartFunc  extends Cluster {
		@Override
		public void Clear(){
			super.Clear();
			group="";
			mi=1;
			lnz=0;
		}

		void setGroup(String arg0){
			group=arg0;
		}

		String getGroup(){
			return group;
		}

		public double calclnZ(double T){
			lnz=-1;

			double Na=Math.log(2.0)-Math.log(mi);
			for(int i=0;i<cTypeMax;i++)
				for(int j=1;j<this.nType[i];j++)	Na+=Math.log(j);

			double E_part=0,Freq_part=0;
			calcZeroE();

            double[] freqs=this.getVibrationData().getFreqs();

			switch(regime){
				case CLASSICAL:
					Freq_part+=-Math.log(cBetaPlanck/T)*freqs.length;
					for(int i=0;i<freqs.length;i++){
						if(freqs[i]>0){
							Freq_part+=- Math.log(freqs[i]);
						}else{
							//cerr<<" PartFuncPop "<<i<<"th has a negative frequency ="<<freq[i]<<". I will stop here."<<endl;
							return -1;
						}
					}
				break;

				case QUANTUM:
					for(int i=0;i<freqs.length;i++){
						if(freqs[i]>0){
							Freq_part+= -Math.log(1.0-Math.exp(-cBetaPlanck*freqs[i]/T));
							//cout<<inv->freq[i]/T<<endl;
						}else{
							//cerr<<" PartFuncPop "<<i<<"th has a negative frequency ="<<freq[i]<<". I will stop here."<<endl;
							return -1;
						}
					}
				break;

				case NOFREQ:
					Freq_part+=-Math.log(cBetaPlanck/T)*freqs.length;
				break;
			}
			Ea = energy - Emin + ZeroE;

			E_part=-beta*Ea/T;

			lnz=Na + E_part + Freq_part;
	//	#if DEBUG
				//System.err.println(sUnit);
				//System.err.printf(" Energy = %f ref = %f at %f K ",energy,Emin,T);
				//System.err.printf(" Ea=  %f  ZP=  %f  Epart= %f Freq_part= %f   lnZ(T)= %f\n",Ea,ZeroE,E_part,Freq_part,lnz);
	//	#endif

			return lnz;

		}

		double calcZeroE(){
			ZeroE=0;
            double[] freqs=this.getVibrationData().getFreqs();
			if(regime==REGIME.QUANTUM){
				for(int i=0;i<freqs.length;i++)
					ZeroE+=freqs[i];
				ZeroE*=planck/2.0;
			}
			return ZeroE;
		}

		String calcPointGroup(){
			Symmetry sym=new Symmetry();
			Symmetry.CalculateSymmetry(sym, this,0);

			PointGroup pg=sym.identifyPointGroup();

	//		Genes_str p;
			String SymCode=sym.getSymmetryCode();
			PointGroupName=pg.getGroupName();
			mi=sym.getPointGroupOrder();
	//		Migrate(p);
	//		SymCalc( p,mi, PointGroupName,SymCode);
			return PointGroupName;
		}

		int mi;
		String PointGroupName="C1";
		public double lnz=0;
		public String group="";

		public double ZeroE;
		public double Ea=0; //energy with zero point energy correction

	};

}
