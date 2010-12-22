/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
import nqcmol.cluster.MolExtra;
import java.io.*;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.potential.GaussianInterfacePotential;
import nqcmol.potential.Potential;
import nqcmol.tools.MTools;
import org.apache.commons.math.FunctionEvaluationException;
import org.apache.commons.math.analysis.polynomials.PolynomialFunction;
import org.apache.commons.math.optimization.OptimizationException;
import org.apache.commons.math.optimization.fitting.ParametricRealFunction;
import org.apache.commons.math.optimization.fitting.PolynomialFitter;
import org.apache.commons.math.optimization.general.LevenbergMarquardtOptimizer;
import org.kohsuke.args4j.Option;


import org.xml.sax.SAXException;

// DOM classes.
import org.w3c.dom.*;
//JAXP 1.1
import javax.xml.parsers.*;
import javax.xml.transform.*;
import javax.xml.transform.stream.*;
import javax.xml.transform.dom.*;
import nqcmol.cluster.VibrationData;
import nqcmol.tools.XmlWriter;


/**
 *
 * @author nqc
 */
public class TaskAnharmonicVibrationAnalysis extends Task {

	@Option(name="-a",usage="parameter for Gaussian. It should be in this format \"method@charge@multiplicity\". For example: \"b3lyp/6-31+G*@0@1\" (note the quote)",metaVar="STRING")
    String sFileParam="hf/3-21G@0@1";

    //@Option(name="-xml",usage="input XML file generated from Gaussian scaling. If specified, the program will perform fitting instead of sanning ",metaVar="STRING")
    //String sFileXML="";

	private String basisSet="";
	private String chargeAndMultiplicity="0 1";
	private String sFileCheckPoint="checkpoint.chk";

	Potential pot;
	GaussianInterfacePotential gau=null;

	@Override
	public String getName(){
		return "AnharmonicVibrationAnalysis";
	}

    static final public String Option="anvib";

    static final public String Descriptions="\t "+Option+" \t - "+ "anharmonic vibration analysis\n";

    public static void main(String[] args) throws IOException  {
        new TaskAnharmonicVibrationAnalysis().Execute(args);
	}

   	@Override
	protected void Initialize() {		
        super.Initialize();
        if(!sFormatIn.contentEquals("xml")){
            pot = MolExtra.SetupPotential("g03", sFileParam, "Hartree", "");
            assert (pot.getEquation().contains("Gaussian"));
            gau = (GaussianInterfacePotential) pot;
            xmllog.writeNormalText(gau.XMLInfo(1));
            basisSet = gau.getBasisSet();
            chargeAndMultiplicity = gau.getChargeAndMultiplicity();

           int dot = sFileIn.lastIndexOf('.');
          // String ext=sFileIn.substring(dot + 1);
           int sep = sFileIn.lastIndexOf('/');
           sFileCheckPoint=sFileIn.substring(sep + 1, dot);
        }
        //String filename = (new File(sFileIn)).getName();
        //String ext = (filename.lastIndexOf(".") == -1) ? "" : filename.substring(filename.lastIndexOf(".") + 1, filename.length());
        //sFileCheckPoint = filename.replace("." + ext, "");
	}
   

	@Override
	protected void Process() {
		try {
            Scanner fileIn = new Scanner(new File(sFileIn));

			xmllog.writeEntity("Note");
			xmllog.writeText(" There are two steps to calculate anharmonic frequencies. Firstly, perform Gaussian scanning  and pipeline the results to an XML file. Then, use XML file generated to approximate frequencies. Remember to specify the output file at the second step.");
			xmllog.endEntity();
            if(!sFormatIn.contentEquals("xml")){//Gaussian scanning
                if (cluster.Read(fileIn, "g03")) {
                    //Reading normal modes
                    cluster.setTag((new File(sFileIn)).getName());
                    xmllog.writeEntity("Cluster").writeAttribute("Tag", cluster.getTag());
                    xmllog.writeAttribute("nAtom", Integer.toString(cluster.getNAtoms()));
                    xmllog.writeAttribute("Energy", Double.toString(cluster.getEnergy()));
                    VibrationData vib=cluster.getVibrationData();
                    if(vib!=null){
                        int iFrom=iMode;
                        int iTo=Math.min(iMode + nModes,vib.getFreqs().length);
                        if(iMode==-1){  iFrom=0; iTo=vib.getFreqs().length;}
                            
                        for(int i=iFrom;i<iTo;i++)  ScanningAlongNormalModes(i);
                    }
                    xmllog.endEntity().flush();
                }
            }else{// perform fitting
                ParseXMLFile();
            }

			fileIn.close();
		} catch (IOException ex) {
			Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	///===================
	@Option(name = "-np", usage = "number of processors using in Gaussian. Default is 1.", metaVar = "INTEGER")
	int numOfProcessors=1;

    @Option(name = "-mem", usage = "Memory allocated for Gaussian.", metaVar = "STRING")
	String mem="500MB";

	@Option(name = "-num", usage = "number of points using in approximation. Default is 5.", metaVar = "INTEGER")
	int nPoints=5;
	

	@Option(name = "-m", usage = "starting mode to be calculated (count from 0), number of calculated modes is determined from -nmodes option. If -1 (default), all normal modes will be calculated.", metaVar = "INTEGER")
	int iMode=-1;

	@Option(name = "-nmodes", usage = "number of modes will be calculated (default is 1). The starting mode is determined from -m option", metaVar = "INTEGER")
	int nModes=1;

	@Option(name = "-degree", usage = "Maximum degree in approximation. Default is 2.", metaVar = "INTEGER")
	int degree=2;

	@Option(name = "-maxX", usage = "maximum displacement in approximation", metaVar = "INTEGER")
	double maxX=0.2;

	Cluster cluster=new Cluster();

	void ScanningAlongNormalModes(int m) throws IOException{
        double[] deltaX=new double[nPoints];
        for(int i=0;i<nPoints;i++){
			//temporary scaling factor
			double beta=(2.0*i/(nPoints-1)-1)*maxX;
			deltaX[i]=beta;
        }

        VibrationData vib=cluster.getVibrationData();

		xmllog.writeEntity("NormalMode").writeAttribute("id", Integer.toString(m));
			if((m>=vib.getFreqs().length)||(m<0))
				xmllog.writeText(String.format("Error. Mode %d is not valid because it is not in between [0,%d)",m,vib.getFreqs().length));
			else{
				long duration =  System.currentTimeMillis();
				xmllog.writeAttribute("Freq", Double.toString(vib.getFreqs(m)));
				xmllog.writeAttribute("ReducedMass", Double.toString(vib.getReducedMass(m)));
				xmllog.writeEntity("Vector");
				for (int k = 0; k < cluster.getNAtoms(); k++) {
					xmllog.writeEntity("Atom").writeAttribute("id", Integer.toString(k));
					xmllog.writeAttribute("x", Double.toString(vib.getNormalModeVectors(m)[k * 3 + 0]));
					xmllog.writeAttribute("y", Double.toString(vib.getNormalModeVectors(m)[k * 3 + 1]));
					xmllog.writeAttribute("z", Double.toString(vib.getNormalModeVectors(m)[k * 3 + 2]));
					xmllog.endEntity().flush();
				}
				xmllog.endEntity().flush();
				//generate gaussian input
				xmllog.writeEntity("GaussianScanning");
				String sGaussianInput = GenGaussianInputFromNormalModes(m,deltaX);
				//System.err.println(sGaussianInput);
				gau.setCluster(cluster);

                double[] energies=null;
				energies = gau.getEnergies(sGaussianInput);
				for (int j = 0; j < energies.length; j++) {
					xmllog.writeEntity("Step").writeAttribute("id", Integer.toString(j+1));
					xmllog.writeAttribute("DeltaX", Double.toString(deltaX[j]));
					xmllog.writeAttribute("Energy", Double.toString(energies[j]));
					xmllog.endEntity().flush();
				}
                xmllog.writeEntity("Duration");
					xmllog.writeNormalText(Double.toString((System.currentTimeMillis()-duration)/1000.0));
				xmllog.endEntity().flush();

				xmllog.endEntity().flush();			
			}
		xmllog.endEntity().flush();
	}

      

    void ParseXMLFile(){
        try {
            File file = new File(sFileIn);
            Document doc = DocumentBuilderFactory.newInstance().newDocumentBuilder().parse(file);
            doc.getDocumentElement().normalize();
            
            //System.out.println("Root element :" + doc.getDocumentElement().getNodeName());
            //System.out.println("-----------------------");

            NodeList listCluster = doc.getElementsByTagName("Cluster");
           
            for (int iCluster = 0; iCluster < listCluster.getLength(); iCluster++) {//for each cluster

               Element nodeCluster = (Element)listCluster.item(iCluster);
               //System.out.println("Tag : "  + nodeCluster.getAttributes().getNamedItem("Tag"));
               //System.out.println("nAtom : "  + nodeCluster.getAttributes().getNamedItem("nAtom"));
               //System.out.println("Energy : "  + nodeCluster.getAttributes().getNamedItem("Energy"));
               
               double energy= XmlWriter.getDoubleAttributeXML("Energy",nodeCluster); //reading energy

               NodeList listNormalMode=nodeCluster.getElementsByTagName("NormalMode");             

               for(int iNormalMode=0;iNormalMode < listNormalMode.getLength();iNormalMode++){//for each normal mode

                  // System.out.println("\tNormalMode : "  + iNormalMode);
                   Element nodeNormalMode=(Element)listNormalMode.item(iNormalMode);
                   double reducedMass= XmlWriter.getDoubleAttributeXML("ReducedMass",nodeNormalMode); //reading reduced mass

                   //System.out.println("\t\tReducedMass : "  + reducedMass);

                   Element nodeGaussianScanning=(Element)nodeNormalMode.getElementsByTagName("GaussianScanning").item(0);
                   NodeList listStep=nodeGaussianScanning.getElementsByTagName("Step");

                   //double[] deltaX=new double[listStep.getLength()+1];
                   //double[] energies=new double[listStep.getLength()+1];

                   double[] deltaX=new double[listStep.getLength()];
                   double[] energies=new double[listStep.getLength()];

                   //deltaX[0]=0; energies[0]=energy;
                   for(int iStep=0;iStep < listStep.getLength();iStep++){//for every steps                        
                        Element nodeStep= (Element)listStep.item(iStep);
                        deltaX[iStep]  = XmlWriter.getDoubleAttributeXML("DeltaX",nodeStep);
                        energies[iStep]= XmlWriter.getDoubleAttributeXML("Energy",nodeStep);
                       //System.out.println("\t\tStep : "  + iStep + " " + deltaX[iStep]+" "+energies[iStep]);
                   }

                  // nodeGaussianScanning.getParentNode().removeChild(nodeGaussianScanning);

                   
                   PolynomialFitting(deltaX,energies,reducedMass,nodeNormalMode);
                   MorseFitting(deltaX,energies,reducedMass,energy,nodeNormalMode);
               }
            }

            //writing out the results
            XmlWriter.writeDocXML(doc, sFileOut);

        } catch (SAXException ex) {
            Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
        } catch (IOException ex) {
            Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
        } catch (ParserConfigurationException ex) {
            Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
        }
    }

    /**
	 * Fitting polynomial potential with observation (x,y)=(delta,energies)
	 * @param params params of Morse potential, will be initialized and updated after fitting
	 * @param calcE values of energies calculated by Morse potential,will be initialized and updated after fitting
	 */
    void PolynomialFitting(double[] deltaX,double[] energies,double reducedMass,Element nodeNormalMode){
        try {
            //Polynomial approximation            

            //start fitting
                PolynomialFitter fitter = new PolynomialFitter(degree, new LevenbergMarquardtOptimizer());
                //fitter.addObservedPoint(1, 0, cluster.getEnergy());
                for (int i = 0; i < deltaX.length; i++) {
                    //System.out.println(" "+i+" "+deltaX[i]+" "+energies[i]);
                    fitter.addObservedPoint(1, deltaX[i], energies[i]);
                }
                PolynomialFunction fittedFunc=null;

                fittedFunc = fitter.fit();
                double[] params=fittedFunc.getCoefficients(); //System.out.println(" length=  "+params.length);
                double[] calcE=new double[deltaX.length];
                for (int j = 0; j < deltaX.length; j++) {
                    calcE[j]=fittedFunc.value(deltaX[j]);
                }
			          
            //double reducedMass = cluster.getReducedMass(m); //in amu
            double forceConst = params[2] * 2; //in Hartree.Angstrom^-2
            double kConvert = Math.sqrt(5.4857990943e-4) * 1e8 * 7.297352569e-3;

            double omega = kConvert * Math.sqrt(forceConst / reducedMass) / (2.0 * Math.PI); //in cm-1


            //remove the old fit if exists
            XmlWriter.removeNodeXML("PolynomialFit",nodeNormalMode);
            //create a new one
            Element nodePoly=XmlWriter.createNodeXML("PolynomialFit",nodeNormalMode);

            //add result
            XmlWriter.addAttributeXML("HarFreq",Double.toString(omega),nodePoly);


            for (int j = 0; j < params.length; j++) {//details of parameters
                Element nodeTerm=XmlWriter.createNodeXML("Term",nodePoly);
                XmlWriter.addAttributeXML("Degree", Integer.toString(j),nodeTerm);
                XmlWriter.addAttributeXML("Coefficient",  Double.toString(params[j]),nodeTerm);
            }
            
//            xmllog.writeEntity("PolynomialFitting");
//            xmllog.writeAttribute("AnFreq", Double.toString(anFreq));//result
//
//            for (int j = 0; j < params.length; j++) {//details of parameters
//                xmllog.writeEntity("Term").writeAttribute("Degree", Integer.toString(j));
//                xmllog.writeAttribute("Coefficient", Double.toString(params[j]));
//                xmllog.endEntity().flush();
//            }
            //writing validity of fitting: rmsE etc
            WriteValidity(deltaX,energies,calcE,nodePoly);
            
            //xmllog.endEntity().flush();
        } catch (OptimizationException ex) {
            Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
//        } catch (IOException ex) {
//            Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
        }
    }
   
	/**
	 * Fitting Morse potential with observation (x,y)=(delta,energies)
	 * @param params params of Morse potential, will be initialized and updated after fitting
	 * @param calcE values of energies calculated by Morse potential,will be initialized and updated after fitting
	 */
	void MorseFitting(double[] deltaX,double energies[],double reducedMass,double energy,Element nodeNormalMode) {
			try {
                //Morse fitting
               

                MorsePotential fittedFunc=new MorsePotential();
                double[] params={1,2.18};
//                CurveFitter fitter = new CurveFitter(new LevenbergMarquardtOptimizer());
//                for (int i = 0; i < deltaX.length; i++) {
//                    fitter.addObservedPoint(1, deltaX[i], energies[i]-energy);
//                }
                
//                
//                params = fitter.fit(fittedFunc, params);

                FileWriter fileTMP=new FileWriter(new File("data.tmp"));
                fileTMP.write("0.0\t%0.0f\n");
                for (int i = 0; i < deltaX.length; i++) {
                    fileTMP.write(String.format("%12.8f\t%12.8f\n", deltaX[i], energies[i]-energy));
                }
                fileTMP.close();

                Runtime r = Runtime.getRuntime();
                //String gnuplotInput="";
                String command = "gnuplot -e \" De=1;alpha=1;f(x)=De*(1-exp(-alpha*x))**2;fit f(x) \\\"data.tmp\\\" using 1:2 via De,alpha; \"";
                Process p = r.exec(new String[]{"/bin/bash", "-c",command});

                //System.out.println(command);

                BufferedReader in = new BufferedReader(new InputStreamReader(p.getErrorStream()));
                String s;
                while ((s = in.readLine()) != null) {
                    if(s.contains("Final set of parameters            Asymptotic Standard Error")){
                        in.readLine();in.readLine();
                        s=in.readLine();
                        int equalSign=s.indexOf("="); int plusSign=s.indexOf("+");
                        params[0]=Double.parseDouble(s.substring(equalSign+1, plusSign-1).trim());
                        s=in.readLine();
                        equalSign=s.indexOf("=");     plusSign=s.indexOf("+");
                        params[1]=Double.parseDouble(s.substring(equalSign+1, plusSign-1).trim());

                    }
                    //System.out.println(s);
                }

                try {
				if (p.waitFor() != 0) {
					System.err.println(" Error termination with GNUPLOT. Exit value = " + p.exitValue());
				}
                } catch (InterruptedException e) {
                    System.err.println(e);
                } finally {
                    // Close the InputStream
                    in.close();
                }


				double[] calcE=new double[energies.length];
				for (int j = 0; j < energies.length; j++) {
					calcE[j]=fittedFunc.value(deltaX[j],params)+energy;
				}

                //=====
                double De=params[0];//in Hartree
                double alpha=params[1];// in AA^-1
                double M=reducedMass;//in a.u

                double forceConst=2*De*Math.pow(alpha,2);
                double kConvert = Math.sqrt(5.4857990943e-4) * 1e8 * 7.297352569e-3;
                double omega = kConvert * Math.sqrt(forceConst / M) / (2.0 * Math.PI); //in cm-1
                
                double ksiomega=alpha*alpha/(M*0.2436);//a*Math.sqrt(2.0*De/reducedMass)/(2*Math.PI);//a/2pi*sqrt(2De/m)
                double anFreq=omega-ksiomega*2.0;

                //remove the old fit if exists
                XmlWriter.removeNodeXML("MorseFit",nodeNormalMode);
                //create a new one
                Element nodeMorse= XmlWriter.createNodeXML("MorseFit",nodeNormalMode);

                //add result
                XmlWriter.addAttributeXML("De",Double.toString(De),nodeMorse);
                XmlWriter.addAttributeXML("alpha",Double.toString(alpha),nodeMorse);
                XmlWriter.addAttributeXML("HarFreq",Double.toString(omega),nodeMorse);
                XmlWriter.addAttributeXML("AnharConstant",Double.toString(ksiomega),nodeMorse);
                XmlWriter.addAttributeXML("AnharFreq",Double.toString(anFreq),nodeMorse);


//                xmllog.writeEntity("MorseFitting");
//                xmllog.writeAttribute("AnFreq", Double.toString(anFreq));//result
//                xmllog.writeAttribute("De", Double.toString(De));
//                xmllog.writeAttribute("alpha", Double.toString(alpha));
//                xmllog.writeAttribute("AnharmonicityConstant", Double.toString(ksiomega));

                //writing validity of fitting: rmsE etc
                WriteValidity(deltaX,energies,calcE,nodeMorse);

                //xmllog.endEntity().flush();
			} catch (IOException ex) {
            Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
        } catch (FunctionEvaluationException ex) {
				Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
			} catch (IllegalArgumentException ex) {
				Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
			}
	}

	class MorsePotential implements ParametricRealFunction{

		@Override
		public double value(double x, double[] param) throws FunctionEvaluationException {
			double De=param[0];
			double a=param[1];
			double re=0;
			return De*Math.pow(1-Math.exp(a*(re-x)),2);
		}

		@Override
		public double[] gradient(double x, double[] param) throws FunctionEvaluationException {
			double De=param[0];
			double a=param[1];
			double re=0;
			double[] gradient=new double[2];

			gradient[0]=Math.pow(1-Math.exp(a*(re-x)),2);
			gradient[1]=2.0*De*(1-Math.exp(a*(re-x)))*(-(re-x)*Math.exp(a*(re-x)));
			return gradient;
		}
		
	}

    double WriteValidity(double[] deltaX,double[] energies,double[] calcE,Element node){
		double rmsE = 0;
		//try {
			for (int j = 0; j < energies.length; j++) {
				double deltaE = Math.abs(calcE[j] - energies[j]);
				rmsE += deltaE * deltaE;
			}
			rmsE = Math.sqrt(rmsE / energies.length);

            //create a new one
            Element nodeValidity=XmlWriter.createNodeXML("Validity",node);

            //add result
            XmlWriter.addAttributeXML("RMSE",Double.toString(rmsE),nodeValidity);

            for (int j = 0; j < energies.length; j++) {
				double deltaE = Math.abs(calcE[j] - energies[j]);
                Element nodeStep=XmlWriter.createNodeXML("Step",nodeValidity);
                XmlWriter.addAttributeXML("DeltaE",Double.toString(deltaE),nodeStep);
                XmlWriter.addAttributeXML("CalcE",Double.toString(calcE[j]),nodeStep);
                XmlWriter.addAttributeXML("ObservedE",Double.toString(energies[j]),nodeStep);
                XmlWriter.addAttributeXML("DeltaX",Double.toString(deltaX[j]),nodeStep);
                XmlWriter.addAttributeXML("id",Integer.toString(j),nodeStep);

			}
           
//            xmllog.writeEntity("Validity");
//			xmllog.writeAttribute("RMSE", Double.toString(rmsE));
//			for (int j = 0; j < energies.length; j++) {
//				double deltaE = Math.abs(calcE[j] - energies[j]);
//				xmllog.writeEntity("Step").writeAttribute("id", Integer.toString(j));
//				xmllog.writeAttribute("DeltaX", Double.toString(deltaX[j]));
//				xmllog.writeAttribute("ObservedE", Double.toString(energies[j]));
//				xmllog.writeAttribute("CalcE", Double.toString(calcE[j]));
//				xmllog.writeAttribute("DeltaE", Double.toString(deltaE));
//				xmllog.endEntity().flush();
//			}
//			xmllog.endEntity().flush();
//		} catch (IOException ex) {
//			Logger.getLogger(TaskAnharmonicVibrationAnalysis.class.getName()).log(Level.SEVERE, null, ex);
//		}
		return rmsE;
	}

	String GenGaussianInputFromNormalModes(int m,double[] delta){
		double[] coord=cluster.getCoords();
		double[] direction=cluster.getVibrationData().getNormalModeVectors()[m];


		String input="";
		
		for(int i=0;i<delta.length;i++){
			//temporary scaling factor
			double beta=delta[i];

			String g03header="%chk="+sFileCheckPoint+String.format("-m%d.chk", m);
			g03header+=" \n%mem="+mem+" \n%nproc="+String.format("%d", numOfProcessors)+" \n#p "+basisSet;
			//g03header+=" \n%nproc="+String.format("%d", numOfProcessors)+" \n#p "+basisSet;


			if(i!=0)	g03header="--Link1--\n"+	g03header + " Guess=Read ";
			input+=g03header+"\n";
			input+=String.format("\nM= %d i= %d beta= %1.5f\n\n",m,i,beta);

			input+=chargeAndMultiplicity+"\n";//Charge and multiplicity

			double[] newCoords=new double[coord.length];

			MTools.VEC_PLUS_VEC(newCoords, coord, direction, 1.0, beta);

			for(int j=0;j<cluster.getNAtoms();j++){
				input+=cluster.getAtomicSymbol(j)+" ";
				input+=String.format("%12.8f %12.8f %12.8f\n",newCoords[j*3+0],newCoords[j*3+1],newCoords[j*3+2]);
			}

			input+="\n";
		}

		return input;
	}

	
}
