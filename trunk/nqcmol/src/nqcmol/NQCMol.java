/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import org.kohsuke.args4j.*;


/**
 *
 * @author nqc
 */
public class NQCMol {
	@Option(name="-t",usage="Please choose one of those operations\n" +
		"\t opt \t - optimize clusters\n" +
		"\t ener \t - calculate energy\n" +
		"\t grad \t - validate gradients. Perform both analytical (if applicant) and numerical gradients and calculate the residual\n" +
		"\t vib \t - calculate harmonic frequencies\n" +
		"\t be \t - calculate binding energies\n" +
		"\t sort \t - sort clusters according to energies\n" +
		"\t class \t - classify clusters according to topologies\n" +
		"\t hsa \t - harmonic superposition approximation\n" +
		"\t equiAtom \t - detect topologically equivalent atoms in clusters\n" +
		"\t sim \t - remove duplicate clusters\n",metaVar="OPER")
    String sTask="";

    @Option(name="-h",usage="Print out the help. Run nqcmol -t [oper] -h for the details of operation.")
    boolean isHelp= false;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {
		 new NQCMol().doMain(args);
	}

	public void doMain(String[] args) throws IOException {
		CmdLineParser parser = new CmdLineParser(this);		 parser.setUsageWidth(120);
		try {
			// parse the arguments.
			parser.parseArgument(args);
		} catch (CmdLineException e) {
			// if there's a problem in the command line, you'll get this exception. this will report an error message.
            System.err.println(e.getMessage());            
            // print the list of available options
            parser.printUsage(System.err);
            System.err.println();

            return;
		}
		//parser.parseArgument("more","args");
	
		if(sTask.contentEquals("ener")){
			Task calc=new TaskCalculateEnergy();
			calc.Execute(args);			return;
		}


		else if(sTask.contentEquals("grad")){
			Task calc=new TaskValidateGradient();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("opt")){
			Task calc=new TaskOptimizeCluster();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("vib")){
			Task calc=new TaskHarmonicVibrationAnalysis();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("anvib")){
			Task calc=new TaskAnharmonicVibrationAnalysis();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("be")){
			Task calc=new TaskCalculateBindingEnergy();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("sym")){
			Task calc=new TaskCalculateSymmetryPointGroup();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("equiAtom")){
			Task calc=new TaskDetectEquivalentAtoms();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("corr")){
			Task calc=new TaskAutoCorrelationAndTransform();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("test")){
			Task calc=new TaskTest();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("gen")){
			Task calc=new TaskRandomGenerateCluster();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("class")){
			Task calc=new TaskClassifyCluster();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("sort")){
			Task calc=new TaskSortCluster();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("sim")){
			Task calc=new TaskRemoveDuplicateCluster();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("fit")){
			TaskFitPotential calc=new TaskFitPotential();
			calc.Execute(args);			return;
		}

		else if(sTask.contentEquals("hsa")){
			Task calc=new TaskHarmonicSuperpositionApproximation();
			calc.Execute(args);			return;
		}


		else if(isHelp){
			System.out.println(" nqcmol - utilities for processing data. Author: Nguyen Quoc Chinh\n");
			//System.out.println(" USAGE\n\t\t nqcmol [OPTION]\n");
			System.out.println(" OPTION\n");
			parser.printUsage(System.out);
			return;
		}

		System.out.print(" You did not choose anything!!\n");

			//BufferedReader buff = new BufferedReader(new FileReader("test/LJlm/lj13.xyz"));
//			FileReader buff = new FileReader("test/LJlm/lj13.xyz");
//			XYZReader reader=new XYZReader(buff);
//			ChemFile chemFile;
//
//			chemFile = (ChemFile) reader.read((ChemObject) new ChemFile());
//			IChemSequence seq = chemFile.getChemSequence(0);
//			IChemModel model = seq.getChemModel(0);
//			IMoleculeSet som = model.getMoleculeSet();
//			IMolecule mol= som.getMolecule(0);
//
//			System.out.print(mol);
//
//
//			FileWriter writer = new FileWriter("test/LJlm/new-lj13.xyz");
//			XYZWriter xyzWriter = new XYZWriter(writer);
//			xyzWriter.write(mol);
//			xyzWriter.close();
//			writer.close();
//
//			System.out.println("\n Now Evaluate the energy \n");
//			LennardJonesFunction pot = new LennardJonesFunction();
//			double energy = pot.energyFunctionOfAMolecule(mol);
//
//			System.out.println(" Energy = "+ Double.toString(energy));
//
//			GVector molecule3Coordinates  = new GVector(mol.getAtomCount()*3);
//			molecule3Coordinates=ForceFieldTools.getCoordinates3xNVector(mol);
			//GeometricMinimizer gm = new GeometricMinimizer();
			//gm.setConvergenceParametersForCGM(100, 0.00001);
	 		//gm.conjugateGradientMinimization(molecule3Coordinates, pot);

    }

	
}