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
	@Option(name="-t",usage="Please choose a task",metaVar="task")
    String sTask="ener";

    @Option(name="-h",usage="Print out the help")
    boolean isHelp= false;

    /**
     * @param args the command line arguments
     */
    public static void main(String[] args) throws Exception {
		 new NQCMol().doMain(args);
	}

	public void doMain(String[] args) throws IOException {
		CmdLineParser parser = new CmdLineParser(this);		 parser.setUsageWidth(80);
		try {
			// parse the arguments.
			parser.parseArgument(args);
			//parser.parseArgument("more","args");
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
			MCalculate calc=new MCalculate();
			calc.CalculateEnergy(args);
			return;
		}


		if(sTask.contentEquals("validgrad")){
			MCalculate calc=new MCalculate();
			calc.ValidateGradients(args);
			return;
		}

		if(sTask.contentEquals("opt")){
			MCalculate calc=new MCalculate();
			calc.Optimize(args);
			return;
		}

		if(sTask.contentEquals("vib")){
			MCalculate calc=new MCalculate();
			calc.HarmonicVibrationAnalysis(args);
			return;
		}

		if(isHelp){
			parser.printUsage(System.out);
		}


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
//			System.out.println("\n Now evaluate the energy \n");
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
