/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Crystal;
import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;



import org.kohsuke.args4j.*;

import nqcmol.tools.MTools;

/**
 *
 * @author nqc
 */
public class TaskCalculateRadialDistributionFunction extends Task{

    @Option(name="-rmin",usage="radial length in &Aring;ngstrom at which the RDF starts",metaVar="DOUBLE")
    double startCutoff= 0;

    @Option(name="-rmax",usage="radial length in &Aring;ngstrom at which the RDF stops",metaVar="DOUBLE")
    double cutoff= 2.0;

    @Option(name="-binsize",usage="width of the bins",metaVar="DOUBLE")
    double resolution= 1e-1;

    @Option(name="-iAtom",usage="index of atom at which RDF starts",metaVar="INTEGER")
    int iAtom= 0;

	@Override
	public String getName(){
		return "CalculateRadialDistributionFunction";
	}

    static final public String Option="rdf";

    static final public String Descriptions="\t "+Option+" \t - "+ "Calculate RDF function\n";

    public static void main(String[] args) throws IOException  {
        new TaskCalculateRadialDistributionFunction().Execute(args);
	}

    private double peakWidth=0;

    int binsToFillOnEachSide ;

    double[] rdf;

    double[] factors;

	@Override
	protected void Process() {
		try {			
			//System.out.println("What the fuck!"+sFileIn);
			int i = 0;
            int length = (int)((cutoff-startCutoff)/resolution) + 1;

        //logger.debug("Creating RDF of length ", length);

            // the next we need for Gaussian smoothing
            binsToFillOnEachSide = (int)(peakWidth*3.0/resolution);
            double sigmaSquare = Math.pow(peakWidth, 2.0);
            factors = new double[binsToFillOnEachSide];
            for (int binCounter=0; binCounter<binsToFillOnEachSide; binCounter++) {
                factors[binCounter] = Math.exp(-1.0*(Math.pow(((double)binCounter)*resolution, 2.0))/sigmaSquare);
            }

            // this we need always
           rdf = new double[length];

			Crystal cryst=new Crystal();
			while (cryst.Read(fileIn, sFormatIn)) {			
				accumulate(cryst);
				xmllog.writeEntity("Cluster");
				xmllog.writeAttribute("id", Integer.toString(i));
                cryst.setCellAngle(90, 90, 90);
                cryst.setCellLength(12.0,12.0, 12.0);
//				xmllog.writeAttribute("nAtoms", Integer.toString(mol.getNAtoms()));
//				xmllog.writeAttribute("SymmetryCode", sym.getSymmetryCode());
//				PointGroup pg=sym.identifyPointGroup();
//				if(pg!=null)
//					xmllog.writeAttribute("PointGroup", pg.getGroupName());
//				else
//					xmllog.writeAttribute("PointGroup", "Unknown");
//
//				xmllog.writeAttribute("PointGroupOrder", Integer.toString(sym.getPointGroupOrder()));
//
				xmllog.endEntity().flush();
				i++;
				break;
			}
			fileIn.close();

//            xmllog.writeEntity("RDF");
//            for(int i=0;i<rdf.length;i++){
//                xmllog.writeEntity("RDF");
//                xmllog.writeAttribute(sFileIn, sFileIn)
//            }
//            xmllog.endEntity().flush();

            if(fileOut!=null){
                for(i=0;i<rdf.length;i++){
                    fileOut.write(String.format("%f\t%f\n",startCutoff+resolution*(i+0.5),rdf[i]));

                }
            }

		} catch (IOException ex) {
			Logger.getLogger(TaskCalculateRadialDistributionFunction.class.getName()).log(Level.SEVERE, null, ex);
		} 
	}



    /**
     * Calculates a RDF for <code>Atom</code> atom in the environment
     * of the atoms in the <code>AtomContainer</code>.
     */
    public double[] accumulate(Crystal cryst) {

        int N1max=(int)(cutoff/cryst.getCellLength(0))+1;
        int N2max=(int)(cutoff/cryst.getCellLength(1))+1;
        int N3max=(int)(cutoff/cryst.getCellLength(2))+1;

        //double[] x0={cryst.getCoords(iAtom*3),cryst.getCoords(iAtom*3+1),cryst.getCoords(iAtom*3+2)};

        for (int i=0; i<cryst.getNAtoms(); i++) if(cryst.getAtomicNumber(i)==8){
            double[] vec=new double[3];
            cryst.dR(iAtom, i, vec);

            for(int l=-N1max;l<=N1max;l++)
                for(int m=-N2max;m<=N2max;m++)
                    for(int n=-N3max;n<=N3max;n++){
                        double distance = Math.sqrt(MTools.DOTPRODUCT(vec, vec));
                        int index = (int)((distance-startCutoff)/this.resolution);

                        double weight=1;
                        if(index<rdf.length){
                            rdf[index] += weight; // unweighted

                            if (this.peakWidth > 0.0) {
                                // apply Gaussian smoothing
                                for (int binCounter=1; binCounter<=binsToFillOnEachSide; binCounter++) {
                                    if ((index - binCounter) >= 0) {
                                        rdf[index - binCounter] += weight*factors[binCounter];
                                    }
                                    if ((index + binCounter) < rdf.length) {
                                        rdf[index + binCounter] += weight*factors[binCounter];
                                    }
                                }
                            }
                        }
                    }
        }
        return rdf;
    }

	
}
