/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import nqcmol.cluster.Cluster;
import java.io.*;
import java.util.Collections;
import java.util.HashMap;
import java.util.Random;
import java.util.Vector;
import java.util.logging.Level;
import java.util.logging.Logger;

import nqcmol.tools.MTools;
import org.kohsuke.args4j.*;



/**
 *
 * @author nqc
 */
public class TaskRandomGenerateCluster extends Task {

	@Option(name = "-fm", usage = "Formula of clusters. It MUST be in quote. For example: (H3O)2(H2O)4(NH3)4", metaVar = "STRING")
	String sFormula = "(H2O)4";

	@Option(name = "-num", usage = "Number of clusters. [1]", metaVar = "INTEGER")
	int num=1;

	@Option(name = "-DMin", usage = "Minimum distance of any neighbour molecules. [2.5]", metaVar = "DOUBLE")
	double DMin=2.5; 

	@Option(name = "-DMax", usage = "Maximum distance of any neighbour molecules. [4.0]", metaVar = "DOUBLE")
	double DMax=4.0;

    @Option(name = "-shape", usage = "Shape of containner where molecules are generated inside. Choose one of (BOX, SPHERE, CYLINDER). [CUBE]", metaVar = "BOOLEAN")
	Shape shape=Shape.BOX;

    @Option(name = "-a1", usage = "1st geometry parameters of container if applicant. Is radius in case of sphere and cylinder. [5]", metaVar = "DOUBLE")
	double a1=5;

    @Option(name = "-a2", usage = "2nd geometry parameters of container if applicant. Is length in case of cylinder[5]", metaVar = "DOUBLE")
	double a2=5;

    @Option(name = "-a3", usage = "3rd geometry parameters of container if applicant. [5]", metaVar = "DOUBLE")
	double a3=5;



    @Option(name = "-nocenter", usage = "Move mass center to origin. [no]", metaVar = "BOOLEAN")
	boolean bNoCenter=false;

    @Option(name = "-noflip", usage = "Donot flip molecules. [yes]", metaVar = "BOOLEAN")
	boolean bNoFlip=false;

	@Override
	public String getName(){
		return "RandomGenerateCluster";
	}

    static final public String Option="gen";

    static final public String Descriptions="\t "+Option+" \t - "+ "Randomly generate clusters\n";

	@Override
	public void ParseArguments(String[] args) {
		try {
			parser = new CmdLineParser(this);
			parser.parseArgument(args);
			String output = "";
			for (int i = 0; i < args.length-1; i++) {
                if(args[i].length()>=2)
				if (args[i].substring(0, 2).contentEquals("-o")) {
					output = args[i];
					sFileOut = args[i + 1];
				}
			}
			//System.out.println(input+" and "+output);
			for (int i = 0; i < Cluster.format.length; i++) {
				String format = "-o" + Cluster.format[i];
				if (format.contentEquals(output)) {
					sFormatOut = Cluster.format[i];
				}
				//System.out.println(format+" "+sFormatIn+" and "+sFormatOut);
			}

			if(!sFileOut.isEmpty()) fileOut=new FileWriter(new File(sFileOut));
			
		} catch (IOException ex) {
			Logger.getLogger(TaskRandomGenerateCluster.class.getName()).log(Level.SEVERE, null, ex);
		} catch (CmdLineException ex) {
			Logger.getLogger(TaskRandomGenerateCluster.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	@Override
	protected void Process() {	
        xmllog.writeAttribute("Formula", sFormula).writeAttribute("NumOfClusters", Integer.toString(num));
        xmllog.writeAttribute("a1",  String.format("%1.3f",a1));
        xmllog.writeAttribute("a2",  String.format("%1.3f",a2));
        xmllog.writeAttribute("a3",  String.format("%1.3f",a3));
        xmllog.writeAttribute("DMin", String.format("%1.3f",DMin));
        xmllog.writeAttribute("DMax",  String.format("%1.3f",DMax));
        ConstructLibrary();
        AnalyzeFormula();


        for(int i=0;i<num;i++){
            xmllog.writeEntity("Cluster").writeAttribute("id", Integer.toString(i));
            ///shuffer the positions
                if(listMol.size()>=1){
                    //for(String s :listMol)
                    //	System.out.printf("%s \n",s);

                    Collections.shuffle(listMol);
                    GenerateMolecularPositions();
                    PlaceMoleculesAndDumpOut();
                }
            xmllog.endEntity().flush();
        }
			
	}
	
	HashMap<String,Cluster> molLib=new HashMap<String,Cluster>();

    enum Shape { BOX, SPHERE, CYLINDER };
	
	protected void ConstructLibrary(){
		Cluster a;
		molLib.clear();
		
		a=new Cluster();		a.setNAtoms(3);
		a.setAtomicNumber(0,8);	a.setAtomicCoords(0,0,0,0);
		a.setAtomicNumber(1,1);	a.setAtomicCoords(1,0.94462,-0.15934,0);
		a.setAtomicNumber(2,1);	a.setAtomicCoords(2,-0.15934,0.94462,0);
		molLib.put("H2O", a);
		
		a=new Cluster();		a.setNAtoms(4);
		a.setAtomicNumber(0,8);	a.setAtomicCoords(0,-0.622681,-0.0364,-0.335512);
		a.setAtomicNumber(1,1);	a.setAtomicCoords(1,-0.990944,-0.427631,-1.16492);
		a.setAtomicNumber(2,1);	a.setAtomicCoords(2,-1.315588,0.076745,0.360097);
		a.setAtomicNumber(3,1);	a.setAtomicCoords(3,0.160787,-0.537889,-0.000937);
		molLib.put("H3O", a);


		a=new Cluster();		a.setNAtoms(2);
		a.setAtomicNumber(0,8);	a.setAtomicCoords(0,0,0,0);
		a.setAtomicNumber(1,1);	a.setAtomicCoords(1,0.94462,-0.15934,0);
		molLib.put("OH", a);

		a=new Cluster();		a.setNAtoms(4);
		a.setAtomicNumber(0,9);	a.setAtomicCoords(0,-0.01525,0.03586637,0.000);
		a.setAtomicNumber(1,1);	a.setAtomicCoords(1,0.31806893,-0.90694672,0.000);
		a.setAtomicNumber(2,1);	a.setAtomicCoords(2,0.31808614,0.50726655,0.81679);
		a.setAtomicNumber(3,1);	a.setAtomicCoords(3,0.31808614,0.50726655,-0.81679);
		molLib.put("NH3", a);

        a=new Cluster();		a.setNAtoms(3);
		a.setAtomicNumber(0,119);	a.setAtomicCoords(0, 0.010,-0.347,-0.170);
		a.setAtomicNumber(1,120);	a.setAtomicCoords(1, 0.580, 0.433,-0.420);
		a.setAtomicNumber(2,121);	a.setAtomicCoords(2,-0.590,-0.087, 0.590);
		molLib.put("spc", a);
        
        a=new Cluster();		a.setNAtoms(5);
		a.setAtomicNumber(0,119);	a.setAtomicCoords(0, 0.006, 0.016,  -0.072);
		a.setAtomicNumber(1,120);	a.setAtomicCoords(1, 0.566, 0.306,   0.648);
		a.setAtomicNumber(2,121);	a.setAtomicCoords(2,-0.624,-0.574,   0.348);
        a.setAtomicNumber(3,122);	a.setAtomicCoords(3, 0.376,-0.314,  -0.562);
		a.setAtomicNumber(4,123);	a.setAtomicCoords(4,-0.324, 0.566,  -0.362);
		molLib.put("tip5p", a);

//		for(String s : molLib.keySet())
//			System.out.printf("%s %s \n",s,molLib.containsKey(s));
	}

	int nAtoms=0;
	Vector<String> listMol;

	private void AnalyzeFormula(){
			int start = 0;
			int c = 0;
			int c1 = 0;
			listMol = new Vector<String>();
			nAtoms = 0;
			String s1 = "";
			String s = sFormula + "(";
			xmllog.writeEntity("AnalyzeFormula");
			for (int i = 0; i < s.length(); i++) {
				if (((s.charAt(i) == '(') || (s.charAt(i) == ')')) && (c > 1)) {					
					String molName = s.substring(start + 1, start+ c );
					//System.out.printf("%s - %d - %d - %s\n",s,start+1,start+ c,molName);
					start = i;
					c = 0;
					//cout<<molName<<endl;
					if (c1 % 2 == 1) {
						if(molLib.containsKey(s1)){
							for (int j = 0; j < Integer.parseInt(molName); j++) {
								listMol.add(s1);
							}
							nAtoms += Integer.parseInt(molName) * molLib.get(s1).getNAtoms();
							xmllog.writeAttribute(s1, molName);
						}else{
							xmllog.writeAttribute(s1, "Unexist");
						}
					} else {
						s1 = molName;
					}
					c1++;
				}
				c++;
			}
			xmllog.writeAttribute("nAtoms", Integer.toString(nAtoms));
			xmllog.endEntity();			
//			for(int i=0;i<listMol.size();i++){
//				System.out.printf("%s\n",listMol.elementAt(i));
//			}
		
	}

	/**
	 * Coordinate of molecules in cluster
	 */
	double[] molCoords;

	Random gen = new Random();

    private void GenerateVectorInsideContainer(double[] x){
        double radius,phi;

        switch(shape){
            case BOX:
                x[0]=a1 * (gen.nextDouble()-0.5);
                x[1]=a2 * (gen.nextDouble()-0.5);
                x[2]=a3 * (gen.nextDouble()-0.5);
            case SPHERE:
                radius=a1 * gen.nextDouble();
                MTools.generateRandomVector(x, gen, radius);
            break;
            case CYLINDER:
                radius=a1 * gen.nextDouble();
                phi=gen.nextDouble()*2*Math.PI;
                //MTools.generateRandomVector(tmp, gen, radius);
                x[0]=radius*Math.cos(phi);
                x[1]=radius*Math.sin(phi);
                x[2]=a2 * (gen.nextDouble()-0.5);
            break;
        }
    }

	/**
	 * Generate random coordinates of each molecule in cluster. Output to molCoords
	 */
	private void GenerateMolecularPositions(){
        ///generate random positions of molecules
        Cluster molTest=new Cluster();
        molTest.setNAtoms(listMol.size());
        double[] p=molTest.getCoords();
        p[0] = 0;		p[1] = 0;		p[2] = 0;
//		
        xmllog.writeEntity("GenerateMolecularPositions");
        for (int i = 1; i < listMol.size(); i++) {
            double dmin = 1000;
            double[] x=new double[3];
            do {
                GenerateVectorInsideContainer(x);
                dmin = 1000;
                for (int j = 0; j < i; j++) {
                    //double molName = sqrt(v.distSq(vlib[j]));
                    double d = molTest.distance(j, x);
                    dmin = Math.min(d, dmin);
                    //cout<<i<<" ++ "<<v<<" "<<vlib[j]<<molName<<" "<<dmin<<" "<<ranLength<<endl;
                }
            } while (!((dmin < DMax) && (dmin > DMin)));
            //cout<<format("%2d   %6.3lf ")%i%dmin<<" : "<<v<<endl;
            p[i*3+0]=x[0];
            p[i*3+1]=x[1];
            p[i*3+2]=x[2];

            xmllog.writeEntity("Molecule").writeAttribute("Type",listMol.get(i));
            xmllog.writeAttribute("Dmin", String.format("%1.3f", dmin));
            xmllog.writeAttribute("x", Double.toString(x[0]));
            xmllog.writeAttribute("y", Double.toString(x[0]));
            xmllog.writeAttribute("z", Double.toString(x[0]));
            xmllog.endEntity();
        }

        xmllog.endEntity();

        molCoords=p.clone();
	}

	/**
	 * Place molecules at calculated coordinates and dumpt out the results.
	 */
	private void PlaceMoleculesAndDumpOut(){
        int c=0;
        Cluster result = new Cluster();
        result.setNAtoms(nAtoms);
                    result.setTag(Integer.toString(c));
                    c++;
        xmllog.writeEntity("PlaceMoleculesAndDumpOut");
        int pos = 0;
        for (int i = 0; i < listMol.size(); i++) {
            Cluster mol = (Cluster) molLib.get(listMol.get(i)).clone();
            //mol.Write(System.err, "xyz");
            //randomly rotate
            double[] axis={gen.nextDouble(), gen.nextDouble(), gen.nextDouble()};
            if(!bNoFlip){
                double angle = 2.0 * Math.PI * gen.nextDouble();
                mol.RotateAxis(angle, axis);
            }

            axis[0] = molCoords[i * 3 + 0];
            axis[1] = molCoords[i * 3 + 1];
            axis[2] = molCoords[i * 3 + 2];

            mol.Translate(axis, 1);
            result.replaceMolecule(pos, mol);
            pos += mol.getNAtoms();
        }

        xmllog.endEntity();

        if(!bNoCenter) result.Center();

        if (fileOut != null) {
            result.Write(fileOut, sFormatOut);
        }
	}
}

