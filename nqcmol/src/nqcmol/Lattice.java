/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.io.*;
import java.text.DateFormat;
import java.text.SimpleDateFormat;
import java.util.Scanner;
import java.util.StringTokenizer;
import java.util.Vector;
import java.util.logging.*;
import nqcmol.tools.MTools;
import nqcmol.tools.XmlWriter;

/**
 *
 * @author nqc
 */
public class Lattice extends Cluster{
	public Lattice(){
	};
		
	public void Get(Lattice src){
		//System.out.println(" Come here");
		setNAtoms(src.getNAtoms());
		//System.out.println(" Come here");
		setCoords(src.getCoords());
		//System.out.println(" Come here");
		setGradient(src.getGradient());
		//System.out.println(" Come here");
		setHessian(src.getHessian());
		//System.out.println(" Come here");
		setFreqs(src.getFreqs());
		//System.out.println(" Come here");
		setAtomicNumber(src.getAtomicNumber());
		//System.out.println(" Come here");
		tag=src.getTag();
		rmsGrad=src.getRmsGrad();	
	}
/*
    int ExpandToSuperCell(int no,int dir){
        ConvertToFract(false);
        Individual cell0,cell;
        Migrate(cell0);

        int d=(no>0)?1:-1;
        for(int kx=1;kx<=fabs(no);kx++){
            cell0.Migrate(cell);
            double v[3];
            int l;		VEC_MUL_NUM(v,A[dir],3,(d*kx));
            cell.Translate(v);
            //cell.Write(cout,F_XYZ);
            AddStructure(cell);
        }
        int l;
        VEC_MUL_NUM(A[dir],A[dir],3,fabs(no));
        a[dir]*=fabs(no);

        return nAtom;
    }
/*
//private functions
    private void  TransformVec(double[] v){
        int l[3],i,j;
        for(i=0;i<3;i++){
            double d=(v[0]*x[(nAtom+i)*3+0] +v[1]*x[(nAtom+i)*3+1]+v[2]*x[(nAtom+i)*3+2]);
            d/=SQR(a[i]);
            //cout<<endl<<i<<" - "<<d<<" "<<a[i]<<endl;
            if(d>0.5) l[i]=-1;
            else if(d<-0.5) l[i]=1;
                    else l[i]=0;
        }

        for(i=0;i<3;i++)
            v[i]+=l[0]*x[(nAtom+0)*3+i]+l[1]*x[(nAtom+1)*3+i]+l[2]*x[(nAtom+2)*3+i];
    }
*/

    public int ConvertToFract(boolean toFrag){
        CalcAB();
        double[][] temp1=new double[nAtoms][3];
        double[][] temp2=new double[nAtoms][3];

        if(!isFract &&toFrag){ //convert from orthogonal coordinated to fractional ones
            MTools.VEC_TO_MAT2(coords,temp1);
            MTools.MAT2_MUL_MAT2(temp2,temp1,B);
            MTools.MAT2_TO_VEC(temp2,coords);
            isFract=toFrag;
        }

        if(isFract && !toFrag){ //and vice visa (to orthogonal coordinates)
            MTools.VEC_TO_MAT2(coords,temp1);
            MTools.MAT2_MUL_MAT2(temp2,temp1,A);
            MTools.MAT2_TO_VEC(temp2,coords);
            isFract=toFrag;
        }

        return 0;
    }
	
	//================== properties of Clusters

    private double[][] A=new double[3][3]; //lattice vectors
	private double[][] B=new double[3][3]; //inverse matric of A

    private double[] a=new double[3];//lattice distances
    private double[] alpha=new double[3];//lattice angle

    private boolean isCalcAB = false;

    /**
     * Get the value of CalcAB
     *
     * @return the value of CalcAB
     */
    public boolean isCalcAB() {
        return isCalcAB;
    }

    /**
     * Set the value of CalcAB
     *
     * @param CalcAB new value of CalcAB
     */
    public void setCalcAB(boolean CalcAB) {
        this.isCalcAB = CalcAB;
    }

    public int CalcAB(){ //calculate A and B matrice
        A[0][0]=a[0];
        A[0][1]=0;
        A[0][2]=0;

        A[1][0]=a[1]*Math.cos(alpha[2]);
        A[1][1]=a[1]*Math.sin(alpha[2]);
        A[1][2]=0;

        A[2][0]=a[2]*Math.cos(alpha[1]);
        A[2][1]=(a[1]*a[2]*Math.cos(alpha[0])-A[2][0]*A[1][0])/A[1][1];
        A[2][2]=Math.sqrt(MTools.SQR(a[2])-MTools.SQR(A[2][0])-MTools.SQR(A[2][1]));


        MTools.InverseMatrix_apache(A, B);

        /*
        for(i=0;i<3;i++,cout<<endl)
            for(j=0;j<3;j++) cout<<A[i][j]<<" ";

        cout<<" Deter of A = "<<det<<endl;
        for(i=0;i<3;i++,cout<<endl)
            for(j=0;j<3;j++) cout<<B[i][j]<<" ";
        */

        isCalcAB=true;

        return 0;
    }

	private boolean isFract=false; //flag to determine x is fractional coordinate or orthogonal;

    @Override
    protected boolean ReadCAR(Scanner scanner){
       Clear();
       tag="0";
       StringTokenizer tokenizer;

       while(scanner.hasNext()){
            String s=scanner.nextLine();
            if(s.contains("PBC")&&!s.contains("=")){
                tokenizer = new StringTokenizer(s, "\t ,;");
                tokenizer.nextToken();
                for(int i=0;i<3;i++)
                    a[i]=Double.parseDouble(tokenizer.nextToken());
                 for(int i=0;i<3;i++)
                    alpha[i]=Math.toRadians(Double.parseDouble(tokenizer.nextToken()));

                CalcAB();

                Vector coords_tmp=new Vector();
                Vector Nz_tmp=new Vector();

                while(scanner.hasNext()){
                    s=scanner.nextLine();
                    if(s.contains("end")) break;

                    String temp=s.substring(7);
                    tokenizer = new StringTokenizer(temp, "\t ,;"); //read no of atoms
					coords_tmp.add(Double.parseDouble(tokenizer.nextToken()));
					coords_tmp.add(Double.parseDouble(tokenizer.nextToken()));
					coords_tmp.add(Double.parseDouble(tokenizer.nextToken()));

                    temp=s.substring(71,73).trim();
                    Nz_tmp.add(getAtomicNumberFromSymbol(temp));

                }
                setNAtoms(Nz_tmp.size());
				for(int i=0;i<nAtoms;i++){
					Nz[i]=(Integer)Nz_tmp.get(i);
					nType[Nz[i]]++;
					coords[i*3+0]=(Double)coords_tmp.get(i*3+0);
					coords[i*3+1]=(Double)coords_tmp.get(i*3+1);
					coords[i*3+2]=(Double)coords_tmp.get(i*3+2);
				}

                isFract=false;
                return true;
            }
        }
       

        return false;
    }

    @Override
    protected void WriteCAR(Writer writer) throws  IOException{        
		if(nAtoms>0){
            DateFormat dateFormat = new SimpleDateFormat("yyyy/MM/dd HH:mm:ss");
            java.util.Date date = new java.util.Date();
            String s1="!BIOSYM archive 3 \nPBC=ON \n\n!DATE     "+ dateFormat.format(date)+"\n";
            writer.append(s1);

            s1=String.format("PBC   %8.4f  %8.4f   %8.4f  %8.4f  %8.4f  %8.4f\n",a[0],a[1],a[2],
                Math.toDegrees(alpha[0]), Math.toDegrees(alpha[1]), Math.toDegrees(alpha[2]));
            writer.append(s1);

           for(int i=0;i<nAtoms;i++){
                s1=cElements[Nz[i]]+String.format("%-3d   %13.9f  %13.9f  %13.9f XXXX 1      xx      ",
                        (i+1),coords[i*3],coords[i*3+1],coords[i*3+2])
                        + cElements[Nz[i]]+ "   0.000\n";
                writer.append(s1);
           }
            writer.append("end \nend\n");
		}
	}

}
