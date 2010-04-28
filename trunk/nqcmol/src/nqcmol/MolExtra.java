/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

//import com.sun.xml.internal.txw2.output.XMLWriter;
//import java.io.StringReader;
//import java.io.StringWriter;
//import javax.xml.transform.OutputKeys;
//import javax.xml.transform.Source;
//import javax.xml.transform.Transformer;
//import javax.xml.transform.TransformerFactory;
//import javax.xml.transform.stream.StreamResult;
//import javax.xml.transform.stream.StreamSource;
import java.util.Vector;
import nqcmol.potential.*;

/**
 *
 * @author nqc
 */
public class MolExtra {
	static public Potential SetupPotential(String sPotential,String sParam,String sUnit,String sMethod){
		Potential pot=null;
		if(sPotential.contentEquals("LJ")){
			pot=new LennardJonesPotential();
		}

		else if(sPotential.contentEquals("OSS2")){
			pot=new OSS2Potential();
		}

                else if(sPotential.contentEquals("FOSS2")){
			pot=new FOSS2Potential();
		}

		else if(sPotential.contentEquals("HF2")){
			pot=new HF2Potential();
		}

		else if(sPotential.contentEquals("TTM21F")){
			pot=new TTM21FPotential();
		}

        else if(sPotential.contentEquals("KJ")){
			pot=new KJPotential();
		}

		else if(sPotential.contentEquals("g03")){
			pot=new GaussianInterfacePotential();
			GaussianInterfacePotential gauss=(GaussianInterfacePotential)pot;
			if(!sParam.isEmpty()){
				String[] tmp=sParam.split("@");
				if(tmp.length>=1)
					gauss.setBasisSet(tmp[0]);
				if(tmp.length>=3)	gauss.setChargeAndMultiplicity(tmp[1]+" "+tmp[2]);
			}
		}

		assert pot==null;
		pot.setParam(sParam);
		pot.setUnit(sUnit);
		pot.setOptMethod(sMethod);
		return pot;
	}


	public static double[] ConvertToDoubleArray(Vector<Double> vec){
		double[] answer=null;
		if(vec.size()>0){
			answer=new double[vec.size()];
			for(int i=0;i<vec.size();i++)
				answer[i]=vec.get(i);
		}
		return answer;
	}


	/**
	 *
	 * @param m input cluster
	 * @param ref
	 * @return
	 */
	static public int indexForBE(Cluster m,Vector<Cluster> ref){
		int nO=m.getNonHydrogenNum();
		int nH=m.getHydrogenNum();

		int remain=nH-ref.get(0).getHydrogenNum()*(nO-1);		//cout<<" Charge =" <<data[i].GetTotalCharge()<<" Mass = "<<data[i].GetMolWt()<<" numNonH= "<<nF<<endl;

		for(int k=0;k<ref.size();k++)
			if(remain==ref.get(k).getHydrogenNum()){
				 return k;
			}
		return 0;
	}

	static public int getIndexOfArguments(String arg0,String[] args){
		for (int i = 0; i < args.length; i++) {
				if (args[i].contentEquals(arg0))
					return i;
			}
		return -1;
	}

//	public static String prettyFormat(String input, int indent) {
//    try {
//        Source xmlInput = new StreamSource(new StringReader(input));
//        StringWriter stringWriter = new StringWriter();
//        StreamResult xmlOutput = new StreamResult(stringWriter);
//        Transformer transformer = TransformerFactory.newInstance().newTransformer();
//        transformer.setOutputProperty(OutputKeys.INDENT, "yes");
//        transformer.setOutputProperty("{http://xml.apache.org/xslt}indent-amount", String.valueOf(indent));
//        transformer.transform(xmlInput, xmlOutput);
//        return xmlOutput.getWriter().toString();
//		} catch (Exception e) {
//			throw new RuntimeException(e); // simple exception handling, please review it
//		}
//	}
//
//	public static String prettyFormat(String input) {
//		return prettyFormat(input, 2);
//	}


}
