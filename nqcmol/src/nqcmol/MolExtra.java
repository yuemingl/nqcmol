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
import nqcmol.potential.*;

/**
 *
 * @author nqc
 */
public class MolExtra {
	static public Potential SetupPotential(String sPotential,String sParam,String sUnit){
		Potential pot=null;
		if(sPotential.contentEquals("LJ")){
			pot=new LennardJonesPotential();
		}

		else if(sPotential.contentEquals("OSS2")){
			pot=new OSS2Potential();
		}

		else if(sPotential.contentEquals("HF_2")){
			pot=new HF2Potential();
		}

		else if(sPotential.contentEquals("TTM21F")){
			pot=new TTM21FPotential();
		}

		assert pot==null;
		pot.setParam(sParam);
		pot.setUnit(sUnit);
		return pot;
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
