/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */

package nqcmol;

import java.util.logging.*;
import org.kohsuke.args4j.*;

/**
 *
 * @author nqc
 */
public class MClassify {
	@Option(name="-h",usage="Print out the help")
    boolean isHelp= false;


	CmdLineParser parser=null;

	public void ParseArguments(String[] args){
		parser = new CmdLineParser(this);
		try {
			parser.parseArgument(args);
		} catch (CmdLineException ex) {
			Logger.getLogger(MCalculate.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	public void ClassifyWaterClusters(String[] args){
		
		
	}
	

}
