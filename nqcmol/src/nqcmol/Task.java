/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nqcmol;

import java.io.*;
import java.util.Scanner;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.tools.XmlWriter;
import org.kohsuke.args4j.*;

/**
 *
 * @author nqc
 */
public class Task {

	@Option(name = "-i", usage = "Input file name. [xyz]", metaVar = "FILE")
	String sFileIn = "";
	String sFormatIn = "xyz";
	@Option(name = "-o", usage = "Output file name. [xyz]", metaVar = "FILE")
	String sFileOut = "";
	String sFormatOut = "xyz";
	
	@Option(name = "-v", usage = "Verbose level. [0]", metaVar = "INTEGER")
	int verbose = 0;
	
	@Option(name = "-h", usage = "Print out the help")
	boolean isHelp = false;
	CmdLineParser parser = null;
	BufferedWriter stdwriter = new BufferedWriter(new OutputStreamWriter(System.out));
	XmlWriter xmllog = new XmlWriter(stdwriter);

	Scanner fileIn=null;
	FileWriter fileOut=null;

	String[] args=null;

	public void ParseArguments(String[] args) {
		this.args=args;
		parser = new CmdLineParser(this);
		try {
			parser.parseArgument(args);

			if(isHelp){
				//parser.printUsage(System.out);
				return ;
			}

			String input = "";
			String output = "";
			for (int i = 0; i < args.length-1; i++) {
				if(args[i].length()>=2){
					if (args[i].substring(0, 2).contentEquals("-i")) {
						input = args[i];
						sFileIn = args[i + 1];
					}
					if (args[i].substring(0, 2).contentEquals("-o")) {
						output = args[i];
						sFileOut = args[i + 1];
					}
				}
			}
			//System.out.println(input+" and "+output);
			for (int i = 0; i < Cluster.format.length; i++) {
				String format = "-i" + Cluster.format[i];
				if (format.contentEquals(input)) {
					sFormatIn = Cluster.format[i];
				}
				format = "-o" + Cluster.format[i];
				if (format.contentEquals(output)) {
					sFormatOut = Cluster.format[i];
				}
				//System.out.println(format+" "+sFormatIn+" and "+sFormatOut);
			}

			fileIn = new Scanner(new File(sFileIn));
			if(!sFileOut.isEmpty()) fileOut=new FileWriter(new File(sFileOut));

		} catch (IOException ex) {
			Logger.getLogger(Task.class.getName()).log(Level.SEVERE, null, ex);
		} catch (CmdLineException ex) {
			Logger.getLogger(Task.class.getName()).log(Level.SEVERE, null, ex);
		}

	}

	public void Execute(String[] args) {
		ParseArguments(args);
		if (isHelp) {
			System.out.println(" nqcmol - utilities for processing data. Author: Nguyen Quoc Chinh\n");
			//System.out.println(" USAGE\n\t\t nqcmol [OPTION]\n");
			System.out.println(" OPTION\n");
			parser.printUsage(System.out);
			return;
		}else{
			try {
				xmllog.writeEntity(getName());
				Initialize();
				Process();
				Finalize();
				xmllog.endEntity();
				xmllog.close();
				stdwriter.close();
			} catch (IOException ex) {
				Logger.getLogger(Task.class.getName()).log(Level.SEVERE, null, ex);
			}
		}
	}

	public String getName() {
		return "None";
	}

	protected void Initialize() {				
		xmllog.writeAttribute("InputFile", sFileIn).writeAttribute("FormatInput", sFormatIn);
		xmllog.writeAttribute("OutputFile", sFileOut).writeAttribute("FormatOutput", sFormatOut);
		xmllog.writeAttribute("Verbose", Integer.toString(verbose));
	}

	protected void Process() {
		throw new UnsupportedOperationException("Not yet implemented");
	}

	protected void Finalize() {
		try {
			if(fileIn!=null) fileIn.close();
			if(fileOut!=null) fileOut.close();
		} catch (IOException ex) {
			Logger.getLogger(Task.class.getName()).log(Level.SEVERE, null, ex);
		}
	}
}
