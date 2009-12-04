/*
 * To change this template, choose Tools | Templates
 * and open the template in the editor.
 */
package nqcmol;

import java.io.*;
import java.util.logging.Level;
import java.util.logging.Logger;
import nqcmol.tools.XmlWriter;
import org.kohsuke.args4j.*;

/**
 *
 * @author nqc
 */
public class Task {

	@Option(name = "-i", usage = "input file name", metaVar = "FILE")
	String sFileIn = "";
	String sFormatIn = "xyz";
	@Option(name = "-o", usage = "output file name", metaVar = "FILE")
	String sFileOut = "";
	String sFormatOut = "xyz";
	@Option(name = "-v", usage = "verbose level, default is 0", metaVar = "INTEGER")
	int verbose = 0;
	@Option(name = "-h", usage = "Print out the help")
	boolean isHelp = false;
	CmdLineParser parser = null;
	BufferedWriter stdwriter = new BufferedWriter(new OutputStreamWriter(System.out));
	XmlWriter xmllog = new XmlWriter(stdwriter);

	public void ParseArguments(String[] args) {
		parser = new CmdLineParser(this);
		try {
			parser.parseArgument(args);
			String input = "";
			String output = "";
			for (int i = 0; i < args.length; i++) {

				if (args[i].substring(0, 2).contentEquals("-i")) {
					input = args[i];
					sFileIn = args[i + 1];
				}
				if (args[i].substring(0, 2).contentEquals("-o")) {
					output = args[i];
					sFileOut = args[i + 1];
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

		} catch (CmdLineException ex) {
			Logger.getLogger(Task.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	public void Execute(String[] args) {
		ParseArguments(args);
		if (isHelp) {
			parser.printUsage(System.out);
			return;
		}

		Initialize();
		Process();
		Finalize();
	}

	public String getName() {
		return "None";
	}

	protected void Initialize() {
		try {
			xmllog.writeEntity(getName());
			xmllog.writeAttribute("InputFile", sFileIn).writeAttribute("FormatInput", sFormatIn);
			xmllog.writeAttribute("OutputFile", sFileOut).writeAttribute("FormatOutput", sFormatOut);
			xmllog.writeAttribute("Verbose", Integer.toString(verbose));
		} catch (IOException ex) {
			Logger.getLogger(Task.class.getName()).log(Level.SEVERE, null, ex);
		}
	}

	protected void Process() {
		throw new UnsupportedOperationException("Not yet implemented");
	}

	protected void Finalize() {
		try {
			xmllog.endEntity();
			xmllog.close();
			stdwriter.close();
		} catch (IOException ex) {
			Logger.getLogger(Task.class.getName()).log(Level.SEVERE, null, ex);
		}
	}
}
