package edu.scripps.yates.censustmt2msstatstmt;

import javax.swing.JFrame;

import org.apache.commons.cli.Options;

import edu.scripps.yates.utilities.swing.AutomaticGUICreator;

public class CensusTMT2MSstatsTMTGui {

	public static void main(String[] args) {
		final Options options = CensusTMT2MSstatsTMT.setupCommandLineOptions();
		final JFrame frame = new AutomaticGUICreator(new CensusTMT2MSstatsTMT(options));
		frame.setVisible(true);
	}

}
