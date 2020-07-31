package edu.scripps.yates.censustmt2msstatstmt.singletons;

import java.io.File;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;

public class UPLR extends UniprotProteinLocalRetriever {
	private static UPLR instance;

	public static UPLR getInstance() {
		if (instance == null) {
			final File runningFolder = FilesManager.getRunningFolder();
			instance = new UPLR(runningFolder);
		}
		return instance;
	}

	private UPLR(File uniprotReleasesFolder) {
		super(uniprotReleasesFolder, true, true, true);
		this.setRetrieveFastaIsoforms(false);
		this.setRetrieveFastaIsoformsFromMainForms(false);
		this.setCacheEnabled(true);
	}

}
