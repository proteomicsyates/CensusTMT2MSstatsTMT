package edu.scripps.yates.censustmt2msstatstmt.singletons;

import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;

/**
 * This class will keep a acc to gene mappping. It will store the mapping from
 * reading the fasta file and then, it will perform a lookup into the
 * UniprotGeneMapping if not present.
 * 
 * @author salvador
 *
 */
public class GeneMapper implements Clearable {
	private final Map<String, String> genesByAcc = new THashMap<String, String>();
	private final UniprotProteinLocalRetriever uplr;
	private static GeneMapper instance;

	public static GeneMapper getInstance() {
		if (instance == null) {
			instance = new GeneMapper(UPLR.getInstance());
		}
		return instance;
	}

	private GeneMapper(UniprotProteinLocalRetriever uplr) {
		this.uplr = uplr;
	}

	public void addMapping(String acc, String gene) {
		this.genesByAcc.put(acc, gene);
	}

	public String getGene(String acc) throws IOException {
		if (!genesByAcc.containsKey(acc)) {
			final Map<String, Entry> annotatedProtein = uplr.getAnnotatedProtein(null, acc);
			if (annotatedProtein.containsKey(acc)) {
				final Entry entry = annotatedProtein.get(acc);
				final List<Pair<String, String>> geneName = UniprotEntryUtil.getGeneName(entry, true, true);
				if (geneName != null && !geneName.isEmpty()) {
					final String gene = geneName.get(0).getFirstelement();
					genesByAcc.put(acc, gene);
				}
			}

		}
		return genesByAcc.get(acc);
	}

	@Override
	public void clearData() {
		genesByAcc.clear();
	}

	public void mapGenesFromProteinDescriptionsInInputFile(Collection<QuantifiedProteinInterface> proteins) {
		for (final QuantifiedProteinInterface protein : proteins) {
			final String gene = FastaParser.getGeneFromFastaHeader(protein.getDescription());
			if (gene != null) {
				addMapping(protein.getAccession(), gene);
			}
		}
	}
}
