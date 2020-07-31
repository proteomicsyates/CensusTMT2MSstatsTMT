package edu.scripps.yates.censustmt2msstatstmt.singletons;

import java.io.File;
import java.io.IOException;
import java.util.Iterator;
import java.util.Map;
import java.util.Set;

import com.compomics.dbtoolkit.io.implementations.FASTADBLoader;
import com.compomics.util.protein.Protein;

import edu.scripps.yates.annotations.uniprot.UniprotProteinLocalRetriever;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.proteomicsmodel.Accession;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class ProteinSequences extends THashMap<String, String> implements Clearable {
	private static ProteinSequences instance;

	public static ProteinSequences getInstance() {
		if (instance == null) {
			instance = new ProteinSequences();
		}
		return instance;
	}

	private ProteinSequences() {
		super();
	}

	/**
	 * It will look for protein sequences of the input acc set. If fasta is
	 * provided, it will look into the fasta file. If not, or not found in the
	 * fasta, it will go to uniprot in case of being uniprot accessions.
	 * 
	 * @param accs accessions of the proteins
	 * @param File fasta file. If not null, it will be used preferentially instead
	 *             of {@link UniprotProteinLocalRetriever} to get the protein
	 *             sequences
	 * @param uplr This will use UniprotKB to get the protein sequences from the
	 *             proteins that are not found in fasta file or when fasta file is
	 *             not provided.
	 * @return
	 * @throws IOException
	 */
	public void getProteinSequences(Set<String> accs, File fastaFile) throws IOException {
		final Set<String> missedAccs = new THashSet<String>();
		if (fastaFile != null) {
			final Map<String, String> sequencesFromFasta = new THashMap<String, String>();
			final FASTADBLoader fastaLoader = new FASTADBLoader();
			if (!fastaLoader.canReadFile(fastaFile)) {
				throw new IllegalArgumentException("Error reading fasta file!");
			}
			System.out.println("Reading fasta file '" + fastaFile.getAbsolutePath() + "'...");

			fastaLoader.load(fastaFile.getAbsolutePath());
			Protein protein = null;
			int numSequences = 0;
			while ((protein = fastaLoader.nextProtein()) != null) {
				final String header = protein.getHeader().getRawHeader();
				final Accession acc = FastaParser.getACC(header);
				final String sequence = protein.getSequence().getSequence();
				sequencesFromFasta.put(acc.getAccession(), sequence);
				numSequences++;
				final String gene = FastaParser.getGeneFromFastaHeader(header);
				if (gene != null) {
					GeneMapper.getInstance().addMapping(acc.getAccession(), gene);
				}
			}
			fastaLoader.close();
			System.out.println(numSequences + " protein sequences read from Fasta file");
			// we look if we missed some accs
			for (final String acc : accs) {
				if (!sequencesFromFasta.containsKey(acc)) {
					missedAccs.add(acc);
				} else {
					put(acc, sequencesFromFasta.get(acc));
				}
			}
		} else {
			accs.stream().forEach(acc -> {
				if (!this.containsKey(acc)) {
					missedAccs.add(acc);
				}
			});
		}
		if (!missedAccs.isEmpty()) {
			// we go to uniprot to get all missed ones
			final Map<String, Entry> annotatedProteins = UPLR.getInstance().getAnnotatedProteins(null, missedAccs);
			final Iterator<String> iterator = missedAccs.iterator();
			while (iterator.hasNext()) {
				final String acc = iterator.next();
				if (annotatedProteins.containsKey(acc)) {
					final Entry entry = annotatedProteins.get(acc);
					final String proteinSequence = UniprotEntryUtil.getProteinSequence(entry);
					if (proteinSequence != null) {
						put(acc, proteinSequence);
						iterator.remove();
					}
				}
			}
		}
		if (!missedAccs.isEmpty()) {
			System.out.println("Protein sequences for " + missedAccs.size() + " proteins were not found");
			final File outputFile = FilesManager.getInstance().printList(missedAccs, "_no_protein_sequences");
			System.out.println("Proteins with no protein sequences written to " + outputFile);
		}

	}

	@Override
	public void clearData() {
		this.clear();
	}
}
