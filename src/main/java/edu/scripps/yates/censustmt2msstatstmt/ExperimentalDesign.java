package edu.scripps.yates.censustmt2msstatstmt;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.utilities.files.FileUtils;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.hash.THashSet;

public class ExperimentalDesign {

	private static final String RUN = "Run";
	private static final String FRACTION = "Fraction";
	private static final String TECH_REP_MIXTURE = "TechRepMixture";
	private static final String CHANNEL = "Channel";
	private static final String CONDITION = "Condition";
	private static final String MIXTURE = "Mixture";
	private static final String BIO_REPLICATE = "BioReplicate";
//	static final String SYMBOL = "@@@@";
	private final String separator;

	private final Map<String, Mixture> mixtureByMixtureName = new THashMap<String, Mixture>();
	private final Map<String, Mixture> mixtureByRun = new THashMap<String, Mixture>();
	private final Map<String, String> techRepMixtureByRun = new THashMap<String, String>();
	private final Map<String, String> bioReplicateByChannelTechRepMixtureAndMixture = new THashMap<String, String>();
	private final Map<String, Set<String>> fractionsByRun = new THashMap<String, Set<String>>();
	private int tmtPlex;

	public ExperimentalDesign(File experimentalDesignFile, String separator) throws IOException {
		try {
			this.separator = separator;
			final List<String> lines = FileUtils.readAllLines(experimentalDesignFile.toPath());
			final String header = lines.get(0);
			final TObjectIntMap<String> indexByHeader = getIndexesByHeaders(header);
			for (int numLine = 1; numLine < lines.size(); numLine++) {
				final String line = lines.get(numLine).trim();
				if ("".equals(line)) {
					continue;
				}
				final String[] split = line.split(separator);
				final String run = split[indexByHeader.get(RUN)].trim();
				final String fraction = split[indexByHeader.get(FRACTION)].trim();
				final String techRepMixture = split[indexByHeader.get(TECH_REP_MIXTURE)].trim();
				Float channel = Float.NaN;
				try {
					channel = Float.valueOf(split[indexByHeader.get(CHANNEL)].trim());
				} catch (final NumberFormatException e) {
					throw new Exception("Channel '" + split[indexByHeader.get(CHANNEL)].trim()
							+ "' is not supported. Only numbers are allowed.");
				}
				final String condition = split[indexByHeader.get(CONDITION)].trim();
				final String mixtureName = split[indexByHeader.get(MIXTURE)].trim();
				Mixture mixture = null;
				if (mixtureByMixtureName.containsKey(mixtureName)) {
					mixture = mixtureByMixtureName.get(mixtureName);
				} else {
					mixture = new Mixture(mixtureName);
					mixtureByMixtureName.put(mixtureName, mixture);
				}
				final String bioReplicate = split[indexByHeader.get(BIO_REPLICATE)].trim();
				//

				mixture.addChannelWithCondition(channel, condition);
				mixture.addRun(run);
				mixtureByRun.put(run, mixture);
				if (!fractionsByRun.containsKey(run)) {
					fractionsByRun.put(run, new THashSet<String>());
				}
				fractionsByRun.get(run).add(fraction);
				techRepMixtureByRun.put(run, techRepMixture);
				bioReplicateByChannelTechRepMixtureAndMixture.put(channel + techRepMixture + mixture.getName(),
						bioReplicate);
			}
		} catch (final Exception e) {
			if (e instanceof IOException) {
				throw (IOException) e;
			}
			throw new IllegalArgumentException(
					"Error reading experimental design from the annotation file: " + e.getMessage());
		}
	}

	private TObjectIntMap<String> getIndexesByHeaders(String header) {
		final TObjectIntMap<String> ret = new TObjectIntHashMap<String>();
		final String[] split = header.split(separator);
		for (int index = 0; index < split.length; index++) {
			ret.put(split[index], index);
		}
		return ret;
	}

	public Set<String> getFractionsByRun(String run) {
		return fractionsByRun.get(run);
	}

	public Mixture getMixtureByRun(String run) {
		return mixtureByRun.get(run);
	}

	public String getTechRepMixtureByRun(String run) {
		return techRepMixtureByRun.get(run);
	}

	public String getBioReplicate(float channel, String techRepMixture, String mixture) {
		return bioReplicateByChannelTechRepMixtureAndMixture.get(channel + techRepMixture + mixture);
	}

	public Collection<Mixture> getMixtures() {
		return mixtureByMixtureName.values();
	}

	/**
	 * Sets the tmtPlex to all mixtures in this experimental design.
	 * 
	 * @param tmtPlex
	 */
	public void setTMTPlex(int tmtPlex) {

		for (final Mixture mixture : getMixtures()) {
			mixture.setTMTPlex(tmtPlex);
		}
	}

}
