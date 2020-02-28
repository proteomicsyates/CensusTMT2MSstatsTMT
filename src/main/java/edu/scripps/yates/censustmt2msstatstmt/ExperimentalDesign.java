package edu.scripps.yates.censustmt2msstatstmt;

import java.io.File;
import java.io.IOException;
import java.nio.file.Files;
import java.util.List;
import java.util.Map;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TDoubleSet;
import gnu.trove.set.hash.TDoubleHashSet;

public class ExperimentalDesign {

	private static final String RUN = "Run";
	private static final String FRACTION = "Fraction";
	private static final String TECH_REP_MIXTURE = "TechRepMixture";
	private static final String CHANNEL = "Channel";
	private static final String CONDITION = "Condition";
	private static final String MIXTURE = "Mixture";
	private static final String BIO_REPLICATE = "BioReplicate";
	static final String SYMBOL = "@@@@";
	private final String separator;
	private final Map<String, TDoubleList> channelsByConditions = new THashMap<String, TDoubleList>();
	private final TDoubleSet totalChannels = new TDoubleHashSet();
	private final Map<String, String> mixtureByRun = new THashMap<String, String>();
	private final TObjectIntMap<String> techRepMixtureByRun = new TObjectIntHashMap<String>();
	private final TObjectIntMap<String> bioReplicateByChannelTechRepMixtureAndMixture = new TObjectIntHashMap<String>();

	public ExperimentalDesign(File experimentalDesignFile, String separator) throws IOException {
		this.separator = separator;
		final List<String> lines = Files.readAllLines(experimentalDesignFile.toPath());
		final String header = lines.get(0);
		final TObjectIntMap<String> indexByHeader = getIndexesByHeaders(header);
		for (int numLine = 1; numLine < lines.size(); numLine++) {
			final String line = lines.get(numLine);
			final String[] split = line.split(separator);
			final String run = split[indexByHeader.get(RUN)];
			final int fraction = Integer.valueOf(split[indexByHeader.get(FRACTION)]);
			final int techRepMixture = Integer.valueOf(split[indexByHeader.get(TECH_REP_MIXTURE)]);
			final double channel = Double.valueOf(split[indexByHeader.get(CHANNEL)]);
			final String condition = split[indexByHeader.get(CONDITION)] + SYMBOL + channel;
			final String mixture = split[indexByHeader.get(MIXTURE)];
			final int bioReplicate = Integer.valueOf(split[indexByHeader.get(BIO_REPLICATE)]);
			//
			TDoubleList channels = null;
			if (channelsByConditions.containsKey(condition)) {
				channels = channelsByConditions.get(condition);
			} else {
				channels = new TDoubleArrayList();
			}
			channels.add(channel);
			totalChannels.add(channel);
			channelsByConditions.put(condition, channels);
			mixtureByRun.put(run, mixture);
			techRepMixtureByRun.put(run, techRepMixture);
			bioReplicateByChannelTechRepMixtureAndMixture.put(channel + techRepMixture + mixture, bioReplicate);
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

	public Map<QuantificationLabel, QuantCondition> getConditionByLabel() {
		final Map<QuantificationLabel, QuantCondition> ret = new THashMap<QuantificationLabel, QuantCondition>();
		for (final String condition : this.channelsByConditions.keySet()) {
			final TDoubleList channels = channelsByConditions.get(condition);
			for (final double channel : channels.toArray()) {
				final QuantificationLabel label = getLabelFromChannel(channel);
				ret.put(label, new QuantCondition(condition));
			}
		}
		return ret;
	}

	private QuantificationLabel getLabelFromChannel(double channel) {
		if (isTMT11Plex()) {
			return getTMT11LabelFromChannel(channel);
		} else if (isTMT6Plex()) {
			return getTMT6LabelFromChannel(channel);
		} else {
			throw new IllegalArgumentException("The number of different channels is " + totalChannels.size()
					+ " which is not supported (only 6 or 11 are supported).");
		}
	}

	public static QuantificationLabel getTMT6LabelFromChannel(double channel) {
		if (Double.compare(channel, 126.0) == 0) {
			return QuantificationLabel.TMT_6PLEX_126;
		}
		if (Double.compare(channel, 127.0) == 0) {
			return QuantificationLabel.TMT_6PLEX_127;
		}
		if (Double.compare(channel, 128.0) == 0) {
			return QuantificationLabel.TMT_6PLEX_128;
		}
		if (Double.compare(channel, 129.0) == 0) {
			return QuantificationLabel.TMT_6PLEX_129;
		}
		if (Double.compare(channel, 130.0) == 0) {
			return QuantificationLabel.TMT_6PLEX_130;
		}
		if (Double.compare(channel, 131.0) == 0) {
			return QuantificationLabel.TMT_6PLEX_131;
		}
		throw new IllegalArgumentException("Channel " + channel + " not recognized for TMT6plex");
	}

	public static QuantificationLabel getTMT11LabelFromChannel(double channel) {
		throw new IllegalArgumentException("TMT11 not yet supported");

	}

	public static double getChannelByLabel(QuantificationLabel label) {
		if (label == QuantificationLabel.TMT_6PLEX_126) {
			return 126.0;
		}
		if (label == QuantificationLabel.TMT_6PLEX_127) {
			return 127.0;
		}
		if (label == QuantificationLabel.TMT_6PLEX_128) {
			return 128.0;
		}
		if (label == QuantificationLabel.TMT_6PLEX_129) {
			return 129.0;
		}
		if (label == QuantificationLabel.TMT_6PLEX_130) {
			return 130.0;
		}
		if (label == QuantificationLabel.TMT_6PLEX_131) {
			return 131.0;
		}
		throw new IllegalArgumentException(label + " not supported yet.");
	}

	public String getMixtureByRun(String run) {
		return mixtureByRun.get(run);
	}

	public int getTechRepMixtureByRun(String run) {
		return techRepMixtureByRun.get(run);
	}

	public int getBioReplicate(double channel, int techRepMixture, String mixture) {
		return bioReplicateByChannelTechRepMixtureAndMixture.get(channel + techRepMixture + mixture);
	}

	public boolean isTMT6Plex() {
		return totalChannels.size() == 6;
	}

	public boolean isTMT11Plex() {
		return totalChannels.size() == 11;
	}

}
