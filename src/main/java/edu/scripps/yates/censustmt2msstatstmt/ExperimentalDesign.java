package edu.scripps.yates.censustmt2msstatstmt;

import java.io.File;
import java.io.IOException;
import java.util.List;
import java.util.Map;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.files.FileUtils;
import gnu.trove.list.TFloatList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.TObjectIntMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectIntHashMap;
import gnu.trove.set.TFloatSet;
import gnu.trove.set.hash.TFloatHashSet;

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
	private final Map<String, TFloatList> channelsByConditions = new THashMap<String, TFloatList>();
	private final TFloatSet totalChannels = new TFloatHashSet();
	private final Map<String, String> mixtureByRun = new THashMap<String, String>();
	private final Map<String, String> techRepMixtureByRun = new THashMap<String, String>();
	private final Map<String, String> bioReplicateByChannelTechRepMixtureAndMixture = new THashMap<String, String>();
	private final Map<String, QuantCondition> conditionsByName = new THashMap<String, QuantCondition>();
	private final Map<QuantificationLabel, Float> channelsByLabel = new THashMap<QuantificationLabel, Float>();

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
				final String run = split[indexByHeader.get(RUN)];
				final String fraction = split[indexByHeader.get(FRACTION)];
				final String techRepMixture = split[indexByHeader.get(TECH_REP_MIXTURE)];
				Float channel = Float.NaN;
				try {
					channel = Float.valueOf(split[indexByHeader.get(CHANNEL)].trim());
				} catch (final NumberFormatException e) {
					throw new Exception("Channel '" + split[indexByHeader.get(CHANNEL)]
							+ "' is not supported. Only numbers are allowed.");
				}
				final String condition = split[indexByHeader.get(CONDITION)] + SYMBOL + channel;
				final String mixture = split[indexByHeader.get(MIXTURE)];
				final String bioReplicate = split[indexByHeader.get(BIO_REPLICATE)];
				//
				TFloatList channels = null;
				if (channelsByConditions.containsKey(condition)) {
					channels = channelsByConditions.get(condition);
				} else {
					channels = new TFloatArrayList();
				}
				channels.add(channel);
				totalChannels.add(channel);
				channelsByConditions.put(condition, channels);
				mixtureByRun.put(run, mixture);
				techRepMixtureByRun.put(run, techRepMixture);
				bioReplicateByChannelTechRepMixtureAndMixture.put(channel + techRepMixture + mixture, bioReplicate);
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

	public Map<QuantificationLabel, QuantCondition> getConditionsByLabel() {
		final Map<QuantificationLabel, QuantCondition> ret = new THashMap<QuantificationLabel, QuantCondition>();
		for (final String condition : this.channelsByConditions.keySet()) {
			QuantCondition quantCondition = null;
			if (conditionsByName.containsKey(condition)) {
				quantCondition = conditionsByName.get(condition);
			} else {
				quantCondition = new QuantCondition(condition);
			}
			final TFloatList channels = channelsByConditions.get(condition);
			for (final float channel : channels.toArray()) {
				final QuantificationLabel label = getLabelFromChannel(channel);
				ret.put(label, quantCondition);
			}
		}
		return ret;
	}

	public Map<QuantCondition, QuantificationLabel> getLabelsByCondition() {
		final Map<QuantCondition, QuantificationLabel> ret = new THashMap<QuantCondition, QuantificationLabel>();
		for (final String condition : this.channelsByConditions.keySet()) {
			QuantCondition quantCondition = null;
			if (conditionsByName.containsKey(condition)) {
				quantCondition = conditionsByName.get(condition);
			} else {
				quantCondition = new QuantCondition(condition);
			}
			final TFloatList channels = channelsByConditions.get(condition);
			for (final float channel : channels.toArray()) {
				final QuantificationLabel label = getLabelFromChannel(channel);
				ret.put(quantCondition, label);
			}
		}
		return ret;
	}

	public QuantificationLabel getLabelFromChannel(float channel) {
		if (isTMT11Plex()) {
			return getTMT11LabelFromChannel(channel);
		} else if (isTMT10Plex()) {
			return getTMT10LabelFromChannel(channel);
		} else if (isTMT6Plex()) {
			return getTMT6LabelFromChannel(channel);
		} else {
			throw new IllegalArgumentException("The number of different channels is " + totalChannels.size()
					+ " which is not supported (only 6 or 11 are supported).");
		}
	}

	private QuantificationLabel getTMT6LabelFromChannel(float channel) {
		final TFloatList sortedChannels = new TFloatArrayList();
		sortedChannels.addAll(this.totalChannels);
		sortedChannels.sort();
		//
		final List<QuantificationLabel> tmt6PlexLabels = QuantificationLabel.getTMT6PlexLabels();
		for (int i = 0; i < tmt6PlexLabels.size(); i++) {
			if (Float.compare(sortedChannels.get(i), channel) == 0) {
				final QuantificationLabel quantificationLabel = tmt6PlexLabels.get(i);
				channelsByLabel.put(quantificationLabel, channel);
				return quantificationLabel;
			}
		}

		throw new IllegalArgumentException("Channel " + channel + " not recognized for TMT6plex");
	}

	private QuantificationLabel getTMT10LabelFromChannel(float channel) {
		final TFloatList sortedChannels = new TFloatArrayList();
		sortedChannels.addAll(this.totalChannels);
		sortedChannels.sort();
		//
		final List<QuantificationLabel> tmt10PlexLabels = QuantificationLabel.getTMT10PlexLabels();
		for (int i = 0; i < tmt10PlexLabels.size(); i++) {

			if (Float.compare(sortedChannels.get(i), channel) == 0) {
				final QuantificationLabel quantificationLabel = tmt10PlexLabels.get(i);
				channelsByLabel.put(quantificationLabel, channel);
				return quantificationLabel;
			}
		}

		throw new IllegalArgumentException("Channel " + channel + " not recognized for TMT6plex");
	}

	private QuantificationLabel getTMT11LabelFromChannel(float channel) {
		final TFloatList sortedChannels = new TFloatArrayList();
		sortedChannels.addAll(this.totalChannels);
		sortedChannels.sort();
		//
		final List<QuantificationLabel> tmt11PlexLabels = QuantificationLabel.getTMT11PlexLabels();
		for (int i = 0; i < tmt11PlexLabels.size(); i++) {
			if (Float.compare(sortedChannels.get(i), channel) == 0) {
				final QuantificationLabel quantificationLabel = tmt11PlexLabels.get(i);
				channelsByLabel.put(quantificationLabel, channel);
				return quantificationLabel;
			}
		}

		throw new IllegalArgumentException("Channel " + channel + " not recognized for TMT6plex");
	}

	public Float getChannelByLabel(QuantificationLabel label) {
		if (channelsByLabel.containsKey(label)) {
			return channelsByLabel.get(label);
		}
		throw new IllegalArgumentException(label + " not supported yet.");
	}

	public String getMixtureByRun(String run) {
		return mixtureByRun.get(run);
	}

	public String getTechRepMixtureByRun(String run) {
		return techRepMixtureByRun.get(run);
	}

	public String getBioReplicate(float channel, String techRepMixture, String mixture) {
		return bioReplicateByChannelTechRepMixtureAndMixture.get(channel + techRepMixture + mixture);
	}

	private boolean isTMT6Plex() {
		return totalChannels.size() == 6;
	}

	private boolean isTMT10Plex() {
		return totalChannels.size() == 10;
	}

	private boolean isTMT11Plex() {
		return totalChannels.size() == 11;
	}

}
