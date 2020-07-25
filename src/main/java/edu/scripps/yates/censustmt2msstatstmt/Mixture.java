package edu.scripps.yates.censustmt2msstatstmt;

import java.util.List;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import gnu.trove.list.TFloatList;
import gnu.trove.list.array.TFloatArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class Mixture {
	private final Map<Float, String> conditionsByChannels = new THashMap<Float, String>();
	private final String name;
	private final Map<String, QuantCondition> conditionsByName = new THashMap<String, QuantCondition>();
	private final Map<QuantificationLabel, Float> channelsByLabel = new THashMap<QuantificationLabel, Float>();
	private final Set<String> runs = new THashSet<String>();

	public Map<QuantificationLabel, QuantCondition> getConditionsByLabels() {
		final Map<QuantificationLabel, QuantCondition> ret = new THashMap<QuantificationLabel, QuantCondition>();

		for (final Float channel : getConditionsByChannels().keySet()) {
			final QuantificationLabel label = getLabelFromChannel(channel);
			final String condition = getConditionsByChannels().get(channel);
			QuantCondition quantCondition = null;
			if (conditionsByName.containsKey(condition)) {
				quantCondition = conditionsByName.get(condition);
			} else {
				quantCondition = new QuantCondition(condition);
			}

			ret.put(label, quantCondition);

		}
		return ret;
	}

	public Mixture(String mixtureName) {
		this.name = mixtureName;
	}

	public Map<Float, String> getConditionsByChannels() {
		return this.conditionsByChannels;
	}

	public String getName() {
		return this.name;
	}

	public QuantificationLabel getLabelFromChannel(float channel) {
		if (isTMT11Plex()) {
			return getTMT11LabelFromChannel(channel);
		} else if (isTMT10Plex()) {
			return getTMT10LabelFromChannel(channel);
		} else if (isTMT6Plex()) {
			return getTMT6LabelFromChannel(channel);
		} else if (isTMT4Plex()) {
			return getTMT4LabelFromChannel(channel);
		} else {
			throw new IllegalArgumentException("The number of different channels is " + getChannels().size()
					+ " which is not supported (only 6 or 11 are supported).");
		}
	}

	private QuantificationLabel getTMT6LabelFromChannel(float channel) {
		final TFloatList sortedChannels = new TFloatArrayList();
		sortedChannels.addAll(this.getChannels());
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

	private QuantificationLabel getTMT4LabelFromChannel(float channel) {
		final TFloatList sortedChannels = new TFloatArrayList();
		sortedChannels.addAll(this.getChannels());
		sortedChannels.sort();
		//
		final List<QuantificationLabel> tmt4PlexLabels = QuantificationLabel.getTMT4PlexLabels();
		for (int i = 0; i < tmt4PlexLabels.size(); i++) {
			if (Float.compare(sortedChannels.get(i), channel) == 0) {
				final QuantificationLabel quantificationLabel = tmt4PlexLabels.get(i);
				channelsByLabel.put(quantificationLabel, channel);
				return quantificationLabel;
			}
		}

		throw new IllegalArgumentException("Channel " + channel + " not recognized for TMT4plex");
	}

	private QuantificationLabel getTMT10LabelFromChannel(float channel) {
		final TFloatList sortedChannels = new TFloatArrayList();
		sortedChannels.addAll(this.getChannels());
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
		sortedChannels.addAll(this.getChannels());
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
		getConditionsByLabels();
		if (channelsByLabel.containsKey(label)) {
			return channelsByLabel.get(label);
		}
		throw new IllegalArgumentException(label + " not supported yet.");
	}

	private boolean isTMT4Plex() {
		return getChannels().size() == 4;
	}

	private boolean isTMT6Plex() {
		return getChannels().size() == 6;
	}

	private boolean isTMT10Plex() {
		return getChannels().size() == 10;
	}

	private boolean isTMT11Plex() {
		return getChannels().size() == 11;
	}

	public void addChannelWithCondition(float channel, String condition) {
		this.conditionsByChannels.put(channel, condition);

	}

	private Set<Float> getChannels() {
		return conditionsByChannels.keySet();
	}

	public void addRun(String run) {
		this.runs.add(run);
	}

	public Set<String> getRuns() {
		return runs;
	}
}
