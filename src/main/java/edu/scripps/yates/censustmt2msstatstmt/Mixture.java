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
	private Integer tmtPlex;

	public Map<QuantificationLabel, QuantCondition> getConditionsByLabels() {
		final Map<QuantificationLabel, QuantCondition> ret = new THashMap<QuantificationLabel, QuantCondition>();

		for (final Float channel : getConditionsByChannels().keySet()) {
			final QuantificationLabel label = getTMTLabelFromChannel(channel);
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

	private QuantificationLabel getTMTLabelFromChannel(float channel) {
		if (tmtPlex == null) {
			throw new IllegalArgumentException(
					"Internal error: We need to set the TMT Plex before calling getLabelFromChannel");
		}

		final TFloatList sortedChannels = new TFloatArrayList();
		sortedChannels.addAll(this.getChannels());
		sortedChannels.sort();
		//

		final List<QuantificationLabel> tmtPlexLabels = QuantificationLabel.getTMTPlexLabels(this.tmtPlex);
		for (int i = 0; i < tmtPlexLabels.size(); i++) {

			if (Float.compare(sortedChannels.get(i), channel) == 0) {
				final QuantificationLabel quantificationLabel = tmtPlexLabels.get(i);
				channelsByLabel.put(quantificationLabel, channel);
				return quantificationLabel;
			}
		}

		throw new IllegalArgumentException("Channel " + channel
				+ " not recognized as a supported channel for TMT. Contact salvador@scripps.edu if you want to support a new type of TMT labeling");

	}

	public Float getChannelByLabel(QuantificationLabel label) {
		getConditionsByLabels();
		if (channelsByLabel.containsKey(label)) {
			return channelsByLabel.get(label);
		}
		return null;
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

	public void setTMTPlex(int tmtPlex) {
		this.tmtPlex = tmtPlex;
	}
}
