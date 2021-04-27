package edu.scripps.yates.censustmt2msstatstmt.util;

import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.EnumMap;
import java.util.List;
import java.util.Map;
import java.util.Optional;
import java.util.Set;
import java.util.stream.Collectors;

import edu.scripps.yates.census.read.CensusOutParser;
import edu.scripps.yates.census.read.model.QuantAmount;
import edu.scripps.yates.census.read.model.interfaces.QuantParser;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.censustmt2msstatstmt.PSMAggregationType;
import edu.scripps.yates.censustmt2msstatstmt.singletons.GeneMapper;
import edu.scripps.yates.censustmt2msstatstmt.singletons.PTMList;
import edu.scripps.yates.censustmt2msstatstmt.singletons.UPLR;
import edu.scripps.yates.utilities.annotations.uniprot.UniprotEntryUtil;
import edu.scripps.yates.utilities.annotations.uniprot.xml.Entry;
import edu.scripps.yates.utilities.fasta.FastaParser;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.ProteinEvidence;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.HasProteins;
import edu.scripps.yates.utilities.proteomicsmodel.HasScores;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class DataUtil {
	private static final double MASS_PRECISION = 0.01;
	public static final String SEPARATOR = "|";

	public static String getIntensitiesString(QuantifiedPeptideInterface peptide, PSMAggregationType aggregation,
			List<QuantificationLabel> labels) {
		final StringBuilder sb = new StringBuilder();
		final Set<Amount> amounts = new THashSet<Amount>();
		peptide.getQuantifiedPSMs().stream().filter(p -> !p.isDiscarded()).map(p -> p.getAmounts())
				.forEach(am -> amounts.addAll(am));

		final Map<QuantificationLabel, TDoubleList> map = new THashMap<QuantificationLabel, TDoubleList>();
		for (final QuantificationLabel label : labels) {

			for (final Amount amount : amounts) {
				final QuantAmount quantAmount = (QuantAmount) amount;
				if (quantAmount.getAmountType() == AmountType.INTENSITY) {
					if (quantAmount.getLabel() == label) {
						if (!map.containsKey(label)) {
							map.put(label, new TDoubleArrayList());
						}
						map.get(label).add(amount.getValue());
						break;
					}
				}
			}
		}

		for (final QuantificationLabel label : labels) {
			final TDoubleList intensities = map.get(label);
			if (intensities == null) {
				continue;
			}
			if (!"".equals(sb.toString())) {
				sb.append("\t");
			}
			if (aggregation == null) {
				sb.append(intensities.get(0));
			} else {
				double value = 0.0;
				switch (aggregation) {
				case AVG:
					value = Maths.mean(intensities);
					break;
				case HIGHEST:
					value = intensities.max();
					break;
				case SUM:
					value = intensities.sum();
				default:
					break;
				}
				sb.append(value);
			}

		}
		return sb.toString();
	}

	public static String getIntensitiesString(QuantifiedPSMInterface psm, List<QuantificationLabel> labels) {
		final StringBuilder sb = new StringBuilder();
		final Set<Amount> amounts = psm.getAmounts();

		for (final QuantificationLabel label : labels) {

			for (final Amount amount : amounts) {
				final QuantAmount quantAmount = (QuantAmount) amount;
				if (amount.getAmountType() == AmountType.INTENSITY) {
					if (quantAmount.getLabel() == label) {
						if (!"".equals(sb.toString())) {
							sb.append("\t");
						}
						sb.append(amount.getValue());
						break;
					}
				}
			}
		}

		return sb.toString();
	}

	public static List<String> getProteinSiteKeys(List<PTMInProtein> ptmsInProtein, boolean mapToGene,
			List<QuantifiedProteinInterface> primaryProteins) {
		Set<String> primaryAccs = null;
		if (primaryProteins != null) {
			primaryAccs = primaryProteins.stream().map(p -> p.getAccession()).collect(Collectors.toSet());
		}
		final List<String> retTMP = new ArrayList<String>();
		for (final PTMInProtein ptmInProtein : ptmsInProtein) {
			if (primaryAccs != null && !primaryAccs.contains(ptmInProtein.getProteinACC())) {
				continue;
			}
			boolean valid = false;
			for (final float deltaMass : PTMList.getInstance()) {
				if (Math.abs(ptmInProtein.getDeltaMass() - deltaMass) < MASS_PRECISION) {
					valid = true;
					break;
				}
			}
			if (valid) {
				String gene = null;
				if (mapToGene) {
					try {
						gene = GeneMapper.getInstance().getGene(ptmInProtein.getProteinACC());
					} catch (final IOException e) {

					}
				}
				if (gene == null) {
					gene = ptmInProtein.getProteinACC();
				}
				retTMP.add(gene + "_" + ptmInProtein.getAa() + ptmInProtein.getPosition());
			}
		}

		// now we check whether we can merge some elements from the same protein/gene
		final List<String> ret = new ArrayList<String>();
		String currentProtein = null;
		for (final String string : retTMP) {
			final String protein = string.substring(0, string.indexOf("_"));
			if (protein.equals(currentProtein)) {
				final String newString = ret.get(ret.size() - 1) + "_" + string.substring(string.indexOf("_") + 1);
				ret.set(ret.size() - 1, newString);
			} else {
				currentProtein = protein;
				ret.add(string);
			}

		}
		return ret;
	}

	public static String getGenes(List<QuantifiedProteinInterface> proteins) {
		final StringBuilder sb = new StringBuilder();
		final List<String> genes = new ArrayList<String>();
		for (final QuantifiedProteinInterface protein : proteins) {
			String gene = FastaParser.getGeneFromFastaHeader(protein.getDescription());
			if (FastaParser.isReverse(protein.getAccession())) {
				gene = null;
			}

			if (gene == null) {
				if (!genes.contains("N/A")) {
					genes.add("N/A");
				}
			} else {
				if (!genes.contains(gene)) {
					genes.add(gene);
				}
			}
		}
		for (final String gene : genes) {
			if ("N/A".equals(gene) && genes.size() > 1) {
				continue;
			}
			if (!"".equals(sb.toString())) {
				sb.append(SEPARATOR);
			}

			sb.append(gene);
		}
		return sb.toString();
	}

	/**
	 * Discards TrEmbl proteins from the groups when at least one Swiss-Prot protein
	 * is found. In other words, just keep SwissProt proteins only, unless there is
	 * none, and in that case, we keep all.
	 * 
	 * @param groups
	 * @return
	 */
	public static List<QuantifiedProteinInterface> getPrimaryProteins(Set<ProteinGroup> groups) {
		final List<QuantifiedProteinInterface> ret = new ArrayList<QuantifiedProteinInterface>();
		final Set<QuantifiedProteinInterface> proteins = new THashSet<QuantifiedProteinInterface>();
		final Set<String> accs = new THashSet<String>();
		for (final ProteinGroup proteinGroup : groups) {
			for (final GroupableProtein groupableProtein : proteinGroup) {
				proteins.add((QuantifiedProteinInterface) groupableProtein);
			}
			accs.addAll(proteinGroup.getAccessions());
		}

		final Map<String, Entry> annotatedProteins = UPLR.getInstance().getAnnotatedProteins(null, accs);

		for (final QuantifiedProteinInterface protein : proteins) {
			final Entry entry = annotatedProteins.get(protein.getAccession());
			if (entry != null && UniprotEntryUtil.isSwissProt(entry)) {
				ret.add(protein);
			}
		}

		if (ret.isEmpty()) {
			// we add all
			ret.addAll(proteins);
		}

		return ret;
	}

	public static String getAccessions(List<QuantifiedProteinInterface> proteins) {
		final StringBuilder sb = new StringBuilder();
		for (final QuantifiedProteinInterface protein : proteins) {
			if (!"".equals(sb.toString())) {

				sb.append(SEPARATOR);
			}
			sb.append(protein.getAccession());
		}
		return sb.toString();
	}

	public static String getGenesString(Set<ProteinGroup> groups) {
		final StringBuilder sb = new StringBuilder();
		for (final ProteinGroup proteinGroup : groups) {
			if (!"".equals(sb.toString())) {
				sb.append(SEPARATOR);
			}
			final StringBuilder sb2 = new StringBuilder();
			for (final GroupableProtein protein : proteinGroup) {
				if (!"".equals(sb2.toString())) {
					sb2.append(",");
				}
				final String description = ((QuantifiedProteinInterface) protein).getDescription();
				String gene = FastaParser.getGeneFromFastaHeader(description);
				if (FastaParser.isReverse(protein.getAccession())) {
					gene = null;
				}
				if (gene != null) {
					sb2.append(gene);
				} else {
					sb2.append("N/A");
				}
			}
			sb.append(sb2.toString());
		}
		return sb.toString();
	}

	public static String getDescriptionString(Set<ProteinGroup> groups) {
		final StringBuilder sb = new StringBuilder();
		for (final ProteinGroup proteinGroup : groups) {
			if (!"".equals(sb.toString())) {
				sb.append(SEPARATOR);
			}
			final StringBuilder sb2 = new StringBuilder();
			for (final GroupableProtein protein : proteinGroup) {
				if (!"".equals(sb2.toString())) {
					sb2.append(",");
				}
				sb2.append(((QuantifiedProteinInterface) protein).getDescription());
			}
			sb.append(sb2.toString());
		}
		return sb.toString();
	}

	public static String getProteinGroupAccessionString(Set<ProteinGroup> groups) {
		final StringBuilder sb = new StringBuilder();
		for (final ProteinGroup proteinGroup : groups) {
			if (!"".equals(sb.toString())) {
				sb.append(SEPARATOR);
			}
			final StringBuilder sb2 = new StringBuilder();
			for (final String acc : proteinGroup.getAccessions()) {
				if (!"".equals(sb2.toString())) {
					sb2.append(",");
				}
				sb2.append(acc);
			}
			sb.append(sb2.toString());
		}
		return sb.toString();
	}

	public static int getIonCount(QuantifiedPSMInterface psm, QuantParser parser) {
		final int reCalculatedIonCount = parser.getReCalculatedIonCount(psm);
		return reCalculatedIonCount;
	}

	public static int getIonCountAfterFilters(QuantifiedPSMInterface psm, QuantParser parser) {
		final String ionKey = KeyUtils.getInstance().getSequenceChargeKey(psm, true, true);
		final Set<QuantifiedPSMInterface> psms = parser.getPSMsByIonKey().get(ionKey);
		if (psms != null) {
			return Long.valueOf(psms.stream().filter(psm2 -> !psm2.isDiscarded()).count()).intValue();
		}
		return 0;
	}

	public static Set<ProteinGroup> getGroups(HasProteins o) {
		final Set<ProteinGroup> collect = o.getProteins().stream().map(p -> p.getProteinGroup())
				.filter(g -> g.getEvidence() != ProteinEvidence.NONCONCLUSIVE).collect(Collectors.toSet());
		return collect;
	}

	public static String getSignalToNoise(QuantifiedPSMInterface psm) {
		final Score score = getScoreObject(psm, CensusOutParser.SIGNAL_TO_NOISE);
		if (score != null) {
			return score.getValue();
		}
		return "";
	}

	public static String getTMTPurity(QuantifiedPSMInterface psm) {
		final Score score = getScoreObject(psm, CensusOutParser.TMT_PURITY);
		if (score != null) {
			return score.getValue();
		}
		return "";
	}

	public static float getScore(QuantifiedPSMInterface psm, String scoreName) {
		final Score score = getScoreObject(psm, scoreName);

		if (score != null) {
			try {
				return Float.valueOf(score.getValue());
			} catch (final NumberFormatException e) {

			}
		}

		return Float.NaN;
	}

	public static Score getScoreObject(HasScores scoreHolder, String scoreName) {
		final Optional<Score> findAny = scoreHolder.getScores().stream().filter(s -> s.getScoreName().equals(scoreName))
				.findAny();
		if (findAny.isPresent()) {
			return findAny.get();
		}
		return null;
	}

	public static boolean isTMTPurityValid(QuantifiedPSMInterface psm, Float minPurity) {
		if (minPurity == null) {
			return true;
		}
		final float purity = getScore(psm, CensusOutParser.TMT_PURITY);
		if (Double.isNaN(purity) || Float.compare(-1.0f, purity) == 0) {
			return true;
		}
		if (purity >= minPurity) {
			return true;
		}
		return false;

	}

	public static List<String> getAccToPrint(List<QuantifiedPSMInterface> psms) {
		final List<String> ret = new ArrayList<String>();
		for (final QuantifiedPSMInterface psm : psms) {
			final Set<QuantifiedProteinInterface> proteins = psm.getQuantifiedProteins();
			for (final QuantifiedProteinInterface protein : proteins) {
				if (!ret.contains(protein.getAccession())) {
					ret.add(protein.getAccession());
				}
			}
		}
		return ret;
	}

	public static int getChargeToPrint(Collection<QuantifiedPSMInterface> psms) {
		int charge = -1;
		for (final QuantifiedPSMInterface psm : psms) {
			if (charge == -1) {
				charge = psm.getChargeState();
			}
			if (psm.getChargeState() != charge) {
				throw new IllegalArgumentException(
						"Incosistency. Groupping by sequence and charge is not properly done!");
			}
		}
		return charge;
	}

	public static Map<QuantificationLabel, Double> getAggregatedIntensities(Collection<QuantifiedPSMInterface> psms,
			PSMAggregationType psmSelection, boolean useRawIntensity) throws IOException {

		switch (psmSelection) {
		case HIGHEST:
			return getHighestIntensitiesToPrint(psms, useRawIntensity);
		case AVG:
			return getAveragedIntensitiesToPrint(psms, useRawIntensity);
		case SUM:
			return getSummedIntensitiesToPrint(psms, useRawIntensity);
		default:
			throw new IllegalArgumentException(psmSelection + " PSM selection is not supported yet");
		}

	}

	private static Map<QuantificationLabel, Double> getHighestIntensitiesToPrint(
			Collection<QuantifiedPSMInterface> psms, boolean useRawIntensity) throws IOException {
		QuantifiedPSMInterface highestPSM = null;
		double highestSum = -Double.MAX_VALUE;
		for (final QuantifiedPSMInterface psm : psms) {
			final Map<QuantificationLabel, Double> ret = getIntensitiesFromPSM(psm, useRawIntensity);
			double sum = 0.0;
			for (final double intensity : ret.values()) {
				sum += intensity;
			}
			if (sum > highestSum) {
				highestPSM = psm;
				highestSum = sum;
			}
		}
		// highest PSm is the one to report
		return getIntensitiesFromPSM(highestPSM, useRawIntensity);
	}

	private static Map<QuantificationLabel, Double> getIntensitiesFromPSM(QuantifiedPSMInterface psm,
			boolean useRawIntensity) throws IOException {
		final Map<QuantificationLabel, Double> ret = new EnumMap<QuantificationLabel, Double>(
				QuantificationLabel.class);
		final Set<Amount> amounts = psm.getAmounts();
		for (final Amount amount : amounts) {
			final QuantAmount quantAmount = (QuantAmount) amount;
			if (quantAmount.getAmountType() == AmountType.NORMALIZED_INTENSITY && useRawIntensity) {
				continue;
			} else if (quantAmount.getAmountType() == AmountType.INTENSITY && !useRawIntensity) {
				continue;
			} else if (quantAmount.getAmountType() == AmountType.INTENSITY
					|| quantAmount.getAmountType() == AmountType.NORMALIZED_INTENSITY) {
				final QuantificationLabel label = quantAmount.getLabel();
				ret.put(label, quantAmount.getValue());
			}
		}
		return ret;
	}

	private static Map<QuantificationLabel, Double> getAveragedIntensitiesToPrint(
			Collection<QuantifiedPSMInterface> psms, boolean useRawIntensity) throws IOException {
		final Map<QuantificationLabel, TDoubleList> toAverage = new EnumMap<QuantificationLabel, TDoubleList>(
				QuantificationLabel.class);
		for (final QuantifiedPSMInterface psm : psms) {
			final Set<Amount> amounts = psm.getAmounts();
			for (final Amount amount : amounts) {
				final QuantAmount quantAmount = (QuantAmount) amount;
				if (quantAmount.getAmountType() == AmountType.NORMALIZED_INTENSITY && useRawIntensity) {
					continue;
				} else if (quantAmount.getAmountType() == AmountType.INTENSITY && !useRawIntensity) {
					continue;
				} else if (quantAmount.getAmountType() == AmountType.INTENSITY
						|| quantAmount.getAmountType() == AmountType.NORMALIZED_INTENSITY) {
					final QuantificationLabel label = quantAmount.getLabel();
					// for the average, don't use the zero
					if (Double.compare(0.0, quantAmount.getValue()) == 0) {
						continue;
					}
					if (!toAverage.containsKey(label)) {
						toAverage.put(label, new TDoubleArrayList());
					}

					toAverage.get(label).add(amount.getValue());
				}
			}
		}
		final Map<QuantificationLabel, Double> ret = new EnumMap<QuantificationLabel, Double>(
				QuantificationLabel.class);
		for (final QuantificationLabel label : toAverage.keySet()) {
			ret.put(label, Maths.mean(toAverage.get(label)));
		}
		return ret;
	}

	private static Map<QuantificationLabel, Double> getSummedIntensitiesToPrint(Collection<QuantifiedPSMInterface> psms,
			boolean useRawIntensity) throws IOException {
		final Map<QuantificationLabel, Double> ret = new EnumMap<QuantificationLabel, Double>(
				QuantificationLabel.class);
		for (final QuantifiedPSMInterface psm : psms) {
			final Set<Amount> amounts = psm.getAmounts();
			for (final Amount amount : amounts) {
				final QuantAmount quantAmount = (QuantAmount) amount;
				if (quantAmount.getAmountType() == AmountType.NORMALIZED_INTENSITY && useRawIntensity) {
					continue;
				} else if (quantAmount.getAmountType() == AmountType.INTENSITY && !useRawIntensity) {
					continue;
				} else if (quantAmount.getAmountType() == AmountType.INTENSITY
						|| quantAmount.getAmountType() == AmountType.NORMALIZED_INTENSITY) {

					final QuantificationLabel label = quantAmount.getLabel();
					if (!ret.containsKey(label)) {
						ret.put(label, amount.getValue());
					} else {
						ret.put(label, ret.get(label) + amount.getValue());
					}
				}
			}
		}
		return ret;
	}

	public static EnumMap<QuantificationLabel, Set<QuantAmount>> getQuantAmountsByLabels(QuantifiedPSMInterface psm) {
		final EnumMap<QuantificationLabel, Set<QuantAmount>> ret = new EnumMap<QuantificationLabel, Set<QuantAmount>>(
				QuantificationLabel.class);
		final Set<Amount> amounts = psm.getAmounts();
		for (final Amount amount : amounts) {
			if (amount instanceof QuantAmount) {
				final QuantAmount qamount = (QuantAmount) amount;
				if (!ret.containsKey(qamount.getLabel())) {
					ret.put(qamount.getLabel(), new THashSet<QuantAmount>());
				}
				ret.get(qamount.getLabel()).add(qamount);
			}
		}
		return ret;
	}

	public static Set<ProteinGroup> getGroups(List<QuantifiedPSMInterface> psms) {
		final Set<ProteinGroup> ret = new THashSet<ProteinGroup>();
		for (final QuantifiedPSMInterface psm : psms) {
			psm.getProteins().forEach(protein -> ret.add(protein.getProteinGroup()));
		}
		return ret;
	}
}
