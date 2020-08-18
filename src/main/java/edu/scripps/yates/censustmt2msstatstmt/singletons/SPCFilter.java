package edu.scripps.yates.censustmt2msstatstmt.singletons;

import java.io.File;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.censustmt2msstatstmt.SPC_FILTER_TYPE;
import edu.scripps.yates.censustmt2msstatstmt.util.DataUtil;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class SPCFilter {
	private final static Logger log = Logger.getLogger(SPCFilter.class);

	/**
	 * Apply the SPC filter depending on the type on the parameter spcFilterType
	 * 
	 * @param psms
	 * @param minSPC
	 * @param spcFilterType
	 * @return
	 */
	public static List<QuantifiedPSMInterface> applySPCFilter(Collection<QuantifiedPSMInterface> psms, Integer minSPC,
			SPC_FILTER_TYPE spcFilterType) {
		if (minSPC == null || minSPC == 0) {
			return psms.stream().distinct().collect(Collectors.toList());
		}

		switch (spcFilterType) {
		case ION:
			return applyMinSPCPerIonFilter(psms, minSPC, spcFilterType);
		case PROTEIN_SITE:
			return applyMinSPCPerProteinSite(psms, minSPC, spcFilterType);
		case SEQUENCE:
			return applyMinSPCPerSequence(psms, minSPC, spcFilterType);
		default:
			break;
		}
		throw new IllegalArgumentException();
	}

	private static List<QuantifiedPSMInterface> applyMinSPCPerSequence(Collection<QuantifiedPSMInterface> psms,
			Integer minSPC, SPC_FILTER_TYPE spcFilterType) {
		final Map<String, Set<QuantifiedPSMInterface>> psmsBySequenceKey = new THashMap<String, Set<QuantifiedPSMInterface>>();
		for (final QuantifiedPSMInterface psm : psms) {
			final String sequenceKey = KeyUtils.getInstance().getSequenceChargeKey(psm, true, false);
			if (!psmsBySequenceKey.containsKey(sequenceKey)) {
				psmsBySequenceKey.put(sequenceKey, new THashSet<QuantifiedPSMInterface>());
			}
			psmsBySequenceKey.get(sequenceKey).add(psm);
		}
		final Set<QuantifiedPSMInterface> ret = new THashSet<QuantifiedPSMInterface>();
		final Set<String> discardedSequences = new THashSet<String>();
		for (final String sequenceKey : psmsBySequenceKey.keySet()) {
			final Set<QuantifiedPSMInterface> psmSet = psmsBySequenceKey.get(sequenceKey);
			if (psmSet.size() >= minSPC) {
				ret.addAll(psmSet);
			} else {
				discardedSequences.add(sequenceKey);
			}
		}
		System.out.println(discardedSequences.size()
				+ " peptide sequences were discarded by minimum number of spec counts " + minSPC);
		final File outFile = FilesManager.getInstance().printDiscardedListBySPC(discardedSequences, spcFilterType, minSPC);
		System.out.println("Discarded peptides sequences written to file " + outFile.getAbsolutePath());

		final List<QuantifiedPSMInterface> collect = ret.stream().distinct().collect(Collectors.toList());
		System.out.println(collect.size()
				+ " PSMs pass the filter of belonging to a peptide sequence with minimum number of spec counts "
				+ minSPC);
		return collect;
	}

	private static List<QuantifiedPSMInterface> applyMinSPCPerIonFilter(Collection<QuantifiedPSMInterface> psms,
			int minSPC, SPC_FILTER_TYPE spcFilterType) {

		final Map<String, Set<QuantifiedPSMInterface>> psmsByIonKey = new THashMap<String, Set<QuantifiedPSMInterface>>();
		for (final QuantifiedPSMInterface psm : psms) {
			final String ionKey = KeyUtils.getInstance().getSequenceChargeKey(psm, true, true);
			if (!psmsByIonKey.containsKey(ionKey)) {
				psmsByIonKey.put(ionKey, new THashSet<QuantifiedPSMInterface>());
			}
			psmsByIonKey.get(ionKey).add(psm);
		}
		final Set<QuantifiedPSMInterface> ret = new THashSet<QuantifiedPSMInterface>();
		final Set<String> discardedIons = new THashSet<String>();
		for (final String ionKey : psmsByIonKey.keySet()) {
			final Set<QuantifiedPSMInterface> psmSet = psmsByIonKey.get(ionKey);
			if (psmSet.size() >= minSPC) {
				ret.addAll(psmSet);
			} else {
				discardedIons.add(ionKey);
			}
		}
		System.out.println(discardedIons.size()
				+ " ions (sequence+charge) were discarded by minimum number of spec counts " + minSPC);
		final File outFile = FilesManager.getInstance().printDiscardedListBySPC(discardedIons, spcFilterType, minSPC);
		System.out.println("Discarded ions written to file " + outFile.getAbsolutePath());

		final List<QuantifiedPSMInterface> collect = ret.stream().distinct().collect(Collectors.toList());
		System.out.println(collect.size()
				+ " PSMs pass the filter of belonging to an ion with minimum number of spec counts " + minSPC);
		return collect;
	}

	private static List<QuantifiedPSMInterface> applyMinSPCPerProteinSite(Collection<QuantifiedPSMInterface> psms,
			Integer minSPC, SPC_FILTER_TYPE spcFilterType) {
		final Map<String, Set<QuantifiedPSMInterface>> psmsByProteinSiteKey = new THashMap<String, Set<QuantifiedPSMInterface>>();
		final Map<QuantifiedPSMInterface, Set<String>> proteinSiteKeysByPSM = new THashMap<QuantifiedPSMInterface, Set<String>>();
		for (final QuantifiedPSMInterface psm : psms) {
			if (psm.isDiscarded()) {
				continue;
			}
			proteinSiteKeysByPSM.put(psm, new THashSet<String>());
			try {
				final List<PTMInProtein> ptmsInProtein = psm.getPTMsInProtein(UPLR.getInstance(),
						ProteinSequences.getInstance());
				final List<String> proteinSiteKeys = DataUtil.getProteinSiteKeys(ptmsInProtein, false, null);
				for (final String proteinSiteKey : proteinSiteKeys) {
					if (!psmsByProteinSiteKey.containsKey(proteinSiteKey)) {
						psmsByProteinSiteKey.put(proteinSiteKey, new THashSet<QuantifiedPSMInterface>());
					}
					psmsByProteinSiteKey.get(proteinSiteKey).add(psm);
					proteinSiteKeysByPSM.get(psm).add(proteinSiteKey);
				}
			} catch (final IllegalArgumentException e) {

			}
		}

		final Set<QuantifiedPSMInterface> ret = new THashSet<QuantifiedPSMInterface>();
		final Set<String> discardedProteinSiteKeys = new THashSet<String>();
		for (final QuantifiedPSMInterface psm : proteinSiteKeysByPSM.keySet()) {
			if ("20181230_P-IGF1_1-19035-AY(+79.966)SSPS(+79.966)TTPEAR-2".equals(psm.toString())) {
				log.info("asdf");
			}
			final Set<String> proteinSiteKeys = proteinSiteKeysByPSM.get(psm);
			boolean someValid = false;
			for (final String proteinSiteKey : proteinSiteKeys) {
				final Set<QuantifiedPSMInterface> psmSet = psmsByProteinSiteKey.get(proteinSiteKey);
				final int numValidPSMs = (int) psmSet.stream().filter(p -> !p.isDiscarded()).count();
				if (numValidPSMs >= minSPC) {
					someValid = true;
					break;
				}
			}
			if (someValid) {
				ret.add(psm);
			} else {
				discardedProteinSiteKeys.addAll(proteinSiteKeys);
			}
		}
		System.out.println(discardedProteinSiteKeys.size()
				+ " protein sites were discarded by minimum number of spec counts " + minSPC);
		final File outFile = FilesManager.getInstance().printDiscardedListBySPC(discardedProteinSiteKeys, spcFilterType,
				minSPC);
		System.out.println("Discarded protein sites written to file " + outFile.getAbsolutePath());

		final List<QuantifiedPSMInterface> collect = ret.stream().distinct().collect(Collectors.toList());
		System.out.println(collect.size()
				+ " PSMs pass the filter of map a PTM of interest to a protein site with minimum number of spec counts "
				+ minSPC);
		return collect;
	}
}
