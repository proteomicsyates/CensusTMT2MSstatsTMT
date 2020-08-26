package edu.scripps.yates.censustmt2msstatstmt.filters;

import java.io.File;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.log4j.Logger;

import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.censustmt2msstatstmt.singletons.FilesManager;
import edu.scripps.yates.censustmt2msstatstmt.singletons.ProteinSequences;
import edu.scripps.yates.censustmt2msstatstmt.singletons.UPLR;
import edu.scripps.yates.censustmt2msstatstmt.util.DataUtil;
import edu.scripps.yates.utilities.grouping.GroupablePeptide;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class IonFilter {
	private final static Logger log = Logger.getLogger(IonFilter.class);

	public static List<QuantifiedPSMInterface> applyMinIonsPerProteinFilter(Collection<ProteinGroup> groups,
			int minIons) {
		final Set<QuantifiedPSMInterface> validPSMs = new THashSet<QuantifiedPSMInterface>();
		final Set<String> discardedProteins = new THashSet<String>();
		for (final ProteinGroup group : groups) {
			final Set<String> ionKeys = new THashSet<String>();
			for (final GroupablePeptide groupablePeptide : group.getPSMs()) {
				final QuantifiedPSMInterface psm = (QuantifiedPSMInterface) groupablePeptide;

				final String ionKey = KeyUtils.getInstance().getSequenceChargeKey(psm, true, true);
				ionKeys.add(ionKey);
			}
			if (ionKeys.size() < minIons) {
				group.getPSMs().forEach(psm -> ((QuantifiedPSMInterface) psm).setDiscarded(true));
				discardedProteins.add(group.getAccessionString(DataUtil.SEPARATOR));
			} else {
				group.getPSMs().forEach(psm -> validPSMs.add((QuantifiedPSMInterface) psm));
			}
		}

		System.out.println(discardedProteins.size() + " proteins were discarded by minimum number of ions " + minIons);
		final File outFile = FilesManager.getInstance().printDiscardedListByIons(discardedProteins, "Protein", minIons);
		System.out.println("Discarded proteins written to file '" + outFile.getAbsolutePath() + "'");

		System.out.println(validPSMs.size()
				+ " PSMs pass the filter of belonging to an ion with minimum number of spec counts " + minIons);
		final List<QuantifiedPSMInterface> collect = validPSMs.stream().distinct().collect(Collectors.toList());
		return collect;
	}

	public static List<QuantifiedPSMInterface> applyMinIonsPerProteinSite(Collection<QuantifiedPSMInterface> psms,
			int minIons) {
		final Map<String, Set<QuantifiedPSMInterface>> psmsByProteinSiteKey = new THashMap<String, Set<QuantifiedPSMInterface>>();
		for (final QuantifiedPSMInterface psm : psms) {
			if (psm.isDiscarded()) {
				continue;
			}

			try {
				final List<PTMInProtein> ptmsInProtein = psm.getPTMsInProtein(UPLR.getInstance(),
						ProteinSequences.getInstance());
				final List<String> proteinSiteKeys = DataUtil.getProteinSiteKeys(ptmsInProtein, false, null);
				for (final String proteinSiteKey : proteinSiteKeys) {
					if (!psmsByProteinSiteKey.containsKey(proteinSiteKey)) {
						psmsByProteinSiteKey.put(proteinSiteKey, new THashSet<QuantifiedPSMInterface>());
					}
					psmsByProteinSiteKey.get(proteinSiteKey).add(psm);
				}
			} catch (final IllegalArgumentException e) {

			}
		}
		final Set<String> discardedProteinSites = new THashSet<String>();
		for (final String proteinSiteKey : psmsByProteinSiteKey.keySet()) {
			if (proteinSiteKey.equals("Q9Y2W1_T874")) {
				System.out.println("asdf");
			}
			final Set<String> ionKeys = new THashSet<String>();
			final Set<QuantifiedPSMInterface> psmsInSite = psmsByProteinSiteKey.get(proteinSiteKey);
			for (final QuantifiedPSMInterface psm : psmsInSite) {
				final String ionKey = KeyUtils.getInstance().getSequenceChargeKey(psm, true, true);
				ionKeys.add(ionKey);
			}
			if (ionKeys.size() < minIons) {
				discardedProteinSites.add(proteinSiteKey);
			}
		}
		final Set<QuantifiedPSMInterface> validPSMs = new THashSet<QuantifiedPSMInterface>();
		for (final QuantifiedPSMInterface psm : psms) {
			if (psm.isDiscarded()) {
				continue;
			}

			final List<PTMInProtein> ptmsInProtein = psm.getPTMsInProtein(UPLR.getInstance(),
					ProteinSequences.getInstance());
			final List<String> proteinSiteKeys = DataUtil.getProteinSiteKeys(ptmsInProtein, false, null);
			// if all protein site keys are discarded, then discard the psm too
			boolean anyValid = false;
			for (final String proteinSiteKey : proteinSiteKeys) {
				if (!discardedProteinSites.contains(proteinSiteKey)) {
					anyValid = true;
					break;
				}
			}
			if (anyValid) {
				validPSMs.add(psm);
			} else {
				psm.setDiscarded(true);
			}
		}

		System.out.println(
				discardedProteinSites.size() + " protein sites were discarded by minimum number of ions " + minIons);
		final File outFile = FilesManager.getInstance().printDiscardedListByIons(discardedProteinSites, "proteinSites",
				minIons);
		System.out.println("Discarded protein sites written to file '" + outFile.getAbsolutePath() + "'");

		final List<QuantifiedPSMInterface> collect = validPSMs.stream().distinct().collect(Collectors.toList());
		System.out.println(collect.size()
				+ " PSMs pass the filter of map a PTM of interest to a protein site with minimum number of ion counts "
				+ minIons);
		return collect;
	}
}
