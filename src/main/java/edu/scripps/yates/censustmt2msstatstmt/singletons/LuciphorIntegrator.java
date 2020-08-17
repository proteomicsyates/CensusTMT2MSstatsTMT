package edu.scripps.yates.censustmt2msstatstmt.singletons;

import java.io.File;
import java.io.IOException;
import java.util.Collection;
import java.util.Map;
import java.util.Set;

import edu.scripps.yates.census.read.QuantParserException;
import edu.scripps.yates.census.read.model.QuantifiedPSM;
import edu.scripps.yates.census.read.model.QuantifiedPeptide;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.censustmt2msstatstmt.util.DataUtil;
import edu.scripps.yates.utilities.luciphor.LuciphorReader;
import edu.scripps.yates.utilities.proteomicsmodel.PSM;
import edu.scripps.yates.utilities.proteomicsmodel.Score;
import edu.scripps.yates.utilities.proteomicsmodel.utils.KeyUtils;
import edu.scripps.yates.utilities.util.Pair;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class LuciphorIntegrator implements Clearable {
	private static LuciphorIntegrator instance;
	private Map<String, Pair<Float, Float>> luciphorFDRsByPSMScan = new THashMap<String, Pair<Float, Float>>();
	private File luciphorFile;
	private final Set<String> modifiedPSMsbyLuciphor = new THashSet<String>();

	public static LuciphorIntegrator getInstance() {
		return instance;
	}

	public static LuciphorIntegrator getInstance(File luciphorFile) {
		if (instance == null) {
			instance = new LuciphorIntegrator(luciphorFile);
		}
		return instance;
	}

	private LuciphorIntegrator(File luciphorFile) {
		this.luciphorFile = luciphorFile;
	}

	/**
	 * Merge census PSMS with luciphor PSMS, by modifying census ones setting new
	 * fullSequence (sequence with PTMs), and by assigning them a new peptide with
	 * the new sequence+PTM
	 * 
	 * @param censusPSMsByScan
	 * @param luciphorPSMsByScan
	 * @param minLuciphorFDR
	 * @return a set of PSM keys that were modified
	 */
	private Set<String> mergeLuciphorPSMs(Map<String, QuantifiedPSM> censusPSMsByScan,
			Map<String, PSM> luciphorPSMsByScan, Float localFDRThreshold, Float globalFDRThreshold,
			boolean chargeStateSensible) {
		final Set<String> ret = new THashSet<String>();
		// also we need the peptideMap
		final Map<String, QuantifiedPeptideInterface> peptideMap = new THashMap<String, QuantifiedPeptideInterface>();
		int newPeptides = 0;
		final Set<String> newPeptideKeys = new THashSet<String>();
		for (final String scanNumber : luciphorPSMsByScan.keySet()) {
			if (censusPSMsByScan.containsKey(scanNumber)) {
				final QuantifiedPSM psm = censusPSMsByScan.get(scanNumber);
				final PSM luciphorPSM = luciphorPSMsByScan.get(scanNumber);
				if (!passLuciphorGlobalFDRThreshold(luciphorPSM, globalFDRThreshold)) {
					continue;
				}
				if (!passLuciphorLocalFDRThreshold(luciphorPSM, localFDRThreshold)) {
					continue;
				}
				if (luciphorPSM.getFullSequence().equals(psm.getFullSequence())) {
					continue;
				}
				ret.add(luciphorPSM.getKey());
				// there has been a change
				final String luciphorPeptideKey = KeyUtils.getInstance().getSequenceChargeKey(luciphorPSM, true,
						chargeStateSensible);
				final QuantifiedPeptideInterface peptide = psm.getQuantifiedPeptide();
				// change full sequence to psm
				psm.setFullSequence(luciphorPSM.getFullSequence());
				// create peptide with new full sequence.
				// psm and newPeptide will be pointing each other after this
				QuantifiedPeptide newPeptide = null;
				if (peptideMap.containsKey(luciphorPeptideKey)) {
					newPeptide = (QuantifiedPeptide) peptideMap.get(luciphorPeptideKey);
					// we add psm to it
					newPeptide.addPSM(psm, true);
				} else {
					newPeptide = new QuantifiedPeptide(psm, true, true, true);
					final QuantifiedPeptide pep = newPeptide;
					// assign Amounts and Conditions to the new peptide from the old one
					peptide.getAmounts().stream().forEach(amount -> pep.addAmount(amount));
					peptide.getConditions().stream().forEach(condition -> pep.addCondition(condition));
				}
				newPeptideKeys.add(newPeptide.getKey());

				// remove previous peptide from peptideMap if doesn't have other psms that are
				// not modified by Luciphor
				boolean remove = true;
				for (final PSM psm2 : peptide.getPSMs()) {
					String scanNumber2 = psm2.getScanNumber();
					if (psm2.getScanNumber().contains("_")) {
						scanNumber2 = psm2.getScanNumber().split("_")[0];
					}
					if (!luciphorPSMsByScan.containsKey(scanNumber2)) {
						remove = false;
					}
				}
				if (remove) {
//					psmList.remove(peptide);
				}
				// add new peptide to peptide list
//				psmList.add(newPeptide);
				peptideMap.put(newPeptide.getKey(), newPeptide);
				newPeptides++;

			}
		}
		System.out.println(ret.size() + " PSMs where modified by Luciphor. ");
		if (newPeptides > 0) {
			System.out.println(newPeptides + " new peptides where created due to Luciphor changes in PTM positions.");
		}
		return ret;
	}

	/**
	 * 
	 * @param parser
	 * @param localFDRThreshold
	 * @param globalFDRThreshold
	 * @throws IOException
	 * @throws QuantParserException
	 */
	public void mergeWithLuciphor(Collection<QuantifiedPSMInterface> psms, Float localFDRThreshold,
			Float globalFDRThreshold) throws QuantParserException {
		if (luciphorFile != null && luciphorFile.exists()) {
			System.out.println("Merging PSMs with Luciphor results...");
			final LuciphorReader luciphorParser = new LuciphorReader(luciphorFile);
			final Map<String, PSM> luciphorPSMs = luciphorParser.getPSMs();
			// merge both luciphor and quant compare
			// we keep a map of PSMs from luciphor using scan as key
			final Map<String, PSM> luciphorPSMsByScan = new THashMap<String, PSM>();
			luciphorPSMs.values().stream().forEach(luciphorPSM -> luciphorPSMsByScan
					.put(luciphorPSM.getScanNumber() + "-" + luciphorPSM.getMSRun().getRunId(), luciphorPSM));
			luciphorFDRsByPSMScan = new THashMap<String, Pair<Float, Float>>();
			for (final PSM luciphorPSM : luciphorPSMs.values()) {
				float globalFDR = Float.NaN;
				float localFDR = Float.NaN;
				final Score globalFDRScore = DataUtil.getScoreObject(luciphorPSM, LuciphorReader.GLOBAL_FDR);
				if (globalFDRScore != null) {
					try {
						globalFDR = Float.valueOf(globalFDRScore.getValue());
					} catch (final NumberFormatException e) {

					}
				}
				final Score localFDRScore = DataUtil.getScoreObject(luciphorPSM, LuciphorReader.LOCAL_FDR);
				if (localFDRScore != null) {
					try {
						localFDR = Float.valueOf(localFDRScore.getValue());
					} catch (final NumberFormatException e) {

					}
				}
				final Pair<Float, Float> pair = new Pair<Float, Float>(localFDR, globalFDR);
				luciphorFDRsByPSMScan.put(luciphorPSM.getScanNumber() + "-" + luciphorPSM.getMSRun().getRunId(), pair);
			}

			// first we keep a map of PSMS from quant compare using scan as key
			final Map<String, QuantifiedPSM> psmsByScan = new THashMap<String, QuantifiedPSM>();
			for (final QuantifiedPSMInterface psm : psms) {
				psmsByScan.put(psm.getScanNumber() + "-" + psm.getMSRun().getRunId(), (QuantifiedPSM) psm);
			}

			modifiedPSMsbyLuciphor.addAll(
					mergeLuciphorPSMs(psmsByScan, luciphorPSMsByScan, localFDRThreshold, globalFDRThreshold, true));
			System.out.println(modifiedPSMsbyLuciphor.size() + " PSMs were modified by Luciphor results");
		}
	}

	private static boolean passLuciphorGlobalFDRThreshold(PSM luciphorPSM, Float minGlobalLuciphorFDR) {
		return passLuciphorFDRThreshold(luciphorPSM, minGlobalLuciphorFDR, LuciphorReader.GLOBAL_FDR);
	}

	private static boolean passLuciphorLocalFDRThreshold(PSM luciphorPSM, Float minLocalLuciphorFDR) {
		return passLuciphorFDRThreshold(luciphorPSM, minLocalLuciphorFDR, LuciphorReader.LOCAL_FDR);
	}

	private static boolean passLuciphorFDRThreshold(PSM luciphorPSM, Float minLuciphorFDR, String fdrTypeName) {
		if (minLuciphorFDR == null) {
			return true;
		}
		for (final Score score : luciphorPSM.getScores()) {
			if (score.getScoreName().equalsIgnoreCase(fdrTypeName)) {
				try {
					final double value = Double.valueOf(score.getValue());
					if (value < minLuciphorFDR) {
						return true;
					}
				} catch (final NumberFormatException e) {

				}
			}
		}
		return false;
	}

	public String getLuciphorFDRStrings(QuantifiedPSMInterface psm) {
		final String key = psm.getScanNumber() + "-" + psm.getMSRun().getRunId();
		if (this.luciphorFDRsByPSMScan != null && this.luciphorFDRsByPSMScan.containsKey(key)) {
			final Pair<Float, Float> pair = luciphorFDRsByPSMScan.get(key);
			return pair.getFirstelement() + "\t" + pair.getSecondElement();
		}
		return "\t";
	}

	public Set<String> getModifiedPSMsbyLuciphor() {
		return modifiedPSMsbyLuciphor;
	}

	public void setLuciphorFile(File luciphorFile) {
		this.luciphorFile = luciphorFile;
		clearData();
	}

	public static void clearDataStatic() {
		if (instance != null) {
			instance.clearData();
		}
	}

	@Override
	public void clearData() {
		luciphorFDRsByPSMScan.clear();
		modifiedPSMsbyLuciphor.clear();
	}

}
