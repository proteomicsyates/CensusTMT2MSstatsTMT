package edu.scripps.yates.censustmt2msstatstmt.singletons;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.Collections;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.apache.tools.ant.DirectoryScanner;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.model.interfaces.QuantParser;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPeptideInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.censustmt2msstatstmt.ExperimentalDesign;
import edu.scripps.yates.censustmt2msstatstmt.Mixture;
import edu.scripps.yates.censustmt2msstatstmt.PSMAggregationType;
import edu.scripps.yates.censustmt2msstatstmt.SPC_FILTER_TYPE;
import edu.scripps.yates.censustmt2msstatstmt.util.DataUtil;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.files.FileUtils;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.sequence.PTMInProtein;
import edu.scripps.yates.utilities.strings.StringUtils;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class FilesManager implements Clearable {
	private final static Logger log = Logger.getLogger(FilesManager.class);
	private static FilesManager instance;
	private final Map<File, String> sheetsByFiles = new THashMap<File, String>();
	private final Map<String, File> filesBySheetNames = new THashMap<String, File>();
	private final List<File> outputfiles = new ArrayList<File>();
	private final List<String> sheetNames = new ArrayList<String>();
	private final List<File> inputFiles;
	private final File inputFolder;
	private final String outputPrefix;

	public static FilesManager getInstance() {
		return instance;
	}

	public static FilesManager getNewInstance(List<File> inputFiles, String outputPrefix) {
		instance = new FilesManager(inputFiles, outputPrefix);
		return instance;
	}

	/**
	 * Creates an storage for the output files based on the inputFiles
	 * 
	 * @param inputFiles
	 */
	private FilesManager(List<File> inputFiles, String outputPrefix) {
		this.inputFiles = inputFiles;
		this.inputFolder = inputFiles.get(0).getParentFile();
		this.outputPrefix = outputPrefix;
	}

	public void addOutputFile(File outputFile, String sheetName) {
		this.sheetsByFiles.put(outputFile, sheetName);
		this.filesBySheetNames.put(sheetName, outputFile);
		this.outputfiles.add(outputFile);
		this.sheetNames.add(sheetName);
	}

	public String getSheetName(File outputFile) {
		return this.sheetsByFiles.get(outputFile);
	}

	private File getExcelOutputFile() {
		return getOutputFile("", "xlsx");
	}

	public void generateExcelFile() throws IOException {
		final long t1 = System.currentTimeMillis();
		final File excelOutputFile = getExcelOutputFile();
		System.out.println(
				"Creating Excel file with all text files as different sheets. This can take a few minutes depending on the size of the data. Please wait...");

		if (excelOutputFile.exists()) {
			System.out.println(
					"Excel file at '" + excelOutputFile.getAbsolutePath() + "' already exists. It will be overriden.");
			final boolean deleted = excelOutputFile.delete();
			if (!deleted) {
				throw new IllegalArgumentException(
						"Excel file at '" + excelOutputFile.getAbsolutePath() + "' couldn't be deleted");
			}
		}
		Collections.sort(sheetNames);
		for (final String sheetName : sheetNames) {
			System.out.println("Adding sheet " + sheetName + "...");
			final File file = filesBySheetNames.get(sheetName);
			try {
				FileUtils.separatedValuesToXLSX(file.getAbsolutePath(), excelOutputFile.getAbsolutePath(), "\t",
						sheetName);
			} catch (final Exception e) {
				e.printStackTrace();
			}
		}
		final long t2 = System.currentTimeMillis() - t1;

		System.out.println("File created at: '" + excelOutputFile.getAbsolutePath() + "' in "
				+ DatesUtil.getDescriptiveTimeFromMillisecs(t2));

	}

	public File getPSMOutputFile() {
		return getOutputFile("_PSMs", "txt");
	}

	public File getPeptideOutputFile() {
		return getOutputFile("_Peptides", "txt");
	}

	/**
	 * 
	 * @param suffix    will be appended directly to the input file name
	 * @param extension without the '.'
	 * @return
	 */
	private File getOutputFile(String suffix, String extension) {

		return new File(inputFolder.getAbsolutePath() + File.separator + outputPrefix + suffix + "." + extension);
	}

	public File getMSStatsOutputFile() {
		return new File(this.inputFolder.getAbsolutePath() + File.separator + outputPrefix + "_msstatsTMT.txt");

	}

	public static File getRunningFolder() {

		final File file = new File(System.getProperty("user.dir"));
		return file;
	}

	public static List<File> getMultipleFiles(File parentFolder, String iOptionValue) {
		final DirectoryScanner scanner = new DirectoryScanner();
		final String[] includes = new String[] { FilenameUtils.getName(iOptionValue) };
		scanner.setIncludes(includes);
		scanner.setBasedir(parentFolder);
		scanner.setCaseSensitive(false);
		scanner.scan();
		final String[] files = scanner.getIncludedFiles();
		final List<File> ret = new ArrayList<File>();
		for (final String fileFullPath : files) {
//			if (FilenameUtils.getName(fileFullPath).startsWith(SUFFIX)) {
//				System.out.println("Ignoring file matching -i parameter: '" + fileFullPath + "'");
//				continue;
//			}
			ret.add(new File(parentFolder.getAbsolutePath() + File.separator + fileFullPath));

		}
		return ret;
	}

	public File printList(Set<String> list, String suffix) {
		final File outputFile = getOutputFile(suffix, "txt");
		FileWriter fw = null;
		try {
			fw = new FileWriter(outputFile);
			for (final String element : list) {
				fw.write(element + "\n");
			}
		} catch (final IOException e) {
			e.printStackTrace();
		} finally {
			if (fw != null) {
				try {
					fw.close();
				} catch (final IOException e) {
					e.printStackTrace();
				}
				String sheetName = suffix;
				if (sheetName.startsWith("_")) {
					sheetName = sheetName.substring(1);
				}
				addOutputFile(outputFile, sheetName);
				return outputFile;
			}
		}

		return null;
	}

	public List<File> getInputFiles() {
		return inputFiles;
	}

	public File printDiscardedList(Set<String> discardedIons, SPC_FILTER_TYPE spcFilterType, int minSPC) {
		final File outputFile = printList(discardedIons, "_discarded_" + minSPC + "SpCper" + spcFilterType.name());
		return outputFile;
	}

	public void printPSMLevelFile(Collection<QuantifiedPSMInterface> psms, QuantParser parser,
			List<QuantificationLabel> labels) throws IOException {
		System.out.println("Printing PSM-level file");
		final File outputFile = getPSMOutputFile();
		final FileWriter fw = new FileWriter(outputFile);
		// header
		fw.write("acc(s)\tgene(s)\tdescription(s)\tprimary acc\tprimary gene\tsequence\tz\tSpC before filters\t");
		if (LuciphorIntegrator.getInstance() != null) {
			fw.write("Modified by Luciphor\t");
		}
		// intensities in header
		for (final QuantificationLabel label : labels) {
			fw.write(label.name() + "\t");
		}
		// localization scores
//		final int maxNumberOfPTMs = maxNumberOfPTMs(psms);
//		for (int i = 1; i <= maxNumberOfPTMs; i++) {
//			fw.write("Localization score " + i + "\t");
//		}
		if (LuciphorIntegrator.getInstance() != null) {
			fw.write("Luciphor localFDR\tLuciphor globalFDR\t");
		}
		fw.write("TMT_purity\tsignal to noise\tion count\tscan number\tfile name\n");
		for (final QuantifiedPSMInterface psm : psms) {
			if (psm.isDiscarded()) {
				continue;
			}

			final StringBuilder sb = new StringBuilder();
			final Set<ProteinGroup> groups = DataUtil.getGroups(psm);
			// accs and genes
			if (!PTMList.getInstance().isEmpty()) {
				String sitesWithAccs = "";
				String sitesWithGenes = "";
				try {
					final List<PTMInProtein> ptMsInProtein = psm.getPTMsInProtein(UPLR.getInstance(),
							ProteinSequences.getInstance());
					final List<String> proteinSiteKeys = DataUtil.getProteinSiteKeys(ptMsInProtein, false, null);
					sitesWithAccs = StringUtils.getSeparatedValueStringFromChars(proteinSiteKeys.toArray(), ",");
					final List<String> proteinSiteKeysWithGenes = DataUtil.getProteinSiteKeys(ptMsInProtein, true,
							null);
					sitesWithGenes = StringUtils.getSeparatedValueStringFromChars(proteinSiteKeysWithGenes.toArray(),
							",");
				} catch (final IllegalArgumentException e) {
					sitesWithAccs = DataUtil.getProteinGroupAccessionString(groups);
					sitesWithGenes = DataUtil.getGenesString(groups);
				}
				sb.append(sitesWithAccs).append("\t");

				sb.append(sitesWithGenes).append("\t");
				if (sitesWithGenes.equals("ATP1B1_Y204,ATP1B1_T214")) {
					log.info("asdf");
				}
			} else {
				sb.append(DataUtil.getProteinGroupAccessionString(groups)).append("\t");
				// genes
				sb.append(DataUtil.getGenesString(groups)).append("\t");
			}

			// descriptions
			sb.append(DataUtil.getDescriptionString(groups)).append("\t");
			// primary acc and primary gene
			final List<QuantifiedProteinInterface> primaryProteins = DataUtil.getPrimaryProtein(groups);

			if (!PTMList.getInstance().isEmpty()) {
				String sitesWithAccs = "";
				String sitesWithGenes = "";
				try {
					if (psm.toString().equals("20181230_P-IGF1_2-67489-Y(+79.966)NPNVLPVQCT(+79.966)GK-2")) {
						log.info("adf");
					}
					final List<PTMInProtein> ptMsInProtein = psm.getPTMsInProtein(UPLR.getInstance(),
							ProteinSequences.getInstance());
					final List<String> proteinSiteKeys = DataUtil.getProteinSiteKeys(ptMsInProtein, false,
							primaryProteins);
					sitesWithAccs = StringUtils.getSeparatedValueStringFromChars(proteinSiteKeys.toArray(), ",");
					final List<String> proteinSiteKeysWithGenes = DataUtil.getProteinSiteKeys(ptMsInProtein, true,
							primaryProteins);
					sitesWithGenes = StringUtils.getSeparatedValueStringFromChars(proteinSiteKeysWithGenes.toArray(),
							",");
				} catch (final IllegalArgumentException e) {
					sitesWithAccs = DataUtil.getAccessions(primaryProteins);
					sitesWithGenes = DataUtil.getGenes(primaryProteins);
				}
				sb.append(sitesWithAccs).append("\t");

				sb.append(sitesWithGenes).append("\t");
			} else {
				sb.append(DataUtil.getAccessions(primaryProteins)).append("\t");

				sb.append(DataUtil.getGenes(primaryProteins)).append("\t");
			}

			// sequence
			sb.append(psm.getFullSequence()).append("\t");
			// charge state
			sb.append(psm.getChargeState()).append("\t");
			// Spc
			sb.append(psm.getPeptide().getPSMs().size()).append("\t");

			if (LuciphorIntegrator.getInstance() != null) {
				// modified by Luciphor
				if (LuciphorIntegrator.getInstance().getModifiedPSMsbyLuciphor() != null) {
					sb.append(LuciphorIntegrator.getInstance().getModifiedPSMsbyLuciphor().contains(psm.getKey()));
				}
				sb.append("\t");
			}
			// intensities
			sb.append(DataUtil.getIntensitiesString(psm, labels)).append("\t");
			// localization scores
//			sb.append(getLocalizationScoresString(psm)).append("\t");
			if (LuciphorIntegrator.getInstance() != null) {
				sb.append(LuciphorIntegrator.getInstance().getLuciphorFDRStrings(psm)).append("\t");
			}
			// tmt purity
			sb.append(DataUtil.getTMTPurity(psm)).append("\t");
			// signal to noise
			sb.append(DataUtil.getSignalToNoise(psm)).append("\t");
			// ion count
			sb.append(DataUtil.getIonCount(psm, parser)).append("\t");
			// scan
			sb.append(psm.getScanNumber()).append("\t");

			// file name
			sb.append(psm.getMSRun().getRunId()).append("\t");
			fw.write(sb.toString() + "\n");
		}
		fw.close();
		addOutputFile(outputFile, "psm-level");
		System.out.println("File created at " + outputFile.getAbsolutePath());
	}

	@Override
	public void clearData() {
		this.outputfiles.clear();
		this.sheetNames.clear();
		this.inputFiles.clear();
		this.outputfiles.clear();
		this.sheetsByFiles.clear();
	}

	public void printPeptideLevelFile(Collection<QuantifiedPeptideInterface> peptides, PSMAggregationType aggregation,
			List<QuantificationLabel> labels) throws IOException {
		System.out.println("Printing peptide-level file...");
		final File outputFile = getPeptideOutputFile();
		final FileWriter fw = new FileWriter(outputFile);
		// header
		fw.write(
				"acc(s)\tgene(s)\tdescription(s)\tprimary acc\tprimary gene\tsequence\tz\tSpC before filters\tSpC after filters");
		// intensities in header
		for (final QuantificationLabel label : labels) {
			fw.write("\t" + label.name());
		}
		fw.write("\n");
		for (final QuantifiedPeptideInterface peptide : peptides) {
			final boolean hasValidPSMs = peptide.getQuantifiedPSMs().stream().filter(psm -> !psm.isDiscarded())
					.findAny().isPresent();
			if (!hasValidPSMs) {
				continue;
			}
			// filter by psms
			final int numPSMs = (int) peptide.getQuantifiedPSMs().stream().filter(p -> !p.isDiscarded()).count();

			final StringBuilder sb = new StringBuilder();
			final Set<ProteinGroup> groups = DataUtil.getGroups(peptide);
			// acc
			if (!PTMList.getInstance().isEmpty()) {
				String sitesWithAccs = "";
				String sitesWithGenes = "";
				try {
					final List<PTMInProtein> ptMsInProtein = peptide.getPTMsInProtein(UPLR.getInstance(),
							ProteinSequences.getInstance());
					final List<String> proteinSiteKeys = DataUtil.getProteinSiteKeys(ptMsInProtein, false, null);
					sitesWithAccs = StringUtils.getSeparatedValueStringFromChars(proteinSiteKeys.toArray(), ",");
					final List<String> proteinSiteKeysWithGenes = DataUtil.getProteinSiteKeys(ptMsInProtein, true,
							null);
					sitesWithGenes = StringUtils.getSeparatedValueStringFromChars(proteinSiteKeysWithGenes.toArray(),
							",");
				} catch (final IllegalArgumentException e) {
					sitesWithAccs = DataUtil.getProteinGroupAccessionString(groups);
					sitesWithGenes = DataUtil.getGenesString(groups);
				}
				sb.append(sitesWithAccs).append("\t");

				sb.append(sitesWithGenes).append("\t");
			} else {
				sb.append(DataUtil.getProteinGroupAccessionString(groups)).append("\t");
				// genes
				sb.append(DataUtil.getGenesString(groups)).append("\t");
			}
			// descriptions
			sb.append(DataUtil.getDescriptionString(groups)).append("\t");
			// primary acc and gene
			final List<QuantifiedProteinInterface> primaryProteins = DataUtil.getPrimaryProtein(groups);
			if (!PTMList.getInstance().isEmpty()) {
				String sitesWithAccs = "";
				String sitesWithGenes = "";
				try {
					final List<PTMInProtein> ptMsInProtein = peptide.getPTMsInProtein(UPLR.getInstance(),
							ProteinSequences.getInstance());
					final List<String> proteinSiteKeys = DataUtil.getProteinSiteKeys(ptMsInProtein, false,
							primaryProteins);
					sitesWithAccs = StringUtils.getSeparatedValueStringFromChars(proteinSiteKeys.toArray(), ",");
					final List<String> proteinSiteKeysWithGenes = DataUtil.getProteinSiteKeys(ptMsInProtein, true,
							primaryProteins);
					sitesWithGenes = StringUtils.getSeparatedValueStringFromChars(proteinSiteKeysWithGenes.toArray(),
							",");
				} catch (final IllegalArgumentException e) {
					sitesWithAccs = DataUtil.getProteinGroupAccessionString(groups);
					sitesWithGenes = DataUtil.getGenesString(groups);
				}
				sb.append(sitesWithAccs).append("\t");
				// primary gene
				sb.append(sitesWithGenes).append("\t");

			} else {
				sb.append(DataUtil.getAccessions(primaryProteins)).append("\t");
				// primary gene
				sb.append(DataUtil.getGenes(primaryProteins)).append("\t");
			}
			// sequence
			sb.append(peptide.getFullSequence()).append("\t");
			// charge state
			sb.append(peptide.getPSMs().get(0).getChargeState()).append("\t");
			// Spc before filter
			sb.append(peptide.getPSMs().size()).append("\t");
			// Spc after filter
			sb.append(numPSMs).append("\t");

			// intensities
			sb.append(DataUtil.getIntensitiesString(peptide, aggregation, labels)).append("\n");
			fw.write(sb.toString());
		}
		fw.close();

		System.out.println("File created at " + outputFile.getAbsolutePath());
		FilesManager.getInstance().addOutputFile(outputFile, "peptide-level");

	}

	public String getInputFileNamesString() {
		final StringBuilder sb = new StringBuilder();
		for (final File file : inputFiles) {
			if (!"".equals(sb.toString())) {
				sb.append(", ");
			}
			sb.append(FilenameUtils.getName(file.getAbsolutePath()));

		}
		return sb.toString();
	}

	public void printMSstatsTMTFFile(List<QuantifiedPSMInterface> psms, ExperimentalDesign experimentalDesign,
			boolean useRawIntensity, PSMAggregationType psmSelection) throws IOException {
		final File outputFile = FilesManager.getInstance().getMSStatsOutputFile();
		FileWriter fw = null;
		try {
			fw = new FileWriter(outputFile);
			System.out.println("Output file: " + outputFile.getAbsolutePath());
			// header
			fw.write(
					"ProteinName\tPeptideSequence\tCharge\tPSM\tMixture\tTechRepMixture\tRun\tChannel\tCondition\tBioReplicate\tIntensity\n");

			// create a map by sequence
			final Map<String, List<QuantifiedPSMInterface>> psmsBySequenceAndCharge = new THashMap<String, List<QuantifiedPSMInterface>>();
			for (final QuantifiedPSMInterface psm : psms) {
				final String seq = psm.getSequence() + "_" + psm.getChargeState();

				if (!psmsBySequenceAndCharge.containsKey(seq)) {
					psmsBySequenceAndCharge.put(seq, new ArrayList<QuantifiedPSMInterface>());
				}
				psmsBySequenceAndCharge.get(seq).add(psm);
			}

			// iterate over it. Sort first the sequences so that they will appear in order
			// in the output
			final Set<String> seqsSet = psmsBySequenceAndCharge.keySet();
			final List<String> seqList = seqsSet.stream().sorted().collect(Collectors.toCollection(ArrayList::new));
			final Set<String> validProteins = new THashSet<String>();
			for (final String seq : seqList) {
				final List<QuantifiedPSMInterface> psmsOfSeq = psmsBySequenceAndCharge.get(seq);

				final List<String> accs = DataUtil.getAccToPrint(psmsOfSeq);
				final String acc = StringUtils.getSortedSeparatedValueStringFromChars(accs, ",");

				validProteins.add(acc);

				final int charge = DataUtil.getChargeToPrint(psmsOfSeq);
				final Set<String> runs = psmsOfSeq.stream().map(psm -> psm.getMSRun().getRunId())
						.collect(Collectors.toSet());
				for (final String run : runs) {
					final List<QuantifiedPSMInterface> psmsOfSeqAndRun = psmsOfSeq.stream()
							.filter(psm -> psm.getMSRun().getRunId().equals(run)).collect(Collectors.toList());
					final Map<QuantificationLabel, Double> intensities = DataUtil.getIntensitiesToPrint(psmsOfSeqAndRun,
							psmSelection, useRawIntensity);
					final QuantifiedPSMInterface firstPSM = psmsOfSeqAndRun.get(0);
					final String peptideSequence = firstPSM.getFullSequence();
					final String psm = peptideSequence + "_" + charge;

					for (final QuantificationLabel label : intensities.keySet().stream().sorted()
							.collect(Collectors.toList())) {
						printMSstatsTMTOutputLine(fw, label, intensities.get(label), acc, charge, peptideSequence, psm,
								run, experimentalDesign);
					}
				}
			}
			System.out.println("\n*******");
			System.out.println("Valid data: " + psms.size() + " PSMs");
			System.out.println(psmsBySequenceAndCharge.size() + " peptides (sequences+charge)");
			System.out.println(validProteins.size() + " proteins\n*******");

		} finally {
			if (fw != null) {
				fw.close();
				log.info("File written at " + outputFile.getAbsolutePath());
				System.out.println("File ready at " + outputFile.getAbsolutePath());
				addOutputFile(outputFile, "MSstatsTMT");
			}
		}
	}

	private void printMSstatsTMTOutputLine(FileWriter fw, QuantificationLabel label, double intensity, String acc,
			int charge, String peptideSequence, String psm, String run, ExperimentalDesign experimentalDesign)
			throws IOException {
		final Mixture mixture = experimentalDesign.getMixtureByRun(run);
		final Float channel = mixture.getChannelByLabel(label);
		final QuantCondition condition = mixture.getConditionsByLabels().get(label);

		final String techRepMixture = experimentalDesign.getTechRepMixtureByRun(run);
		final String bioReplicate = experimentalDesign.getBioReplicate(channel, techRepMixture, mixture.getName());
		final String channelString = Float.valueOf(channel).toString();

		fw.write(acc + "\t" + peptideSequence + "\t" + charge + "\t" + psm + "\t" + mixture.getName() + "\t"
				+ techRepMixture + "\t" + run + "\t" + channelString + "\t" + condition.getName() + "\t" + bioReplicate
				+ "\t" + intensity + "\n");

	}

	public static void clearStaticData() {
		if (instance != null) {
			instance.clearData();
		}

	}
}
