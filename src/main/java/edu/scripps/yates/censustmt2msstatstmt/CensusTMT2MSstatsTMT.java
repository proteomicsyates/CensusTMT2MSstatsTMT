package edu.scripps.yates.censustmt2msstatstmt;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.EnumMap;
import java.util.Iterator;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.Option;
import org.apache.commons.cli.ParseException;
import org.apache.commons.io.FilenameUtils;
import org.apache.log4j.Logger;
import org.apache.tools.ant.DirectoryScanner;
import org.springframework.core.io.ClassPathResource;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.CensusOutParser;
import edu.scripps.yates.census.read.QuantParserException;
import edu.scripps.yates.census.read.WrongTMTLabels;
import edu.scripps.yates.census.read.model.QuantAmount;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.utilities.appversion.AppVersion;
import edu.scripps.yates.utilities.maths.Maths;
import edu.scripps.yates.utilities.progresscounter.ProgressCounter;
import edu.scripps.yates.utilities.progresscounter.ProgressPrintingType;
import edu.scripps.yates.utilities.properties.PropertiesUtil;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.swing.CommandLineProgramGuiEnclosable;
import edu.scripps.yates.utilities.swing.DoNotInvokeRunMethod;
import edu.scripps.yates.utilities.swing.SomeErrorInParametersOcurred;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class CensusTMT2MSstatsTMT extends CommandLineProgramGuiEnclosable {
	private final static Logger log = Logger.getLogger(CensusTMT2MSstatsTMT.class);
	private static final String SUFFIX = "_msstatsTMT.txt";
	private static AppVersion version;
	private List<File> inputFiles = new ArrayList<File>();
	private boolean uniquePeptides;
	private boolean useRawIntensity;
	private PSMSelectionType psmSelection;
	private File experimentalDesignFile;
	private String experimentalDesignSeparator;
	private String decoyPrefix;
	private ExperimentalDesign ed;
	private int minNumPeptides;

//	public CensusTMT2MSstatsTMT(File inputFile, File experimentalDesignFile, String experimentalDesignSeparator,
//			boolean uniquePeptides, boolean useRawIntensity, String decoyPrefix, PSMSelectionType psmSelection) {
//		super(null);// we dont use that functionality to create GUI here
//		this.inputFile = inputFile;
//		this.experimentalDesignFile = experimentalDesignFile;
//		this.experimentalDesignSeparator = experimentalDesignSeparator;
//		this.uniquePeptides = uniquePeptides;
//		this.useRawIntensity = useRawIntensity;
//		this.decoyPrefix = decoyPrefix;
//		this.psmSelection = psmSelection;
//	}

	public CensusTMT2MSstatsTMT(String[] args) throws SomeErrorInParametersOcurred, ParseException {
		super(args);

	}

	public static void main(String[] args) {

		CensusTMT2MSstatsTMT c = null;
		try {
			c = new CensusTMT2MSstatsTMT(args);
			c.safeRun();
		} catch (final DoNotInvokeRunMethod e) {
			// do nothing
		} catch (final Exception e) {
			e.printStackTrace();
		}
	}

	private String getInputFileNamesString() {
		final StringBuilder sb = new StringBuilder();
		for (final File file : inputFiles) {
			if (!"".equals(sb.toString())) {
				sb.append(", ");
			}
			sb.append(FilenameUtils.getName(file.getAbsolutePath()));

		}
		return sb.toString();
	}

	private void runConversion() throws Exception {
		FileWriter fw = null;
		try {
			ed = null;
			System.out.println("Starting conversion of file " + getInputFileNamesString());
			final File outputFile = getOutputFile();
			fw = new FileWriter(outputFile);
			System.out.println("Output file: " + outputFile.getAbsolutePath());
			// header
			printHeader(fw);
			//
			final Map<QuantificationLabel, QuantCondition>[] conditionsByLabelsList = getConditionsByLabelsArray(
					this.inputFiles, getExperimentalDesign());
			String plural = "";
			if (this.inputFiles.size() > 1) {
				plural = "s";
			}
			System.out.println("Reading input file" + plural + "...");
			final CensusOutParser parser = new CensusOutParser(this.inputFiles.toArray(new File[0]),
					conditionsByLabelsList, null, null);
			if (parser.isTMT10().values().iterator().next()) {
				System.out.println("TMT 10-Plex detected.");
			} else if (parser.isTMT6().values().iterator().next()) {
				System.out.println("TMT 6-Plex detected.");
			} else if (parser.isTMT11().values().iterator().next()) {
				System.out.println("TMT 11-Plex detected.");
			}
			final List<QuantifiedPSMInterface> psms = new ArrayList<QuantifiedPSMInterface>();
			psms.addAll(parser.getPSMMap().values());
			System.out.println(psms.size() + " PSMs read from input file.");
			// filter by min number of peptides per protein

			final ProgressCounter counter = new ProgressCounter(psms.size(), ProgressPrintingType.PERCENTAGE_STEPS, 0);
			Iterator<QuantifiedPSMInterface> psmsIterator = psms.iterator();
			if (decoyPrefix != null && !"".equals(decoyPrefix)) {
				final Set<QuantifiedProteinInterface> discardedProteins = new THashSet<QuantifiedProteinInterface>();
				final int initialPSMs = psms.size();
				while (psmsIterator.hasNext()) {
					final QuantifiedPSMInterface psm = psmsIterator.next();

					for (final QuantifiedProteinInterface protein : psm.getQuantifiedProteins()) {
						// //decoy
						if (decoyPrefix != null && protein.getAccession().startsWith(decoyPrefix)) {
							psmsIterator.remove();
							discardedProteins.add(protein);
							break;
						}
					}
				}
				System.out.println(discardedProteins.size() + " proteins discarded as DECOYs.");
				System.out.println(
						"Now working with " + psms.size() + " PSMs (" + (initialPSMs - psms.size()) + " discarded).");
			}
			if (minNumPeptides > 1 || uniquePeptides) {
				int psmsDiscardedByUniqueness = 0;
				final Set<QuantifiedProteinInterface> discardedProteins = new THashSet<QuantifiedProteinInterface>();
				final int initialPSMs = psms.size();
				// filtering by minNumPeptides
				psmsIterator = psms.iterator();
				while (psmsIterator.hasNext()) {
					counter.increment();
					final String printIfNecessary = counter.printIfNecessary();
					if (!"".equals(printIfNecessary)) {
						System.out.println("Filtering data..." + printIfNecessary);
					}
					final QuantifiedPSMInterface psm = psmsIterator.next();
					if (uniquePeptides && psm.getQuantifiedProteins().size() > 1) {
						psmsIterator.remove();
						psmsDiscardedByUniqueness++;
						continue;
					}
					boolean anyProteinIsValid = false;
					for (final QuantifiedProteinInterface protein : psm.getQuantifiedProteins()) {
						final Set<String> peptideCharges = new THashSet<String>();
						protein.getQuantifiedPSMs().stream()
								.forEach(p -> peptideCharges.add(p.getFullSequence() + p.getChargeState()));
						if (peptideCharges.size() >= minNumPeptides) {
							anyProteinIsValid = true;
						} else {
							discardedProteins.add(protein);
						}
					}
					if (!anyProteinIsValid) {
						psmsIterator.remove();
					}
				}
				System.out.println(discardedProteins.size() + " proteins discarded for not having at least "
						+ minNumPeptides + " peptides (sequence+charge).");
				if (uniquePeptides) {
					System.out.println(psmsDiscardedByUniqueness + " PSMs discarded because they are not unique.");
				}
				System.out.println(
						" Now working with " + psms.size() + " PSMs (" + (initialPSMs - psms.size()) + " discarded).");

			}

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

				final List<String> accs = getAccToPrint(psmsOfSeq);
				final String acc = StringUtils.getSortedSeparatedValueStringFromChars(accs, ",");

				validProteins.add(acc);

				final int charge = getChargeToPrint(psmsOfSeq);
				final Set<String> runs = psmsOfSeq.stream().map(psm -> psm.getMSRun().getRunId())
						.collect(Collectors.toSet());
				for (final String run : runs) {
					final List<QuantifiedPSMInterface> psmsOfSeqAndRun = psmsOfSeq.stream()
							.filter(psm -> psm.getMSRun().getRunId().equals(run)).collect(Collectors.toList());
					final Map<QuantificationLabel, Double> intensities = getIntensitiesToPrint(psmsOfSeqAndRun,
							this.psmSelection, this.useRawIntensity);
					final QuantifiedPSMInterface firstPSM = psmsOfSeqAndRun.get(0);
					final String peptideSequence = firstPSM.getFullSequence();
					final String psm = peptideSequence + "_" + charge;

					for (final QuantificationLabel label : intensities.keySet().stream().sorted()
							.collect(Collectors.toList())) {
						printOutputLine(fw, label, intensities.get(label), acc, charge, peptideSequence, psm, run);
					}
				}
			}
			System.out.println("\n*******");
			System.out.println("Valid data: " + psms.size() + " PSMs");
			System.out.println(psmsBySequenceAndCharge.size() + " peptides (sequences+charge)");
			System.out.println(validProteins.size() + " proteins\n*******");

			log.info("File written at " + outputFile.getAbsolutePath());
			System.out.println("File ready at " + outputFile.getAbsolutePath());
		} catch (final Exception e) {
			if (e instanceof WrongTMTLabels) {
				throw new WrongTMTLabels(
						"Error in annotation file. The TMT channels in annotation file don't match with TMT channels in input file: "
								+ e.getMessage());
			} else {
				throw e;
			}
		} finally {
			if (fw != null) {
				fw.close();
			}
		}
	}

	/**
	 * It returns a Map<QuantificationLabel, QuantCondition> per one if the
	 * inputFiles. <br>
	 * To do that, it reads each of the files and looks for the run names, so it
	 * figures which mixture belongs to, and returns the mixture's Map
	 * 
	 * @param inputFiles2
	 * @param experimentalDesign
	 * @return
	 */
	private Map<QuantificationLabel, QuantCondition>[] getConditionsByLabelsArray(List<File> inputFiles2,
			ExperimentalDesign experimentalDesign) {
		final Map<QuantificationLabel, QuantCondition>[] ret = new Map[inputFiles2.size()];
		int index = 0;
		if (inputFiles2.size() == 1) {
			ret[0] = experimentalDesign.getMixtures().iterator().next().getConditionsByLabels();
			return ret;
		}
		for (final File inputFile : inputFiles2) {

			try {
				final CensusOutParser parser = new CensusOutParser(inputFile, null);
				parser.setIgnoreACCFormat(true);
				parser.setIgnoreTaxonomies(true);
				final Set<String> runs = parser.getPSMMap().values().stream().map(psm -> psm.getMSRun().getRunId())
						.collect(Collectors.toSet());
				for (final Mixture mixture : experimentalDesign.getMixtures()) {
					final Set<String> runsInMixture = mixture.getRuns();
					boolean valid = true;
					for (final String runInMixture : runsInMixture) {
						if (!runs.contains(runInMixture)) {
							valid = false;
							break;
						}
					}
					if (valid) {
						ret[index] = mixture.getConditionsByLabels();
					}
				}
				if (ret[index] == null) {
					throw new IllegalArgumentException("Runs in input file '" + inputFile.getAbsolutePath() + "' ("
							+ StringUtils.getSeparatedValueStringFromChars(runs.toArray(), ",")
							+ ") are not found in the annotation file.");
				}
			} catch (final FileNotFoundException e) {
			} catch (final QuantParserException e) {
			} finally {
				index++;
			}
		}

		return ret;

	}

	private void printHeader(FileWriter fw) throws IOException {
		fw.write(
				"ProteinName\tPeptideSequence\tCharge\tPSM\tMixture\tTechRepMixture\tRun\tChannel\tCondition\tBioReplicate\tIntensity\n");

	}

	private Map<QuantificationLabel, Double> getIntensitiesToPrint(List<QuantifiedPSMInterface> psms,
			PSMSelectionType psmSelection, boolean useRawIntensity) throws IOException {

		switch (psmSelection) {
		case HIGHEST:
			return getHighestIntensitiesToPrint(psms, useRawIntensity);
		case AVERAGE:
			return getAveragedIntensitiesToPrint(psms, useRawIntensity);
		case SUM:
			return getSummedIntensitiesToPrint(psms, useRawIntensity);
		default:
			throw new IllegalArgumentException(psmSelection + " PSM selection is not supported yet");
		}

	}

	private Map<QuantificationLabel, Double> getHighestIntensitiesToPrint(List<QuantifiedPSMInterface> psms,
			boolean useRawIntensity) throws IOException {
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

	private Map<QuantificationLabel, Double> getIntensitiesFromPSM(QuantifiedPSMInterface psm, boolean useRawIntensity)
			throws IOException {
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

	private Map<QuantificationLabel, Double> getAveragedIntensitiesToPrint(List<QuantifiedPSMInterface> psms,
			boolean useRawIntensity) throws IOException {
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

	private Map<QuantificationLabel, Double> getSummedIntensitiesToPrint(List<QuantifiedPSMInterface> psms,
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

	private int getChargeToPrint(List<QuantifiedPSMInterface> psms) {
		int charge = -1;
		for (final QuantifiedPSMInterface psm : psms) {
			if (charge == -1) {
				charge = psm.getChargeState();
			}
			if (psm.getChargeState() != charge) {
				throw new IllegalArgumentException(
						"Incosistency. GRoupping by sequence and charge is not properly done!");
			}
		}
		return charge;
	}

	private List<String> getAccToPrint(List<QuantifiedPSMInterface> psms) {
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

	private void printOutputLine(FileWriter fw, QuantificationLabel label, double intensity, String acc, int charge,
			String peptideSequence, String psm, String run) throws IOException {
		final ExperimentalDesign experimentalDesign = getExperimentalDesign();
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

	private File getOutputFile() {
		return new File(this.inputFiles.get(0).getParent() + File.separator
				+ FilenameUtils.getBaseName(this.inputFiles.get(0).getAbsolutePath()) + SUFFIX);
	}

	private ExperimentalDesign getExperimentalDesign() throws IOException {
		if (ed == null) {
			ed = new ExperimentalDesign(experimentalDesignFile, experimentalDesignSeparator);

		}
		return ed;
	}

	public static AppVersion getVersion() {
		if (version == null) {
			try {
				final String tmp = PropertiesUtil
						.getProperties(new ClassPathResource(AppVersion.APP_PROPERTIES).getInputStream())
						.getProperty("assembly.dir");
				if (tmp.contains("v")) {
					version = new AppVersion(tmp.split("v")[1]);
				} else {
					version = new AppVersion(tmp);
				}
			} catch (final Exception e) {
				e.printStackTrace();
			}
		}
		return version;

	}

	@Override
	protected List<Option> defineCommandLineOptions() {
		// create Options object
		final List<Option> options = new ArrayList<Option>();

		// add t option
		final Option option1 = new Option("i", "input", true,
				"Path to the input file(s). It can refer to multiple files by using wildcard '*', i.e: '/path/to/my/files/census*.out'");
		option1.setRequired(true);
		options.add(option1);
		final Option option2 = new Option("an", "annotation", true, "Path to the experimental design file.");
		option2.setRequired(true);
		options.add(option2);

		final Option option3 = new Option("ps", "psm_selection", true,
				"[OPTIONAL] What to do with multiple PSMs of the same peptide (SUM, AVERAGE or HIGHEST). If not provided, HIGHEST will be choosen.");
		options.add(option3);

		final Option option4 = new Option("u", "unique", false,
				"[OPTIONAL] Use only unique peptides. If not provided, all peptides will be used.");
		options.add(option4);

		final Option option5 = new Option("r", "raw", false,
				"[OPTIONAL] Use of raw intensity. If not provided, normalized intensity will be used.");
		options.add(option5);

		final Option option6 = new Option("d", "decoy", true,
				"[OPTIONAL] Remove decoy hits. Decoys hits will have this prefix in their accession number. If not provided, no decoy filtering will be used.");
		options.add(option6);

		final Option option7 = new Option("m", "minPeptides", true,
				"[OPTIONAL] Minimum number of peptides+charge per protein. If not provided, even proteins with 1 peptide will be quantified");
		options.add(option7);
		return options;
	}

	@Override
	public void run() {
		try {
			runConversion();
		} catch (final Exception e) {
			e.printStackTrace();
			throw new RuntimeException(e);
		}
	}

	@Override
	public String getTitleForFrame() {
		return "CensusTMT2MSstatsTMT converter";
	}

	@Override
	protected void initToolFromCommandLineOptions(CommandLine cmd) throws SomeErrorInParametersOcurred {
		File inputFile = null;
		File experimentalDesignFile = null;
		final String experimentalDesignSeparator = ",";
		boolean uniquePeptides = false;
		boolean useRawIntensity = false;
		String decoyPrefix = null;
		PSMSelectionType psmSelection = PSMSelectionType.HIGHEST;
		if (cmd.hasOption("i")) {
			this.inputFiles = new ArrayList<File>();
			final String iOptionValue = cmd.getOptionValue("i");
			inputFile = new File(iOptionValue);
			final File parentFile = inputFile.getParentFile();
			if (parentFile == null && !inputFile.exists()) {
				inputFile = new File(System.getProperty("user.dir") + File.separator + iOptionValue);
			}

			if (!inputFile.exists()) {
				// check whether is using a wildcard in the name
				if (iOptionValue.matches(".*\\*.*")) {
					final List<File> inputFiles = getMultipleFiles(parentFile, iOptionValue);
					if (inputFiles.isEmpty()) {
						errorInParameters("Input files '-i " + iOptionValue + "' don't exist or are not found");
					} else {
						this.inputFiles.addAll(inputFiles);
					}
				} else {
					errorInParameters("Input file '-i " + iOptionValue + "' doesn't exist or is not found");
				}
			} else {
				inputFiles.add(inputFile);
			}
		} else {
			errorInParameters("Input file is missing");
		}
		String plural = "";
		if (this.inputFiles.size() > 1) {
			plural += "s";
		}
		System.out.println("Input file" + plural + ":");
		int i = 1;
		for (final File file : inputFiles) {
			System.out.println(i++ + "- " + file.getAbsolutePath() + "'");
		}
		if (cmd.hasOption("an")) {
			final String anOptionValue = cmd.getOptionValue("an");
			experimentalDesignFile = new File(anOptionValue);
			final File parentFile = experimentalDesignFile.getParentFile();
			if (parentFile == null || !experimentalDesignFile.exists()) {
				experimentalDesignFile = new File(System.getProperty("user.dir") + File.separator + anOptionValue);
			}

			if (!experimentalDesignFile.exists()) {
				errorInParameters("Annotation file (experimental design) '-an " + anOptionValue
						+ "' doesn't exist or is not found");
			}
		} else {
			errorInParameters("Experimental design file is missing");
		}
		if (cmd.hasOption("ps")) {
			try {
				psmSelection = PSMSelectionType.valueOf(cmd.getOptionValue("ps"));
			} catch (final Exception e) {
				e.printStackTrace();
				errorInParameters("option 'ps' ('" + cmd.getOptionValue("ps") + "') is not valid. Valid values are "
						+ PSMSelectionType.printValidValues());
			}

		}
		log.info("Using ps option " + psmSelection);

		if (cmd.hasOption("u")) {
			uniquePeptides = true;
			log.info("Using only unique peptides.");
		} else {
			log.info("Using all peptides regardless of their uniqueness.");
		}

		if (cmd.hasOption("r")) {
			useRawIntensity = true;
			log.info("Using raw intensities.");
		} else {
			log.info("Using normalized intensities.");
		}

		if (cmd.hasOption("d")) {
			decoyPrefix = cmd.getOptionValue("d");
			if ("".equals(decoyPrefix)) {
				decoyPrefix = null;
			}
		}
		int minNumPeptides = 0;
		if (cmd.hasOption("m")) {
			minNumPeptides = Integer.valueOf(cmd.getOptionValue("m"));
		}
		log.info("Using d option for decoy prefix='" + decoyPrefix + "'");
		this.experimentalDesignFile = experimentalDesignFile;
		this.experimentalDesignSeparator = experimentalDesignSeparator;
		this.uniquePeptides = uniquePeptides;
		this.useRawIntensity = useRawIntensity;
		this.decoyPrefix = decoyPrefix;
		this.psmSelection = psmSelection;
		this.minNumPeptides = minNumPeptides;
	}

	private List<File> getMultipleFiles(File parentFolder, String iOptionValue) {
		final DirectoryScanner scanner = new DirectoryScanner();
		final String[] includes = new String[] { FilenameUtils.getName(iOptionValue) };
		scanner.setIncludes(includes);
		scanner.setBasedir(parentFolder);
		scanner.setCaseSensitive(false);
		scanner.scan();
		final String[] files = scanner.getIncludedFiles();
		final List<File> ret = new ArrayList<File>();
		for (final String fileFullPath : files) {
			if (FilenameUtils.getName(fileFullPath).endsWith(SUFFIX)) {
				System.out.println("Ignoring file matching -i parameter: '" + fileFullPath + "'");
				continue;
			}
			ret.add(new File(parentFolder.getAbsolutePath() + File.separator + fileFullPath));

		}
		return ret;
	}

	@Override
	public String printCommandLineSintax() {
		return "CensusTMT2MSstatsTMT -i [input file] -an [annotation file]";
	}

}
