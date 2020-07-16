package edu.scripps.yates.censustmt2msstatstmt;

import java.io.File;
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
import org.springframework.core.io.ClassPathResource;

import edu.scripps.yates.census.analysis.QuantCondition;
import edu.scripps.yates.census.read.CensusOutParser;
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
	private static AppVersion version;
	private File inputFile;
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

	public CensusTMT2MSstatsTMT(String[] args)
			throws ParseException, DoNotInvokeRunMethod, SomeErrorInParametersOcurred {
		super(args);

	}

	public static void main(String[] args) {

		CensusTMT2MSstatsTMT c = null;
		try {
			c = new CensusTMT2MSstatsTMT(args);
			c.run();
		} catch (final DoNotInvokeRunMethod e) {
			// do nothing
		} catch (final Exception e) {
			e.printStackTrace();
		}
	}

	private void runConversion() throws IOException {
		System.out.println("Starting conversion of file " + inputFile.getAbsolutePath());
		final File outputFile = getOutputFile();
		final FileWriter fw = new FileWriter(outputFile);
		System.out.println("Output file: " + outputFile.getAbsolutePath());
		// header
		printHeader(fw);
		final Map<QuantificationLabel, QuantCondition> conditionsByLabels = getExperimentalDesign()
				.getConditionByLabel();
		System.out.println("Reading input file...");
		final CensusOutParser parser = new CensusOutParser(this.inputFile, conditionsByLabels);
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
			System.out.print("\n-" + discardedProteins.size() + " proteins discarded as DECOYs.");
			System.out.println(
					" Now working with " + psms.size() + " PSMs (" + (initialPSMs - psms.size()) + " discarded).\n");
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
					if (protein.getAccession().equals("A0A096MJY1")) {
						log.info("asdf");
					}
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
			System.out.print("\n " + discardedProteins.size() + " proteins discarded for not having at least "
					+ minNumPeptides + " peptides (sequence+charge).");
			if (uniquePeptides) {
				System.out.print(" " + psmsDiscardedByUniqueness + " PSMs discarded because they are not unique.");
			}
			System.out.println(
					" Now working with " + psms.size() + " PSMs (" + (initialPSMs - psms.size()) + " discarded).");

		}

		// create a map by sequence
		final Map<String, List<QuantifiedPSMInterface>> psmsBySequenceAndCharge = new THashMap<String, List<QuantifiedPSMInterface>>();
		for (final QuantifiedPSMInterface psm : psms) {
			final String seq = psm.getBeforeSeq() + "." + psm.getSequence() + "." + psm.getAfterSeq() + "_"
					+ psm.getChargeState();
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
			final Map<QuantificationLabel, Double> intensities = getIntensitiesToPrint(psmsOfSeq, this.psmSelection,
					this.useRawIntensity);
			final QuantifiedPSMInterface firstPSM = psmsOfSeq.get(0);
			final String peptideSequence = firstPSM.getBeforeSeq() + "." + firstPSM.getSequence() + "."
					+ firstPSM.getAfterSeq();
			final String psm = peptideSequence + "_" + charge;
			final String run = firstPSM.getMSRun().getRunId();
			for (final QuantificationLabel label : intensities.keySet().stream().sorted()
					.collect(Collectors.toList())) {
				printOutputLine(fw, label, intensities.get(label), acc, charge, peptideSequence, psm, run);
			}

		}
		System.out.print("\n*******\nValid data: " + psms.size() + " PSMs");
		System.out.print(", " + psmsBySequenceAndCharge.size() + " peptides (sequences+charge)");
		System.out.println(", " + validProteins.size() + " proteins\n*******");
		fw.close();
		log.info("File written at " + outputFile.getAbsolutePath());
		System.out.println("File ready at " + outputFile.getAbsolutePath());
	}

	private void printHeader(FileWriter fw) throws IOException {
		fw.write(
				"ProteinName\tPeptideSequence\tCharge\tPSM\tMixture\tTechRepMixture\tRun\tChannel\tCondition\tBioReplicate\tIntensity\n");

	}

	private Map<QuantificationLabel, Double> getIntensitiesToPrint(List<QuantifiedPSMInterface> psms,
			PSMSelectionType psmSelection, boolean useRawIntensity) {

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
			boolean useRawIntensity) {
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

	private Map<QuantificationLabel, Double> getIntensitiesFromPSM(QuantifiedPSMInterface psm,
			boolean useRawIntensity) {
		final Map<QuantificationLabel, Double> ret = new EnumMap<QuantificationLabel, Double>(
				QuantificationLabel.class);
		final Set<Amount> amounts = psm.getAmounts();
		for (final Amount amount : amounts) {
			if (amount.getAmountType() == AmountType.NORMALIZED_INTENSITY && useRawIntensity) {
				continue;
			} else if (amount.getAmountType() == AmountType.INTENSITY && !useRawIntensity) {
				continue;
			} else if (amount.getAmountType() == AmountType.INTENSITY
					|| amount.getAmountType() == AmountType.NORMALIZED_INTENSITY) {
				final String name = amount.getCondition().getName();
				final String[] split = name.split(ExperimentalDesign.SYMBOL);
				final String conditionName = split[0];
				final double channel = Double.valueOf(split[1]);
				final QuantificationLabel label = ExperimentalDesign.getTMT6LabelFromChannel(channel);
				ret.put(label, amount.getValue());
			}
		}
		return ret;
	}

	private Map<QuantificationLabel, Double> getAveragedIntensitiesToPrint(List<QuantifiedPSMInterface> psms,
			boolean useRawIntensity) {
		final Map<QuantificationLabel, TDoubleList> toAverage = new EnumMap<QuantificationLabel, TDoubleList>(
				QuantificationLabel.class);
		for (final QuantifiedPSMInterface psm : psms) {
			final Set<Amount> amounts = psm.getAmounts();
			for (final Amount amount : amounts) {
				if (amount.getAmountType() == AmountType.NORMALIZED_INTENSITY && useRawIntensity) {
					continue;
				} else if (amount.getAmountType() == AmountType.INTENSITY && !useRawIntensity) {
					continue;
				} else if (amount.getAmountType() == AmountType.INTENSITY
						|| amount.getAmountType() == AmountType.NORMALIZED_INTENSITY) {
					final String name = amount.getCondition().getName();
					final String[] split = name.split(ExperimentalDesign.SYMBOL);
					final String conditionName = split[0];
					final double channel = Double.valueOf(split[1]);
					final QuantificationLabel label = ExperimentalDesign.getTMT6LabelFromChannel(channel);
					// for the average, don't use the zero
					if (Double.compare(0.0, amount.getValue()) == 0) {
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
			boolean useRawIntensity) {
		final Map<QuantificationLabel, Double> ret = new EnumMap<QuantificationLabel, Double>(
				QuantificationLabel.class);
		for (final QuantifiedPSMInterface psm : psms) {
			final Set<Amount> amounts = psm.getAmounts();
			for (final Amount amount : amounts) {
				if (amount.getAmountType() == AmountType.NORMALIZED_INTENSITY && useRawIntensity) {
					continue;
				} else if (amount.getAmountType() == AmountType.INTENSITY && !useRawIntensity) {
					continue;
				} else if (amount.getAmountType() == AmountType.INTENSITY
						|| amount.getAmountType() == AmountType.NORMALIZED_INTENSITY) {
					final String name = amount.getCondition().getName();
					final String[] split = name.split(ExperimentalDesign.SYMBOL);
//					final String conditionName = split[0];
					final double channel = Double.valueOf(split[1]);
					final QuantificationLabel label = ExperimentalDesign.getTMT6LabelFromChannel(channel);
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
		final double channel = Double.valueOf(ExperimentalDesign.getChannelByLabel(label));
		final QuantCondition condition = experimentalDesign.getConditionByLabel().get(label);
		final String conditionName = condition.getName().split(ExperimentalDesign.SYMBOL)[0];
		final String mixture = experimentalDesign.getMixtureByRun(run);
		final String techRepMixture = experimentalDesign.getTechRepMixtureByRun(run);
		final String bioReplicate = experimentalDesign.getBioReplicate(channel, techRepMixture, mixture);
		String channelString = Double.valueOf(channel).toString();
		if (experimentalDesign.isTMT6Plex()) {
			channelString = String.valueOf(Double.valueOf(channel).intValue());
		}
		fw.write(acc + "\t" + peptideSequence + "\t" + charge + "\t" + psm + "\t" + mixture + "\t" + techRepMixture
				+ "\t" + run + "\t" + channelString + "\t" + conditionName + "\t" + bioReplicate + "\t" + intensity
				+ "\n");

	}

	private File getOutputFile() {
		return new File(this.inputFile.getParent() + File.separator
				+ FilenameUtils.getBaseName(this.inputFile.getAbsolutePath()) + "_msstatsTMT.txt");
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
		final Option option1 = new Option("i", "input", true, "Path to the input file.");
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
	public void run() throws Exception {
		runConversion();
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
			final String iOptionValue = cmd.getOptionValue("i");
			inputFile = new File(iOptionValue);
			final File parentFile = inputFile.getParentFile();
			if (parentFile == null || !inputFile.exists()) {
				inputFile = new File(System.getProperty("user.dir") + File.separator + iOptionValue);
			}

			if (!inputFile.exists()) {
				errorInParameters("Input file '-i " + iOptionValue + "' doesn't exist or is not found");
			}
		} else {
			errorInParameters("Input file is missing");
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
		this.inputFile = inputFile;
		this.experimentalDesignFile = experimentalDesignFile;
		this.experimentalDesignSeparator = experimentalDesignSeparator;
		this.uniquePeptides = uniquePeptides;
		this.useRawIntensity = useRawIntensity;
		this.decoyPrefix = decoyPrefix;
		this.psmSelection = psmSelection;
		this.minNumPeptides = minNumPeptides;
	}

	@Override
	public String printCommandLineSintax() {
		return "CensusTMT2MSstatsTMT -i [input file] -an [annotation file]";
	}

}
