package edu.scripps.yates.censustmt2msstatstmt;

import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.List;
import java.util.Map;
import java.util.Set;
import java.util.stream.Collectors;

import org.apache.commons.cli.CommandLine;
import org.apache.commons.cli.CommandLineParser;
import org.apache.commons.cli.DefaultParser;
import org.apache.commons.cli.HelpFormatter;
import org.apache.commons.cli.Options;
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
import edu.scripps.yates.utilities.properties.PropertiesUtil;
import edu.scripps.yates.utilities.proteomicsmodel.Amount;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType;
import gnu.trove.list.TDoubleList;
import gnu.trove.list.array.TDoubleArrayList;
import gnu.trove.map.TObjectDoubleMap;
import gnu.trove.map.custom_hash.TObjectDoubleCustomHashMap;
import gnu.trove.map.hash.THashMap;
import gnu.trove.map.hash.TObjectDoubleHashMap;

public class CensusTMT2MSstatsTMT {
	private final static Logger log = Logger.getLogger(CensusTMT2MSstatsTMT.class);
	private static AppVersion version;
	private static Options options;
	private final File inputFile;
	private final boolean uniquePeptides;
	private final boolean useRawIntensity;
	private final PSMSelectionType psmSelection;
	private final File experimentalDesignFile;
	private final String experimentalDesignSeparator;
	private ExperimentalDesign ed;

	public CensusTMT2MSstatsTMT(File inputFile, File experimentalDesignFile, String experimentalDesignSeparator,
			boolean uniquePeptides, boolean useRawIntensity, PSMSelectionType psmSelection) {
		this.inputFile = inputFile;
		this.experimentalDesignFile = experimentalDesignFile;
		this.experimentalDesignSeparator = experimentalDesignSeparator;
		this.uniquePeptides = uniquePeptides;
		this.useRawIntensity = useRawIntensity;
		this.psmSelection = psmSelection;
	}

	public static void main(String[] args) {
		setupCommandLineOptions();
		File inputFile = null;
		File experimentalDesignFile = null;
		final String experimentalDesignSeparator = ",";
		boolean uniquePeptides = false;
		boolean useRawIntensity = false;
		PSMSelectionType psmSelection = PSMSelectionType.HIGHEST;
		final CommandLineParser parser = new DefaultParser();
		try {
			final CommandLine cmd = parser.parse(options, args);
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
			if (cmd.hasOption("ed")) {
				final String edOptionValue = cmd.getOptionValue("ed");
				experimentalDesignFile = new File(edOptionValue);
				final File parentFile = experimentalDesignFile.getParentFile();
				if (parentFile == null || !experimentalDesignFile.exists()) {
					experimentalDesignFile = new File(System.getProperty("user.dir") + File.separator + edOptionValue);
				}

				if (!experimentalDesignFile.exists()) {
					errorInParameters(
							"Experimental design file '-ed " + edOptionValue + "' doesn't exist or is not found");
				}
			} else {
				errorInParameters("Experimental design file is missing");
			}
			if (cmd.hasOption("ps")) {
				try {
					psmSelection = PSMSelectionType.valueOf(cmd.getOptionValue("ps"));
				} catch (final Exception e) {
					e.printStackTrace();
					errorInParameters("option 'ps' ('" + options.getOption("ps").getArgName()
							+ "') is not valid. Valid values are " + PSMSelectionType.printValidValues());
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

			final CensusTMT2MSstatsTMT c = new CensusTMT2MSstatsTMT(inputFile, experimentalDesignFile,
					experimentalDesignSeparator, uniquePeptides, useRawIntensity, psmSelection);
			c.run();
		} catch (final ParseException e) {
			errorInParameters(e.getMessage());
		} catch (final Exception e) {
			e.printStackTrace();
			System.err.println(e.getMessage());
		}
	}

	private void run() throws IOException {
		final FileWriter fw = new FileWriter(getOutputFile());
		// header
		printHeader(fw);
		final Map<QuantificationLabel, QuantCondition> conditionsByLabels = getExperimentalDesign()
				.getConditionByLabel();
		final CensusOutParser parser = new CensusOutParser(this.inputFile, conditionsByLabels);
		final Collection<QuantifiedPSMInterface> psms = parser.getPSMMap().values();
		log.info(psms.size() + " psms read");
		// create a map by sequence
		final Map<String, List<QuantifiedPSMInterface>> psmsBySequence = new THashMap<String, List<QuantifiedPSMInterface>>();
		for (final QuantifiedPSMInterface psm : psms) {
			final String seq = psm.getBeforeSeq() + "." + psm.getSequence() + "." + psm.getAfterSeq() + "_"
					+ psm.getChargeState();
			if (!psmsBySequence.containsKey(seq)) {
				psmsBySequence.put(seq, new ArrayList<QuantifiedPSMInterface>());
			}
			psmsBySequence.get(seq).add(psm);
		}
		// iterate over it. Sort first the sequences so that they will appear in order
		// in the output
		final Set<String> seqsSet = psmsBySequence.keySet();
		final List<String> seqList = seqsSet.stream().sorted().collect(Collectors.toCollection(ArrayList::new));
		for (final String seq : seqList) {
			final List<QuantifiedPSMInterface> psmsOfSeq = psmsBySequence.get(seq);

			final List<String> accs = getAccToPrint(psmsOfSeq);
			if (accs.size() > 1 && this.uniquePeptides) {
				log.debug("Skipping PSMs from peptide " + seq + " because not uniqueness");
				continue;
			}
			final int charge = getChargeToPrint(psmsOfSeq);
			final TObjectDoubleMap<QuantificationLabel> intensities = getIntensitiesToPrint(psmsOfSeq,
					this.psmSelection, this.useRawIntensity);
			final QuantifiedPSMInterface firstPSM = psmsOfSeq.get(0);
			final String peptideSequence = firstPSM.getBeforeSeq() + "." + firstPSM.getSequence() + "."
					+ firstPSM.getAfterSeq();
			final String psm = peptideSequence + "_" + charge;
			final String run = firstPSM.getMSRun().getRunId();
			for (final QuantificationLabel label : intensities.keySet().stream().sorted()
					.collect(Collectors.toList())) {
				printOutputLine(fw, label, intensities.get(label), accs.get(0), charge, peptideSequence, psm, run);
			}

		}
		fw.close();
		log.info("File written at " + getOutputFile().getAbsolutePath());
	}

	private void printHeader(FileWriter fw) throws IOException {
		fw.write(
				"ProteinName\tPeptideSequence\tCharge\tPSM\tMixture\tTechRepMixture\tRun\tChannel\tCondition\tBioReplicate\tIntensity\n");

	}

	private TObjectDoubleMap<QuantificationLabel> getIntensitiesToPrint(List<QuantifiedPSMInterface> psms,
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

	private TObjectDoubleMap<QuantificationLabel> getHighestIntensitiesToPrint(List<QuantifiedPSMInterface> psms,
			boolean useRawIntensity) {
		QuantifiedPSMInterface highestPSM = null;
		double highestSum = -Double.MAX_VALUE;
		for (final QuantifiedPSMInterface psm : psms) {
			final TObjectDoubleMap<QuantificationLabel> ret = getIntensitiesFromPSM(psm, useRawIntensity);
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

	private TObjectDoubleMap<QuantificationLabel> getIntensitiesFromPSM(QuantifiedPSMInterface psm,
			boolean useRawIntensity) {
		final TObjectDoubleMap<QuantificationLabel> ret = new TObjectDoubleHashMap<QuantificationLabel>();
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

	private TObjectDoubleMap<QuantificationLabel> getAveragedIntensitiesToPrint(List<QuantifiedPSMInterface> psms,
			boolean useRawIntensity) {
		final Map<QuantificationLabel, TDoubleList> toAverage = new THashMap<QuantificationLabel, TDoubleList>();
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
					final int channel = Integer.valueOf(split[1]);
					final QuantificationLabel label = ExperimentalDesign.getTMT6LabelFromChannel(channel);
					if (!toAverage.containsKey(label)) {
						toAverage.put(label, new TDoubleArrayList());
					}
					toAverage.get(label).add(amount.getValue());
				}
			}
		}
		final TObjectDoubleMap<QuantificationLabel> ret = new TObjectDoubleCustomHashMap<QuantificationLabel>();
		for (final QuantificationLabel label : toAverage.keySet()) {
			ret.put(label, Maths.mean(toAverage.get(label)));
		}
		return ret;
	}

	private TObjectDoubleMap<QuantificationLabel> getSummedIntensitiesToPrint(List<QuantifiedPSMInterface> psms,
			boolean useRawIntensity) {
		final TObjectDoubleMap<QuantificationLabel> ret = new TObjectDoubleCustomHashMap<QuantificationLabel>();
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
					final int channel = Integer.valueOf(split[1]);
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
		final int techRepMixture = experimentalDesign.getTechRepMixtureByRun(run);
		final int bioReplicate = experimentalDesign.getBioReplicate(channel, techRepMixture, mixture);
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

	private static void setupCommandLineOptions() {
		// create Options object
		options = new Options();

		// add t option
		options.addRequiredOption("i", "input", true, "Path to the input file.");
		options.addRequiredOption("ed", "experimental_design", true, "Path to the experimental design file.");
		options.addOption("ps", "psm_selection", true,
				"[OPTIONAL] What to do with multiple PSMs of the same peptide (SUM, AVERAGE or HIGHEST). If not provided, HIGHEST will be choosen.");
		options.addOption("u", "unique", false,
				"[OPTIONAL] Use only unique peptides. If not provided, all peptides will be used.");
		options.addOption("r", "raw", false,
				"[OPTIONAL] Use of raw intensity. If not provided, normalized intensity will be used.");

	}

	private static void errorInParameters(String header) {
		// automatically generate the help statement
		final HelpFormatter formatter = new HelpFormatter();
		if (header == null) {
			formatter.printHelp("CensusTMT2MSstatsTMT -i [input file] -ed [experimental design file]", options);
		} else {
			formatter.printHelp("CensusTMT2MSstatsTMT -i [input file] -ed [experimental design file]",
					"\n\n************\n" + header + "\n************\n\n", options,
					"Contact Salvador Martinez-Bartolome at salvador@scripps.edu for more help");
		}
		System.exit(-1);
	}

}
