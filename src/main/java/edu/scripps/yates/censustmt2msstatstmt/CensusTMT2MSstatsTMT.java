package edu.scripps.yates.censustmt2msstatstmt;

import java.io.File;
import java.io.FileNotFoundException;
import java.io.IOException;
import java.util.ArrayList;
import java.util.Collection;
import java.util.EnumMap;
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
import edu.scripps.yates.census.read.QuantParserException;
import edu.scripps.yates.census.read.WrongTMTLabels;
import edu.scripps.yates.census.read.model.QuantAmount;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedPSMInterface;
import edu.scripps.yates.census.read.model.interfaces.QuantifiedProteinInterface;
import edu.scripps.yates.census.read.util.QuantificationLabel;
import edu.scripps.yates.censustmt2msstatstmt.singletons.FilesManager;
import edu.scripps.yates.censustmt2msstatstmt.singletons.GeneMapper;
import edu.scripps.yates.censustmt2msstatstmt.singletons.IonFilter;
import edu.scripps.yates.censustmt2msstatstmt.singletons.LuciphorIntegrator;
import edu.scripps.yates.censustmt2msstatstmt.singletons.PTMList;
import edu.scripps.yates.censustmt2msstatstmt.singletons.ProteinSequences;
import edu.scripps.yates.censustmt2msstatstmt.singletons.SPCFilter;
import edu.scripps.yates.censustmt2msstatstmt.singletons.UPLR;
import edu.scripps.yates.censustmt2msstatstmt.util.DataUtil;
import edu.scripps.yates.utilities.appversion.AppVersion;
import edu.scripps.yates.utilities.dates.DatesUtil;
import edu.scripps.yates.utilities.grouping.GroupableProtein;
import edu.scripps.yates.utilities.grouping.PAnalyzer;
import edu.scripps.yates.utilities.grouping.ProteinGroup;
import edu.scripps.yates.utilities.properties.PropertiesUtil;
import edu.scripps.yates.utilities.proteomicsmodel.PTM;
import edu.scripps.yates.utilities.proteomicsmodel.enums.AmountType;
import edu.scripps.yates.utilities.sequence.PositionInProtein;
import edu.scripps.yates.utilities.strings.StringUtils;
import edu.scripps.yates.utilities.swing.CommandLineProgramGuiEnclosable;
import edu.scripps.yates.utilities.swing.DoNotInvokeRunMethod;
import edu.scripps.yates.utilities.swing.SomeErrorInParametersOcurred;
import gnu.trove.map.hash.THashMap;
import gnu.trove.set.hash.THashSet;

public class CensusTMT2MSstatsTMT extends CommandLineProgramGuiEnclosable {
	private final static Logger log = Logger.getLogger(CensusTMT2MSstatsTMT.class);

	private static AppVersion version;
	private boolean uniquePeptides;
	private boolean useRawIntensity;
	private PSMAggregationType psmSelection;
	private File experimentalDesignFile;
	private String decoyPrefix;
	private Integer minSPC;
	private boolean msstatsoutput;
	private File fastaFile;
	private Float luciphorLocalFDRThreshold;
	private Float luciphorGlobalFDRThreshold;
	private String outputPrefix;
	private DecoyFilter decoyFilter;

	private CensusOutParser parser;

	private edu.scripps.yates.censustmt2msstatstmt.SPC_FILTER_TYPE spcFilterType;

	private File luciphorFile;

	private Float minPurity;

	private boolean createExcelFile;

	private boolean simplifyProteinGroups;

	private Integer minIONs;

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

	private void runConversionToMSstatsTMT() throws Exception {

		if (decoyPrefix != null) {
			decoyFilter = new DecoyFilter(decoyPrefix);
		}
		ExperimentalDesign ed = null;
		try {

			System.out.println("Starting conversion of file " + FilesManager.getInstance().getInputFileNamesString());

			if (this.msstatsoutput) {
				ed = new ExperimentalDesign(experimentalDesignFile, ",");
			}
			// 1- get conditions and labels
			final Map<QuantificationLabel, QuantCondition>[] conditionsByLabelsList = getConditionsByLabelsArray(
					FilesManager.getInstance().getInputFiles(), ed);

			// check whether multiple mixtures have different experimental designs or not
			final boolean multipleMixturesWithDifferentExperimentalDesigns = checkMixtureExperimentalDesigns(
					conditionsByLabelsList);

			// 2- read input files with parser

			parser = new CensusOutParser(FilesManager.getInstance().getInputFiles().toArray(new File[0]),
					conditionsByLabelsList, null, null);
			parser.setChargeSensible(true);
			parser.setDistinguishModifiedSequences(true);
			parser.setDecoyPattern(this.decoyPrefix);

			List<QuantificationLabel> labels = null;
			if (parser.isTMT10().values().iterator().next()) {
				System.out.println("10-Plex detected.");
				labels = QuantificationLabel.getTMT10PlexLabels();
			} else if (parser.isTMT6().values().iterator().next()) {
				System.out.println("6-Plex detected.");
				labels = QuantificationLabel.getTMT6PlexLabels();
			} else if (parser.isTMT11().values().iterator().next()) {
				System.out.println("11-Plex detected.");
				labels = QuantificationLabel.getTMT11PlexLabels();
			} else if (parser.isTMT4().values().iterator().next()) {
				System.out.println("4-Plex detected.");
				labels = QuantificationLabel.getTMT4PlexLabels();
			}
			String plural = "";
			if (FilesManager.getInstance().getInputFiles().size() > 1) {
				plural = "s";
			}
			System.out.println("Reading input file" + plural + "...");
			// 3- map genes to proteins from protein descriptions in proteins
			final Collection<QuantifiedProteinInterface> proteins = parser.getProteinMap().values();
			GeneMapper.getInstance().mapGenesFromProteinDescriptionsInInputFile(proteins);

			// get the PSMs
			List<QuantifiedPSMInterface> psms = new ArrayList<QuantifiedPSMInterface>();
			psms.addAll(parser.getPSMMap().values());
			System.out.println(psms.size() + " PSMs read from input file.");

			// setting missing intensities to 1
			setMissingIntensityValuesTo1(psms, labels);

			// 4- merge psms with luciphor ones
			if (luciphorFile != null) {
				final LuciphorIntegrator luciphorIntegrator = LuciphorIntegrator.getInstance(luciphorFile);
				luciphorIntegrator.mergeWithLuciphor(psms, luciphorLocalFDRThreshold, luciphorGlobalFDRThreshold);
			}

			// 5- group proteins so that we can remove non conclusive ones
			final List<GroupableProtein> proteinsToGroup = new ArrayList<GroupableProtein>();
			proteinsToGroup.addAll(proteins);
			final boolean separateNonConclusiveProteins = true;
			final PAnalyzer panalyzer = new PAnalyzer(separateNonConclusiveProteins);
			final List<ProteinGroup> groups = panalyzer.run(proteinsToGroup);

			// 6- retrieve annotations from uniprot
			final long t1 = System.currentTimeMillis();
			System.out.println("Retrieving information from UniprotKB for " + proteins.size()
					+ " proteins. This can take a few minutes (2-3) for the first time. Next time information will be ready to use in few seconds. Please wait...");
			UPLR.getInstance().getAnnotatedProteins(null, parser.getProteinMap().keySet());
			final long t2 = System.currentTimeMillis() - t1;
			System.out.println(
					"Information retrieved from UniprotKB in " + DatesUtil.getDescriptiveTimeFromMillisecs(t2));

			// 7- if ptms is not empty, we have to get the protein sequences
			if (!PTMList.getInstance().isEmpty() || (this.fastaFile != null && this.fastaFile.exists())) {
				ProteinSequences.getInstance().getProteinSequences(parser.getProteinMap().keySet(), this.fastaFile);
			}

			final Set<String> validPSMKeys = new THashSet<String>();
			int discardedByTMTPurity = 0;
			int discardedByPTMs = 0;
			int discardedAsDecoy = 0;
			int discardedByUniqueness = 0;

			// 8- Iterate over psms
			final Set<String> psmsWronglymapped = new THashSet<String>();
			for (final QuantifiedPSMInterface psm : psms) {
				// 8.1- check whether the peptides are wrongly mapped in the input file
				for (final QuantifiedProteinInterface protein : psm.getQuantifiedProteins()) {
					final String acc = protein.getAccession();
					final List<PositionInProtein> startingPositionsInProtein = psm.getStartingPositionsInProtein(acc,
							UPLR.getInstance(), ProteinSequences.getInstance());
					if (startingPositionsInProtein.isEmpty()) {
						psmsWronglymapped.add(psm.getScanNumber() + "\t" + psm.getFullSequence() + "\t" + acc + "\t"
								+ protein.getDescription());
					}
				}
				// 8.2- check decoy
				if (decoyFilter != null && decoyFilter.isDecoy(psm)) {
					psm.setDiscarded(true);
					discardedAsDecoy++;
					continue;
				}
				// 8.3- TMT purity
				if (minPurity != null && !DataUtil.isTMTPurityValid(psm, minPurity)) {
					psm.setDiscarded(true);
					discardedByTMTPurity++;
					continue;
				}

				// 8.4 if ptmList is not empty, check if contains that PTM
				if (!PTMList.getInstance().isEmpty()) {
					if (psm.getPTMs().isEmpty()) {
						psm.setDiscarded(true);
						discardedByPTMs++;
						continue;
					}
					boolean ptmFound = false;
					for (final float deltaMass : PTMList.getInstance()) {
						for (final PTM ptm : psm.getPTMs()) {
							if (Math.abs(ptm.getMassShift() - deltaMass) < PTMList.MASS_PRECISION) {
								ptmFound = true;
								break;
							}
						}
					}
					if (!ptmFound) {
						psm.setDiscarded(true);
						discardedByPTMs++;
						continue;
					}
				}
				// 8.5- discard non-unique psms
				if (uniquePeptides) {
					if (psm.getQuantifiedProteins().size() > 1) {
						psm.setDiscarded(true);
						discardedByUniqueness++;
					}
				}
				// if we are here, the psm is fine.
				validPSMKeys.add(psm.getKey());
			}

			if (!psmsWronglymapped.isEmpty()) {
				String message = "WARNING!! There are " + psmsWronglymapped.size()
						+ " that are wrongly mapped to the proteins in your input file!";
				if (!PTMList.getInstance().isEmpty()) {
					message += " For the PTM-based analysis, these peptides were not able to be mapped to any site in the protein.";
				}
				if (fastaFile != null) {
					message += " This means that either you used the wrong fasta file (not the one used in the search of your data), or there is a problem in your input files";
				} else {
					message += " This means that sequences retrieved in UniprotKB have changed compared to the ones you used in your search. Try to use the fasta file you actually used for search if you want to avoid this problem.";
				}
				System.out.println(message);
				FilesManager.getInstance().printList(psmsWronglymapped, "wrong_mapped");
			}
			if (discardedByUniqueness > 0) {
				System.out.println(discardedByUniqueness + " PSMs discarded because they are not unique");
			}
			if (discardedAsDecoy > 0) {
				System.out.println(discardedAsDecoy + " PSMs discarded as decoys");
			}
			if (minPurity != null) {
				System.out.println(discardedByTMTPurity + " PSMs discarded by TMT purity threshold " + minPurity);
			}
			if (!PTMList.getInstance().isEmpty()) {
				System.out.println(discardedByPTMs + " PSMs discarded for not having any PTM with delta mass: ["
						+ StringUtils.getSeparatedValueStringFromChars(PTMList.getInstance().toArray(), ",") + "]");
			}
			System.out.println(psms.size() - validPSMKeys.size() + " PSMs discarded in total");
			System.out.println(validPSMKeys.size() + " PSMs pass the filters (if any)");

			// 9- SPC filter
			if (minSPC != null && minSPC > 1) {
				psms = SPCFilter.applySPCFilter(psms, this.minSPC, this.spcFilterType);
			}

			// 9.5- SPC filter
			if (minIONs != null && minIONs > 1) {
				if (PTMList.getInstance().isEmpty()) {
					psms = IonFilter.applyMinIonsPerProteinFilter(groups, this.minIONs);
				} else {
					psms = IonFilter.applyMinIonsPerProteinSite(psms, minIONs);
				}
			}
			// check whether some psm is discarded
			psms = psms.stream().filter(psm -> !psm.isDiscarded()).collect(Collectors.toList());

			// 10- print PSM level file
			FilesManager.getInstance().printPSMLevelFile(psms, parser, labels);

			// 11- print MSstatsTMT file
			if (this.msstatsoutput) {
				final boolean aggregateByPTMs = !PTMList.getInstance().isEmpty();
				FilesManager.getInstance().printMSstatsTMTFFile(psms, ed, this.useRawIntensity, this.psmSelection,
						aggregateByPTMs, simplifyProteinGroups);
			}

			// 12- print peptide level file
			FilesManager.getInstance().printPeptideLevelFile(parser.getPeptideMap().values(), psmSelection, labels);

			//
			if (createExcelFile) {
				FilesManager.getInstance().generateExcelFile();
			}
			//

		} catch (final Exception e) {
			if (e instanceof WrongTMTLabels) {
				throw new WrongTMTLabels(
						"Error in annotation file. The TMT channels in annotation file don't match with TMT channels in input file: "
								+ e.getMessage());
			} else {
				throw e;
			}
		}
	}

	private boolean checkMixtureExperimentalDesigns(Map<QuantificationLabel, QuantCondition>[] conditionsByLabelsList) {
		// first we check if the plexes are equal, if not we throw exception
		int n = -1;
		for (final Map<QuantificationLabel, QuantCondition> map : conditionsByLabelsList) {
			if (n == -1) {
				n = map.size();
			} else {
				if (n != map.size()) {
					throw new IllegalArgumentException("Different Plexes detected in input files");
				}
			}
		}

		return true;
	}

	private void setMissingIntensityValuesTo1(List<QuantifiedPSMInterface> psms, List<QuantificationLabel> labels) {
		int intensitiesChanged = 0;
		for (final QuantifiedPSMInterface psm : psms) {

			final EnumMap<QuantificationLabel, Set<QuantAmount>> amountsByLabels = DataUtil
					.getQuantAmountsByLabels(psm);
			for (final QuantificationLabel label : labels) {
				final Set<QuantAmount> quantAmounts = amountsByLabels.get(label);
				if (quantAmounts != null && !quantAmounts.isEmpty()) {
					boolean someChanged = false;
					for (final QuantAmount amount : quantAmounts) {
						if (Double.compare(amount.getValue(), 0.0) == 0) {
							amount.setValue(1.0);
							someChanged = true;
						}
					}
					if (someChanged) {
						intensitiesChanged++;
					}
				} else {
					// create empty amounts
					final QuantAmount quantAmount = new QuantAmount(1.0, AmountType.INTENSITY, null, label);
					psm.addAmount(quantAmount);
					final QuantAmount normQuantAmount = new QuantAmount(1.0, AmountType.NORMALIZED_INTENSITY, null,
							label);
					psm.addAmount(normQuantAmount);
					intensitiesChanged++;
				}
			}
		}
		if (intensitiesChanged > 0) {
			System.out.println(intensitiesChanged + "  missing intensities were set to 1.0");
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
	 * @throws IOException
	 */
	private Map<QuantificationLabel, QuantCondition>[] getConditionsByLabelsArray(List<File> inputFiles2,
			ExperimentalDesign experimentalDesign) throws IOException {
		final Map<QuantificationLabel, QuantCondition>[] ret = new Map[inputFiles2.size()];
		if (experimentalDesign != null) {
			if (inputFiles2.size() == 1) {
				ret[0] = experimentalDesign.getMixtures().iterator().next().getConditionsByLabels();
				return ret;
			}
			int index = 0;
			for (final File inputFile : inputFiles2) {

				try {
					final CensusOutParser parser = new CensusOutParser(inputFile, null);
					parser.setIgnoreACCFormat(true);
					parser.setIgnoreTaxonomies(true);
					parser.setChargeSensible(true);
					parser.setDistinguishModifiedSequences(true);
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
		} else {

			final CensusOutParser parser = new CensusOutParser(inputFiles2, null, null, null);
			parser.setIgnoreACCFormat(true);
			parser.setIgnoreTaxonomies(true);
			parser.setChargeSensible(true);
			parser.setDistinguishModifiedSequences(true);
			int index = 0;
			for (final File inputFile : inputFiles2) {
				try {
					final Map<QuantificationLabel, QuantCondition> conditionsByLabels = new THashMap<QuantificationLabel, QuantCondition>();
					final boolean istmt4 = parser.isTMT4Files().get(inputFile);
					final boolean istmt6 = parser.isTMT6Files().get(inputFile);
					final boolean istmt10 = parser.isTMT10Files().get(inputFile);
					final boolean istmt11 = parser.isTMT11Files().get(inputFile);
					if (istmt4) {
						final List<QuantificationLabel> tmt4PlexLabels = QuantificationLabel.getTMT4PlexLabels();
						int cond = 1;
						for (final QuantificationLabel quantificationLabel : tmt4PlexLabels) {
							conditionsByLabels.put(quantificationLabel, new QuantCondition("cond" + cond));
							cond++;
						}
					} else if (istmt6) {
						final List<QuantificationLabel> tmt6PlexLabels = QuantificationLabel.getTMT6PlexLabels();
						int cond = 1;
						for (final QuantificationLabel quantificationLabel : tmt6PlexLabels) {
							conditionsByLabels.put(quantificationLabel, new QuantCondition("cond" + cond));
							cond++;
						}
					} else if (istmt10) {
						final List<QuantificationLabel> tmt10PlexLabels = QuantificationLabel.getTMT10PlexLabels();
						int cond = 1;
						for (final QuantificationLabel quantificationLabel : tmt10PlexLabels) {
							conditionsByLabels.put(quantificationLabel, new QuantCondition("cond" + cond));
							cond++;
						}
					} else if (istmt11) {
						final List<QuantificationLabel> tmt11PlexLabels = QuantificationLabel.getTMT11PlexLabels();
						int cond = 1;
						for (final QuantificationLabel quantificationLabel : tmt11PlexLabels) {
							conditionsByLabels.put(quantificationLabel, new QuantCondition("cond" + cond));
							cond++;
						}
					} else {
						throw new IllegalArgumentException(
								"The file is not detected as any TMT type (TMT6, TMT10 or TMT11");
					}
					ret[index] = conditionsByLabels;

				} catch (final FileNotFoundException e) {
				} finally {
					index++;
				}
			}

		}

		return ret;

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
		final Option optionInputFile = new Option("input", true,
				"Path to the input file(s). It can refer to multiple files by using wildcards ('*'), i.e: '/path/to/my/files/census*.out'");
		optionInputFile.setRequired(true);
		options.add(optionInputFile);

		// in case of wanting MSstatsTMT output
		final Option optionMSstats = new Option("msstats", false, "Create output to use in MSstatsTMT.");
		options.add(optionMSstats);

		final Option optionAnnotation = new Option("annotation", true,
				"Path to the experimental design file. Option required if 'msstats' is selected. "
						+ "Wildcards ('*') can be used, but only when refering to 1 file, i.e: '/path/to/my/files/annotation*.txt");
		options.add(optionAnnotation);

		final Option optionPTM = new Option("ptms", true,
				"Comma separated list of masses corresponding to the PTMs that we want to use to aggregate quant values (instead of peptides)");
		options.add(optionPTM);

		final Option optionFasta = new Option("fasta", true,
				"Full path to a fasta file used for getting the protein sequences for a ptm-based analysis. "
						+ "This option will be ignored if 'ptms' is empty. "
						+ "If not provided, or if proteins in input file are not found in this fasta file, the tool will try to get the protein sequence from UniprotKB in case of having Uniprot accessions."
						+ " Wildcards ('*') can be used, but only when refering to 1 file, i.e: '/path/to/my/files/*.fasta");
		options.add(optionFasta);

		final Option optionLuciphor = new Option("luciphor", true,
				"Full path to the luciphor (PTM localization algorithm) results file. "
						+ "Luciphor results will be merged with input data after applying some FDR cutoffs."
						+ " Wildcards ('*') can be used, but only when refering to 1 file, i.e: '/path/to/my/files/luciphor*.txt");
		optionLuciphor.setRequired(false);
		options.add(optionLuciphor);

		final Option optionLFDR = new Option("luc_lfdr", true,
				"In case of providing a luciphor file, this paremeter specifies the local FDR threshold to apply. "
						+ "All Luciphor results above that threshold will be ignored.");
		options.add(optionLFDR);

		final Option optionGFDR = new Option("luc_gfdr", true,
				"In case of providing a luciphor file, this paremeter specifies the global FDR threshold to apply. "
						+ "All Luciphor results above that threshold will be ignored.");
		options.add(optionGFDR);

		final Option optionPSMAggregation = new Option("aggregation", true,
				"How to aggregate PSM intensities to peptide level. Valid values are: ["
						+ PSMAggregationType.printValidValues() + "]. If not provided, " + PSMAggregationType.HIGHEST
						+ " will be choosen.");
		options.add(optionPSMAggregation);

		final Option optionUnique = new Option("unique", false, "If selected, it will discard any non-unique peptide.");
		options.add(optionUnique);

		final Option optionRaw = new Option("raw", false,
				"Use of raw intensity. If not provided, normalized intensity will be used.");
		options.add(optionRaw);

		final Option optionDecoy = new Option("decoy", true,
				"Remove decoy proteins. Decoys proteins will have this prefix in their accession number. If not provided, no decoy filtering will be used.");
		options.add(optionDecoy);

		final Option optionSPC = new Option("min_spc", true,
				"Minimum number of spectral counts. See 'spc_type' option to specify the behaviour of this filter.");
		options.add(optionSPC);

		final Option optionSPCType = new Option("spc_type", true,
				"This parameter specifies at which level to apply the minimum number of spectral counts filter. "
						+ "Valid values are: " + SPC_FILTER_TYPE.SEQUENCE.name() + " (peptide sequence with PTMs), "
						+ SPC_FILTER_TYPE.ION.name() + " (peptide sequence with PTMs plus charge) or "
						+ SPC_FILTER_TYPE.PROTEIN_SITE.name()
						+ " (protein site in case of being a PTM-based analysis). " + "If not provided, "
						+ SPC_FILTER_TYPE.ION.name()
						+ " will be used. If 'min_spc' is empty, this parameter will be ignored.");
		options.add(optionSPCType);

		final Option optionMinION = new Option("min_ions", true,
				"Minimum number of ions (sequence+charge) per protein (or PTM site in case of aggregating by PTMs).");
		options.add(optionMinION);

		final Option optionTMTPurity = new Option("tmt_purity", true, "Minimum TMT purity value per PSM.");
		options.add(optionTMTPurity);

		final Option optionOutputPrefix = new Option("output_prefix", true,
				"Prefix added to all the output files generated. If not provided, the name of the first input input will be used as prefix.");
		options.add(optionOutputPrefix);

		final Option optionCreateExcel = new Option("excel", false,
				"If selected, an Excel file will be created compiling all the other text output files in different sheets.");
		options.add(optionCreateExcel);

		final Option simplyProteinGroups = new Option("simplify_protein_groups", false,
				"If selected, protein accessions of a group will be simplified (if possible) to show only protein accession from UniprotKB SwissProt, ignoring accession from UniprotKB TrEmBML.");
		options.add(simplyProteinGroups);
		return options;
	}

	@Override
	public void run() {
		try {
			runConversionToMSstatsTMT();
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
		FilesManager.clearStaticData();
		LuciphorIntegrator.clearDataStatic();
		PTMList.getInstance().clearData();
		experimentalDesignFile = null;
		uniquePeptides = false;
		useRawIntensity = false;
		decoyPrefix = null;
		psmSelection = PSMAggregationType.HIGHEST;

		fastaFile = null;

		List<File> inputFiles = null;
		if (cmd.hasOption("input")) {
			inputFiles = new ArrayList<File>();
			final String iOptionValue = cmd.getOptionValue("input");
			File inputFile = new File(iOptionValue);
			final File parentFile = inputFile.getParentFile();
			if (parentFile == null && !inputFile.exists()) {
				inputFile = new File(System.getProperty("user.dir") + File.separator + iOptionValue);
			}

			if (!inputFile.exists()) {
				// check whether is using a wildcard in the name
				if (iOptionValue.matches(".*\\*.*")) {
					inputFiles = FilesManager.getMultipleFiles(parentFile, iOptionValue);
					if (inputFiles.isEmpty()) {
						errorInParameters("Input files '-input " + iOptionValue + "' don't exist or are not found");
					}
				} else {
					errorInParameters("Input file '-input " + iOptionValue + "' doesn't exist or is not found");
				}
			} else {
				inputFiles.add(inputFile);
			}
		} else {
			errorInParameters("Input file is missing");
		}
		String plural = "";
		if (inputFiles.size() > 1) {
			plural += "s";
		}
		System.out.println("Input file" + plural + ":");
		int i = 1;
		for (final File file : inputFiles) {
			System.out.println(i++ + "- " + file.getAbsolutePath() + "'");
		}
		// output prefix
		this.outputPrefix = null;
		if (cmd.hasOption("output_prefix")) {
			this.outputPrefix = cmd.getOptionValue("output_prefix").trim();
			if ("".equals(this.outputPrefix)) {
				this.outputPrefix = null;
			}
		}
		if (outputPrefix == null) {
			this.outputPrefix = FilenameUtils.getBaseName(inputFiles.get(0).getAbsolutePath());
			System.out.println("Option '-output_prefix' not provided. Using prefix '" + this.outputPrefix
					+ "' for all output files.");
		}
		FilesManager.getNewInstance(inputFiles, outputPrefix);
		//
		this.msstatsoutput = false;
		if (cmd.hasOption("msstats")) {
			msstatsoutput = true;
		}
		//
		if (this.msstatsoutput) {
			if (cmd.hasOption("annotation")) {
				final String anOptionValue = cmd.getOptionValue("annotation");
				experimentalDesignFile = new File(anOptionValue);
				final File parentFile = experimentalDesignFile.getParentFile();
				if (parentFile == null || !experimentalDesignFile.exists()) {
					experimentalDesignFile = new File(System.getProperty("user.dir") + File.separator + anOptionValue);
				}

				if (!experimentalDesignFile.exists()) {
					errorInParameters("Annotation file (experimental design) '-annotation " + anOptionValue
							+ "' doesn't exist or is not found");
				}
			} else {
				errorInParameters(
						"Experimental design file is missing and it is required if option -msstats is selected");
			}
		} else {
			if (cmd.hasOption("annotation")) {
				System.out.println("Option -annotation ignored");
			}
		}
		//

		if (cmd.hasOption("ptms")) {
			final String ptmsString = cmd.getOptionValue("ptms").trim();
			try {
				final List<Float> ptms = parsePTMsString(ptmsString);
				PTMList.getInstance().addAll(ptms);
			} catch (final NumberFormatException e) {
				errorInParameters("Error in '-ptms' option. Only a comma separated list of real numbers is valid");
			}
		}
//		if (!PTMList.getInstance().isEmpty()) {
		if (cmd.hasOption("fasta")) {
			final String fastaFilePath = cmd.getOptionValue("fasta").trim();
			fastaFile = new File(fastaFilePath);
			final File parentFile = fastaFile.getParentFile();
			if (parentFile == null && !fastaFile.exists()) {
				fastaFile = new File(System.getProperty("user.dir") + File.separator + fastaFilePath);
			}
			if (!fastaFile.exists()) {
				if (fastaFilePath.matches(".*\\*.*")) {
					final List<File> fastaFiles = FilesManager.getMultipleFiles(parentFile, fastaFilePath);
					if (fastaFiles.isEmpty()) {
						errorInParameters("Fasta file '-fasta " + fastaFilePath + "' don't exist or is not found");
					} else if (fastaFiles.size() > 1) {
						errorInParameters(
								"Option '-fasta " + fastaFilePath + "' refers to more than one file in the folder");
					} else {
						this.fastaFile = fastaFiles.get(0);
					}
				} else {
					errorInParameters("Fasta file '" + fastaFile.getAbsolutePath() + "' is not found");
				}
			}
		} else {
			if (!PTMList.getInstance().isEmpty()) {
				System.out.println(
						"Fasta file not provided when using '-ptms' option. Protein sequences will be retrieved from Uniprot if possible...");
			}
		}
//		} else {
//			System.out.println("Option -fasta ignored");
//		}
		//
		if (cmd.hasOption("aggregation")) {
			try {
				psmSelection = PSMAggregationType.valueOf(cmd.getOptionValue("aggregation"));
			} catch (final Exception e) {
				e.printStackTrace();
				errorInParameters("option 'aggregation' ('" + cmd.getOptionValue("aggregation")
						+ "') is not valid. Valid values are " + PSMAggregationType.printValidValues());
			}

		}
		log.info("Using ps option " + psmSelection);
		this.luciphorFile = null;
		//
		if (cmd.hasOption("luciphor")) {
			final String luciphorFilePath = cmd.getOptionValue("luciphor");
			luciphorFile = new File(luciphorFilePath);
			final File parentFile = luciphorFile.getParentFile();
			if (parentFile == null && !luciphorFile.exists()) {
				luciphorFile = new File(System.getProperty("user.dir") + File.separator + luciphorFilePath);
			}
			if (!luciphorFile.exists()) {
				if (luciphorFilePath.matches(".*\\*.*")) {
					final List<File> luciphorFiles = FilesManager.getMultipleFiles(parentFile, luciphorFilePath);
					if (luciphorFiles.isEmpty()) {
						errorInParameters(
								"Luciphor file '-luciphor " + luciphorFilePath + "' don't exist or is not found");
					} else if (luciphorFiles.size() > 1) {
						errorInParameters("Option '-luciphor " + luciphorFilePath
								+ "' refers to more than one file in the folder");
					} else {
						luciphorFile = luciphorFiles.get(0);
					}
				} else {
					errorInParameters("Luciphor file '" + luciphorFile.getAbsolutePath() + "' not found");
				}
			}
		}
		//
		if (luciphorFile != null && luciphorFile.exists()) {

			if (cmd.hasOption("luc_lfdr")) {
				try {
					this.luciphorLocalFDRThreshold = Float.valueOf(cmd.getOptionValue("luc_lfdr"));
				} catch (final NumberFormatException e) {
					errorInParameters("Local FDR threshold has to be a real number");
				}
			}
			//
			if (cmd.hasOption("luc_gfdr")) {
				try {
					this.luciphorGlobalFDRThreshold = Float.valueOf(cmd.getOptionValue("luc_gfdr"));
				} catch (final NumberFormatException e) {
					errorInParameters("Global FDR threshold has to be a real number");
				}
			}
		} else {
			if (cmd.hasOption("luc_lfdr")) {

				System.out.println("Option '-luc_lfdr' ignored.");
			}
			if (cmd.hasOption("luc_gfdr")) {
				System.out.println("Option '-luc_gfdr' ignored.");
			}

		}
		//
		if (cmd.hasOption("unique")) {
			uniquePeptides = true;
			log.info("Using only unique peptides.");
		} else {
			log.info("Using all peptides regardless of their uniqueness.");
		}
		//
		if (cmd.hasOption("raw")) {
			useRawIntensity = true;
			log.info("Using raw intensities.");
		} else {
			log.info("Using normalized intensities.");
		}

		if (cmd.hasOption("decoy")) {
			decoyPrefix = cmd.getOptionValue("decoy");
			if ("".equals(decoyPrefix)) {
				decoyPrefix = null;
			}
			log.info("Using d option for decoy prefix='" + decoyPrefix + "'");
		}
		minSPC = null;
		if (cmd.hasOption("min_spc")) {
			try {
				minSPC = Integer.valueOf(cmd.getOptionValue("min_spc"));
				if (minSPC < 0) {
					throw new NumberFormatException();
				}
			} catch (final NumberFormatException e) {
				errorInParameters("Option '-min_spc' should have a positive integer");
			}
		}
		if (cmd.hasOption("spc_type")) {
			spcFilterType = SPC_FILTER_TYPE.getValueOf(cmd.getOptionValue("spc_type"));
			if (spcFilterType == null) {
				errorInParameters("Spectral count filter type can only be one of these values: "
						+ SPC_FILTER_TYPE.getValuesString(", "));
			}

		} else {
			if (this.minSPC != null) {
				errorInParameters("Option '-spc_type' is missing");
			} else {
				this.spcFilterType = SPC_FILTER_TYPE.ION;
			}
		}
		//
		minIONs = null;
		if (cmd.hasOption("min_ions")) {
			try {
				minIONs = Integer.valueOf(cmd.getOptionValue("min_ions"));
				if (minIONs < 0) {
					throw new NumberFormatException();
				}
			} catch (final NumberFormatException e) {
				errorInParameters("Option '-min_ions' should have a positive integer");
			}
		}
		//
		if (cmd.hasOption("pur")) {
			try {
				minPurity = Float.valueOf(cmd.getOptionValue("pur"));

			} catch (final NumberFormatException e) {
				errorInParameters("Minimum TMT purity has to be a number");
			}
		} else {
			minPurity = null;
		}
		this.createExcelFile = false;
		if (cmd.hasOption("excel")) {
			this.createExcelFile = true;
		}
		this.simplifyProteinGroups = false;
		if (cmd.hasOption("simplify_protein_groups")) {
			this.simplifyProteinGroups = true;
		}
	}

	private List<Float> parsePTMsString(String ptmsString) {
		final List<Float> ret = new ArrayList<Float>();
		if (ptmsString.contains(",")) {
			final String[] split = ptmsString.split(",");
			for (final String ptmString : split) {
				ret.add(Float.valueOf(ptmString.trim()));
			}
		} else {
			ret.add(Float.valueOf(ptmsString.trim()));
		}
		if (ret.isEmpty()) {
			return null;
		}
		return ret;
	}

	@Override
	public String printCommandLineSintax() {
		return "CensusTMT2MSstatsTMT -i [input file] -an [annotation file]";
	}

}
