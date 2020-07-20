# CensusTMT2MSstatsTMT
Converter from Census TMT output file to the input of MSstatsTMT

This tool can be executed in command line, or with a graphical interface, when using parameter *-gui* (included in batch files *START_win.bat* and *START_mac_linux.sh*).  
The GUI version is built automatically from the command line version to give graphical support to the options in the command line (implementing class [*CommandLineProgramGuiEnclosable*](https://github.com/proteomicsyates/utilities/blob/master/src/main/java/edu/scripps/yates/utilities/swing/CommandLineProgramGuiEnclosable.java)).

Both versions are available to download at: http://sealion.scripps.edu/CensusTMT2MSstatsTMT/

Command line options:
```
usage: CensusTMT2MSstatsTMT -i [input file] -an [annotation file]

 -an,--annotation <arg>      Path to the experimental design file.
 -d,--decoy <arg>            [OPTIONAL] Remove decoy hits. Decoys hits
                             will have this prefix in their accession
                             number. If not provided, no decoy filtering
                             will be used.
 -i,--input <arg>            Path to the input file.
 -m,--minPeptides <arg>      [OPTIONAL] Minimum number of peptides+charge
                             per protein. If not provided, even proteins
                             with 1 peptide will be quantified
 -ps,--psm_selection <arg>   [OPTIONAL] What to do with multiple PSMs of
                             the same peptide (SUM, AVERAGE or HIGHEST).
                             If not provided, HIGHEST will be choosen.
 -r,--raw                    [OPTIONAL] Use of raw intensity. If not
                             provided, normalized intensity will be used.
 -u,--unique                 [OPTIONAL] Use only unique peptides. If not
                             provided, all peptides will be used.
Contact Salvador Martinez-Bartolome at salvador at scripps.edu for more help
```
To know more about the annotation file, go to http://msstats.org/msstatstmt/  

The annotation file is a COMMA-SEPARATED file (CSV) containing the information about the experimental design.  
The file should have the following columns:

| Column | Explanation | 
| ------ | :---------- |
| Run | MS run ID. It should be the same as Spectrum.File info in output of Proteome Discoverer.|
| Channel | Labeling information (126, … 131). It should be consistent with the channel columns in output of Proteome Discoverer.|
| Condition | Condition (ex. Healthy, Cancer, Time0). If the channal doesn’t have sample, please add Empty under Condition.|
| Mixture | Mixture of samples labeled with different TMT reagents, which can be analyzed in a single mass spectrometry experiment.|
| TechRepMixture | Technical replicate of one mixture. One mixture may have multiple technical replicates. For example, if TechRepMixture = 1, 2 are the two technical replicates of one mixture, then they should match with same Mixture value.|
| Fraction | Fraction ID. One technical replicate of one mixture may be fractionated into multiple fractions to increase the analytical depth. Then one technical replicate of one mixture should correspond to multuple fractions. For example, if Fraction = 1, 2, 3 are three fractions of the first technical replicate of one TMT mixture of biological subjects, then they should have same TechRepMixture and Mixture value.|
| BioReplicate | Unique ID for biological subject. If the channal doesn’t have sample, please add Empty under BioReplicate.|
