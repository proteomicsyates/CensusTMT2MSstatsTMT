# CensusTMT2MSstatsTMT
Converter from Census TMT output file to the input of MSstatsTMT

Here we developed a command line version and a GUI version. The GUI version is built automatically from the command line version to give graphical support to the options in the command line (implementing class [*CommandLineProgramGuiEnclosable*](https://github.com/proteomicsyates/utilities/blob/master/src/main/java/edu/scripps/yates/utilities/swing/CommandLineProgramGuiEnclosable.java)).

Both versions are available to download at: http://sealion.scripps.edu/Census2MSstatsTMT/

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
