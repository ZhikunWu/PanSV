

### install interproscan
```
mamba install bioconda::interproscan
```



### download database
```

######################################
# First time usage please README !!! #
######################################

The databases are huge and consequently not shipped within this installation.
Please download and install the Databases manually by following the commands below:
!!! /!\ Edit the 2 first lines to match the wished version of the DB /!\ !!!

Commands:
=========
# See here for latest db available: https://github.com/ebi-pf-team/interproscan or http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/
# Set versions
version_major=5.59
version_minor=91.0
CONDA_PREFIX=/the/path/to/your/interproscan/conda/env/

# get the md5 of the databases
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${version_major}-${version_minor}/interproscan-${version_major}-${version_minor}-64-bit.tar.gz.md5
# get the databases (with core because much faster to download)
wget http://ftp.ebi.ac.uk/pub/software/unix/iprscan/5/${version_major}-${version_minor}/interproscan-${version_major}-${version_minor}-64-bit.tar.gz
# checksum
md5sum -c interproscan-${version_major}-${version_minor}-64-bit.tar.gz.md5
# untar gz
tar xvzf interproscan-${version_major}-${version_minor}-64-bit.tar.gz
# remove the sample DB bundled by default
rm -rf $CONDA_PREFIX/share/InterProScan/data/
# copy the new db
cp -r interproscan-${version_major}-${version_minor}/data $CONDA_PREFIX/share/InterProScan/

```

/home/wuzhikun/anaconda3/envs/PanSV/share/InterProScan/data




### parameters

```

$ /home/wuzhikun/anaconda3/envs/PanSV/bin/interproscan.sh --help
29/05/2024 13:22:56:997 Welcome to InterProScan-5.59-91.0
29/05/2024 13:22:56:999 Running InterProScan v5 in STANDALONE mode... on Linux
usage: java -XX:+UseParallelGC -XX:ParallelGCThreads=2 -XX:+AggressiveOpts -XX:+UseFastAccessorMethods -Xms128M
            -Xmx2048M -jar interproscan-5.jar


Please give us your feedback by sending an email to

interhelp@ebi.ac.uk

 -appl,--applications <ANALYSES>                           Optional, comma separated list of analyses.  If this option
                                                           is not set, ALL analyses will be run.
 -b,--output-file-base <OUTPUT-FILE-BASE>                  Optional, base output filename (relative or absolute path).
                                                           Note that this option, the --output-dir (-d) option and the
                                                           --outfile (-o) option are mutually exclusive.  The
                                                           appropriate file extension for the output format(s) will be
                                                           appended automatically. By default the input file path/name
                                                           will be used.
 -cpu,--cpu <CPU>                                          Optional, number of cores for inteproscan.
 -d,--output-dir <OUTPUT-DIR>                              Optional, output directory.  Note that this option, the
                                                           --outfile (-o) option and the --output-file-base (-b) option
                                                           are mutually exclusive. The output filename(s) are the same
                                                           as the input filename, with the appropriate file extension(s)
                                                           for the output format(s) appended automatically .
 -dp,--disable-precalc                                     Optional.  Disables use of the precalculated match lookup
                                                           service.  All match calculations will be run locally.
 -dra,--disable-residue-annot                              Optional, excludes sites from the XML, JSON output
 -etra,--enable-tsv-residue-annot                          Optional, includes sites in TSV output
 -exclappl,--excl-applications <EXC-ANALYSES>              Optional, comma separated list of analyses you want to
                                                           exclude.
 -f,--formats <OUTPUT-FORMATS>                             Optional, case-insensitive, comma separated list of output
                                                           formats. Supported formats are TSV, XML, JSON, and GFF3.
                                                           Default for protein sequences are TSV, XML and GFF3, or for
                                                           nucleotide sequences GFF3 and XML.
 -goterms,--goterms                                        Optional, switch on lookup of corresponding Gene Ontology
                                                           annotation (IMPLIES -iprlookup option)
 -help,--help                                              Optional, display help information
 -i,--input <INPUT-FILE-PATH>                              Optional, path to fasta file that should be loaded on Master
                                                           startup. Alternatively, in CONVERT mode, the InterProScan 5
                                                           XML file to convert.
 -incldepappl,--incl-dep-applications <INC-DEP-ANALYSES>   Optional, comma separated list of deprecated analyses that
                                                           you want included.  If this option is not set, deprecated
                                                           analyses will not run.
 -iprlookup,--iprlookup                                    Also include lookup of corresponding InterPro annotation in
                                                           the TSV and GFF3 output formats.
 -ms,--minsize <MINIMUM-SIZE>                              Optional, minimum nucleotide size of ORF to report. Will only
                                                           be considered if n is specified as a sequence type. Please be
                                                           aware of the fact that if you specify a too short value it
                                                           might be that the analysis takes a very long time!
 -o,--outfile <EXPLICIT_OUTPUT_FILENAME>                   Optional explicit output file name (relative or absolute
                                                           path).  Note that this option, the --output-dir (-d) option
                                                           and the --output-file-base (-b) option are mutually
                                                           exclusive. If this option is given, you MUST specify a single
                                                           output format using the -f option.  The output file name will
                                                           not be modified. Note that specifying an output file name
                                                           using this option OVERWRITES ANY EXISTING FILE.
 -pa,--pathways                                            Optional, switch on lookup of corresponding Pathway
                                                           annotation (IMPLIES -iprlookup option)
 -t,--seqtype <SEQUENCE-TYPE>                              Optional, the type of the input sequences (dna/rna (n) or
                                                           protein (p)).  The default sequence type is protein.
 -T,--tempdir <TEMP-DIR>                                   Optional, specify temporary file directory (relative or
                                                           absolute path). The default location is temp/.
 -verbose,--verbose                                        Optional, display more verbose log output
 -version,--version                                        Optional, display version number
 -vl,--verbose-level <VERBOSE-LEVEL>                       Optional, display verbose log output at level specified.
 -vtsv,--output-tsv-version                                Optional, includes a TSV version file along with any TSV
                                                           output (when TSV output requested)
Copyright © EMBL European Bioinformatics Institute, Hinxton, Cambridge, UK. (http://www.ebi.ac.uk) The InterProScan
software itself is provided under the Apache License, Version 2.0 (http://www.apache.org/licenses/LICENSE-2.0.html).
Third party components (e.g. member database binaries and models) are subject to separate licensing - please see the
individual member database websites for details.

Available analyses:
                      TIGRFAM (15.0) : TIGRFAMs are protein families based on hidden Markov models (HMMs).
                       FunFam (4.3.0) : Prediction of functional annotations for novel, uncharacterized sequences.
                         SFLD (4) : SFLD is a database of protein families based on hidden Markov models (HMMs).
                  SUPERFAMILY (1.75) : SUPERFAMILY is a database of structural and functional annotations for all proteins and genomes.
                       Gene3D (4.3.0) : Structural assignment for whole genes and genomes using the CATH domain structure database.
                        Hamap (2021_04) : High-quality Automated and Manual Annotation of Microbial Proteomes.
                        Coils (2.2.1) : Prediction of coiled coil regions in proteins.
              ProSiteProfiles (2022_01) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them.
                        SMART (7.1) : SMART allows the identification and analysis of domain architectures based on hidden Markov models (HMMs). 
                          CDD (3.18) : CDD predicts protein domains and families based on a collection of well-annotated multiple sequence alignment models.
                       PRINTS (42.0) : A compendium of protein fingerprints - a fingerprint is a group of conserved motifs used to characterise a protein family.
                        PIRSR (2021_05) : PIRSR is a database of protein families based on hidden Markov models (HMMs) and Site Rules.
              ProSitePatterns (2022_01) : PROSITE consists of documentation entries describing protein domains, families and functional sites as well as associated patterns and profiles to identify them.
                      AntiFam (7.0) : AntiFam is a resource of profile-HMMs designed to identify spurious protein predictions.
                         Pfam (35.0) : A large collection of protein families, each represented by multiple sequence alignments and hidden Markov models (HMMs).
                   MobiDBLite (2.0) : Prediction of intrinsically disordered regions in proteins.
                        PIRSF (3.10) : The PIRSF concept is used as a guiding principle to provide comprehensive and non-overlapping clustering of UniProtKB sequences into a hierarchical order to reflect their evolutionary relationships.

Deactivated analyses:
                  SignalP_EUK (4.1) : Analysis SignalP_EUK is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
                        TMHMM (2.0c) : Analysis TMHMM is deactivated, because the resources expected at the following paths do not exist: bin/tmhmm/2.0c/decodeanhmm, data/tmhmm/2.0c/TMHMM2.0c.model
                      Phobius (1.01) : Analysis Phobius is deactivated, because the resources expected at the following paths do not exist: bin/phobius/1.01/phobius.pl
        SignalP_GRAM_POSITIVE (4.1) : Analysis SignalP_GRAM_POSITIVE is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp
                      PANTHER (17.0) : Analysis Panther is deactivated, because the resources expected at the following paths do not exist: data/panther/17.0/famhmm/binHmm
        SignalP_GRAM_NEGATIVE (4.1) : Analysis SignalP_GRAM_NEGATIVE is deactivated, because the resources expected at the following paths do not exist: bin/signalp/4.1/signalp



```


