
## [EDTA](https://github.com/oushujun/EDTA): The Extensive de novo TE Annotator


### install

```
git clone https://github.com/oushujun/EDTA.git
conda env create -f EDTA.yml

```


### parameters

```
(EDTA) wuzhikun@mu03 16:58:51 ^_^ /home/wuzhikun/software/EDTA 
$ perl /home/wuzhikun/software/EDTA/EDTA.pl

########################################################
##### Extensive de-novo TE Annotator (EDTA) v2.1.0  ####
##### Shujun Ou (shujun.ou.1@gmail.com)             ####
########################################################



At least 1 parameter is required:
1) Input fasta file: --genome

This is the Extensive de-novo TE Annotator that generates a high-quality
structure-based TE library. Usage:

perl EDTA.pl [options]
    --genome [File]     The genome FASTA file. Required.
    --species [Rice|Maize|others]   Specify the species for identification of TIR
                    candidates. Default: others
    --step [all|filter|final|anno]  Specify which steps you want to run EDTA.
                    all: run the entire pipeline (default)
                    filter: start from raw TEs to the end.
                    final: start from filtered TEs to finalizing the run.
                    anno: perform whole-genome annotation/analysis after
                        TE library construction.
    --overwrite [0|1]   If previous raw TE results are found, decide to overwrite
                (1, rerun) or not (0, default).
    --cds [File]    Provide a FASTA file containing the coding sequence (no introns,
            UTRs, nor TEs) of this genome or its close relative.
    --curatedlib [File] Provided a curated library to keep consistant naming and
                classification for known TEs. TEs in this file will be
                trusted 100%, so please ONLY provide MANUALLY CURATED ones.
                This option is not mandatory. It's totally OK if no file is
                provided (default).
    --sensitive [0|1]   Use RepeatModeler to identify remaining TEs (1) or not (0,
                default). This step is slow but MAY help to recover some TEs.
    --anno [0|1]    Perform (1) or not perform (0, default) whole-genome TE annotation
            after TE library construction.
    --rmout [File]  Provide your own homology-based TE annotation instead of using the
            EDTA library for masking. File is in RepeatMasker .out format. This
            file will be merged with the structural-based TE annotation. (--anno 1
            required). Default: use the EDTA library for annotation.
    --evaluate [0|1]    Evaluate (1) classification consistency of the TE annotation.
                (--anno 1 required). Default: 0. This step is slow and does
                not change the annotation result.
    --exclude [File]    Exclude bed format regions from TE annotation. Default: undef.
                (--anno 1 required).
    --force [0|1]   When no confident TE candidates are found: 0, interrupt and exit
            (default); 1, use rice TEs to continue.
    --u [float] Neutral mutation rate to calculate the age of intact LTR elements.
            Intact LTR age is found in this file: *EDTA_raw/LTR/*.pass.list.
            Default: 1.3e-8 (per bp per year, from rice).
    --repeatmodeler [path]  The directory containing RepeatModeler (default: read from ENV)
    --repeatmasker [path]   The directory containing RepeatMasker (default: read from ENV)
    --check_dependencies Check if dependencies are fullfiled and quit
    --threads|-t [int]  Number of theads to run this script (default: 4)
    --debug  [0|1]  Retain intermediate files (default: 0)
    --help|-h   Display this help info

```


### run EDTA

```
(EDTA) wuzhikun@fat01 17:03:03 O_O /home/wuzhikun/Project/Vigna/EDTA 
$ perl /home/wuzhikun/software/EDTA/EDTA.pl --step all  --genome  /home/wuzhikun/Project/Vigna/scaffold/final3/Vigna_unguiculata_assembly.fasta  --species others --cds /home/wuzhikun/Project/Vigna/GenePrediction/gene/Vigna_unguiculata_assembly.all.cds.fasta --overwrite 1 --sensitive 1 --anno 1 --evaluate 1 --threads 30

```


### pan-EDTA

```
(EDTA) wuzhikun@mu03 17:08:39 ^_^ /home/wuzhikun/software/EDTA 
$ sh /home/wuzhikun/software/EDTA/panEDTA.sh
ERROR: The genomes  file is not found or is empty

Usage: /home/wuzhikun/software/EDTA/panEDTA.sh -g genome_list.txt -c cds.fasta -t 10
    -g  A list of genome files with paths accessible from the working directory.
            Required: You can provide only a list of genomes in this file (one column, one genome each row).
            Optional: You can also provide both genomes and CDS files in this file (two columns, one genome and one CDS each row).
                  Missing of CDS files (eg, for some or all genomes) is allowed.
    -c      Optional. Coding sequence files in fasta format.
            The CDS file provided via this parameter will fill in the missing CDS files in the genome list.
            If no CDS files are provided in the genome list, then this CDS file will be used on all genomes.
    -l  Optional. A manually curated, non-redundant library following the RepeatMasker naming format.
    -f  Minimum number of full-length TE copies in individual genomes to be kept as candidate TEs for the pangenome.
            Lower is more inclusive, and will ↑ library size, ↑ sensitivity, and ↑ inconsistency.
            Higher is more stringent, and will ↓ library size, ↓ sensitivity, and ↓ inconsistency.
            Default: 3.
    -t  Number of CPUs to run panEDTA. Default: 10.

```


