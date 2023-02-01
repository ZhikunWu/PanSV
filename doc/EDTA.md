
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


########################################################
##### Extensive de-novo TE Annotator (EDTA) v2.1.0  ####
##### Shujun Ou (shujun.ou.1@gmail.com)             ####
########################################################



Wed Jan  4 17:06:50 CST 2023	Dependency checking:
				All passed!

	A CDS file /home/wuzhikun/Project/Vigna/GenePrediction/gene/Vigna_unguiculata_assembly.all.cds.fasta is provided via e sure this is the DNA sequence of coding regions only.

Wed Jan  4 17:07:18 CST 2023	Obtain raw TE libraries using various structure-based programs: 
Wed Jan  4 17:07:18 CST 2023	EDTA_raw: Check dependencies, prepare working directories.

Wed Jan  4 17:07:21 CST 2023	Start to find LTR candidates.

Wed Jan  4 17:07:21 CST 2023	Identify LTR retrotransposon candidates from scratch.

Wed Jan  4 19:20:17 CST 2023	Finish finding LTR candidates.

Wed Jan  4 19:20:17 CST 2023	Start to find TIR candidates.

Wed Jan  4 19:20:17 CST 2023	Identify TIR candidates from scratch.

Species: others
Wed Jan  4 22:03:13 CST 2023	Finish finding TIR candidates.

Wed Jan  4 22:03:13 CST 2023	Start to find Helitron candidates.

Wed Jan  4 22:03:13 CST 2023	Identify Helitron candidates from scratch.

Thu Jan  5 03:23:24 CST 2023	Finish finding Helitron candidates.

Thu Jan  5 03:23:24 CST 2023	Execution of EDTA_raw.pl is finished!

Thu Jan  5 03:23:25 CST 2023	Obtain raw TE libraries finished.
				All intact TEs found by EDTA: 
					Vigna_unguiculata_assembly.fasta.mod.EDTA.intact.fa
					Vigna_unguiculata_assembly.fasta.mod.EDTA.intact.gff3

Thu Jan  5 03:23:25 CST 2023	Perform EDTA advance filtering for raw TE candidates and generate the stage 1 library: 

Thu Jan  5 03:37:05 CST 2023	EDTA advance filtering finished.

Thu Jan  5 03:37:05 CST 2023	Perform EDTA final steps to generate a non-redundant comprehensive TE library:

				Use RepeatModeler to identify any remaining TEs that are missed by structure-based methods.

2023-01-06 11:38:35,850 -WARNING- Grid computing is not available because DRMAA not configured properly: Could not find drmaa specify its full path using the environment variable DRMAA_LIBRARY_PATH
2023-01-06 11:38:35,939 -INFO- VARS: {'sequence': 'Vigna_unguiculata_assembly.fasta.mod.RM.consensi.fa', 'hmm_database': 'rex'nucl', 'prefix': 'Vigna_unguiculata_assembly.fasta.mod.RM.consensi.fa.rexdb', 'force_write_hmmscan': False, 'processors': 30mp', 'min_coverage': 20, 'max_evalue': 0.001, 'disable_pass2': False, 'pass2_rule': '80-80-80', 'no_library': False, 'no_revecleanup': False, 'p2_identity': 80.0, 'p2_coverage': 80.0, 'p2_length': 80.0}
2023-01-06 11:38:35,939 -INFO- checking dependencies:
2023-01-06 11:38:36,168 -INFO- hmmer	3.3.1	OK
2023-01-06 11:38:36,269 -INFO- blastn	2.10.0+	OK
2023-01-06 11:38:36,271 -INFO- check database rexdb
2023-01-06 11:38:36,271 -INFO- db path: /home/wuzhikun/anaconda3/envs/EDTA/lib/python3.6/site-packages/TEsorter/database
2023-01-06 11:38:36,271 -INFO- db file: REXdb_protein_database_viridiplantae_v3.0_plus_metazoa_v3.hmm
2023-01-06 11:38:36,273 -INFO- REXdb_protein_database_viridiplantae_v3.0_plus_metazoa_v3.hmm	OK
2023-01-06 11:38:36,273 -INFO- Start classifying pipeline
2023-01-06 11:38:36,342 -INFO- total 1705 sequences
2023-01-06 11:38:36,342 -INFO- translating `Vigna_unguiculata_assembly.fasta.mod.RM.consensi.fa` in six frames
/home/wuzhikun/anaconda3/envs/EDTA/lib/python3.6/site-packages/Bio/Seq.py:2338: BiopythonWarning: Partial codon, len(sequenceof three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.
  BiopythonWarning,
2023-01-06 11:38:39,541 -INFO- HMM scanning against `/home/wuzhikun/anaconda3/envs/EDTA/lib/python3.6/site-packages/TEsorter/otein_database_viridiplantae_v3.0_plus_metazoa_v3.hmm`
2023-01-06 11:38:39,763 -INFO- Creating server instance (pp-1.6.4.4)
2023-01-06 11:38:39,763 -INFO- Running on Python 3.6.12 linux
2023-01-06 11:38:42,868 -INFO- pp local server started with 30 workers
2023-01-06 11:38:43,301 -INFO- Task 0 started
2023-01-06 11:38:43,303 -INFO- Task 1 started
2023-01-06 11:38:43,304 -INFO- Task 2 started
2023-01-06 11:38:43,305 -INFO- Task 3 started
2023-01-06 11:38:43,306 -INFO- Task 4 started
2023-01-06 11:38:43,306 -INFO- Task 5 started
2023-01-06 11:38:43,307 -INFO- Task 6 started
2023-01-06 11:38:43,308 -INFO- Task 7 started
2023-01-06 11:38:43,309 -INFO- Task 8 started
2023-01-06 11:38:43,309 -INFO- Task 9 started
2023-01-06 11:38:43,310 -INFO- Task 10 started
2023-01-06 11:38:43,311 -INFO- Task 11 started
2023-01-06 11:38:43,311 -INFO- Task 12 started
2023-01-06 11:38:43,312 -INFO- Task 13 started
2023-01-06 11:38:43,313 -INFO- Task 14 started
2023-01-06 11:38:43,313 -INFO- Task 15 started
2023-01-06 11:38:43,314 -INFO- Task 16 started
2023-01-06 11:38:43,315 -INFO- Task 17 started
2023-01-06 11:38:43,315 -INFO- Task 18 started
2023-01-06 11:38:43,316 -INFO- Task 19 started
2023-01-06 11:38:43,317 -INFO- Task 20 started
2023-01-06 11:38:43,318 -INFO- Task 21 started
2023-01-06 11:38:43,319 -INFO- Task 22 started
2023-01-06 11:38:43,320 -INFO- Task 23 started
2023-01-06 11:38:43,321 -INFO- Task 24 started
2023-01-06 11:38:43,321 -INFO- Task 25 started
2023-01-06 11:38:43,322 -INFO- Task 26 started
2023-01-06 11:38:43,323 -INFO- Task 27 started
2023-01-06 11:38:43,323 -INFO- Task 28 started
2023-01-06 11:38:43,324 -INFO- Task 29 started
2023-01-06 11:38:51,558 -INFO- generating gene anntations
2023-01-06 11:38:51,729 -INFO- 77 sequences classified by HMM
2023-01-06 11:38:51,729 -INFO- see protein domain sequences in `Vigna_unguiculata_assembly.fasta.mod.RM.consensi.fa.rexdb.domion gff3 file in `Vigna_unguiculata_assembly.fasta.mod.RM.consensi.fa.rexdb.dom.gff3`
2023-01-06 11:38:51,729 -INFO- classifying the unclassified sequences by searching against the classified ones
2023-01-06 11:38:51,814 -INFO- using the 80-80-80 rule
2023-01-06 11:38:51,814 -INFO- run CMD: `makeblastdb -in ./tmp/pass1_classified.fa -dbtype nucl`
2023-01-06 11:38:51,932 -INFO- run CMD: `blastn -query ./tmp/pass1_unclassified.fa -db ./tmp/pass1_classified.fa -out ./tmp/p.fa.blastout -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs -num_threads 30`
2023-01-06 11:38:52,981 -INFO- 1 sequences classified in pass 2
2023-01-06 11:38:52,982 -INFO- total 78 sequences classified.
2023-01-06 11:38:52,982 -INFO- see classified sequences in `Vigna_unguiculata_assembly.fasta.mod.RM.consensi.fa.rexdb.cls.tsv
2023-01-06 11:38:52,982 -INFO- writing library for RepeatMasker in `Vigna_unguiculata_assembly.fasta.mod.RM.consensi.fa.rexdb
2023-01-06 11:38:53,072 -INFO- writing classified protein domains in `Vigna_unguiculata_assembly.fasta.mod.RM.consensi.fa.rex
2023-01-06 11:38:53,078 -INFO- Summary of classifications:
Order           Superfamily      # of Sequences# of Clade Sequences    # of Clades# of full Domains
LTR             Copia                        39             38              5              0
LTR             Gypsy                        10             10              6              0
pararetrovirus  unknown                       3              0              0              0
LINE            unknown                       8              0              0              0
TIR             EnSpm_CACTA                   3              0              0              0
TIR             MuDR_Mutator                  5              0              0              0
TIR             hAT                           4              0              0              0
Helitron        unknown                       4              0              0              0
Maverick        unknown                       1              0              0              0
mixture         mixture                       1              0              0              0
2023-01-06 11:38:53,078 -INFO- Pipeline done.
2023-01-06 11:38:53,079 -INFO- cleaning the temporary directory ./tmp
Fri Jan  6 11:46:24 CST 2023	Clean up TE-related sequences in the CDS file with TEsorter:

2023-01-06 11:46:26,462 -WARNING- Grid computing is not available because DRMAA not configured properly: Could not find drmaa specify its full path using the environment variable DRMAA_LIBRARY_PATH
2023-01-06 11:46:26,473 -INFO- VARS: {'sequence': 'Vigna_unguiculata_assembly.all.cds.fasta.code', 'hmm_database': 'rexdb', ', 'prefix': 'Vigna_unguiculata_assembly.all.cds.fasta.code.rexdb', 'force_write_hmmscan': False, 'processors': 30, 'tmp_dir':verage': 20, 'max_evalue': 0.001, 'disable_pass2': False, 'pass2_rule': '80-80-80', 'no_library': False, 'no_reverse': False,lse, 'p2_identity': 80.0, 'p2_coverage': 80.0, 'p2_length': 80.0}
2023-01-06 11:46:26,473 -INFO- checking dependencies:
2023-01-06 11:46:26,493 -INFO- hmmer	3.3.1	OK
2023-01-06 11:46:26,580 -INFO- blastn	2.10.0+	OK
2023-01-06 11:46:26,581 -INFO- check database rexdb
2023-01-06 11:46:26,581 -INFO- db path: /home/wuzhikun/anaconda3/envs/EDTA/lib/python3.6/site-packages/TEsorter/database
2023-01-06 11:46:26,581 -INFO- db file: REXdb_protein_database_viridiplantae_v3.0_plus_metazoa_v3.hmm
2023-01-06 11:46:26,581 -INFO- REXdb_protein_database_viridiplantae_v3.0_plus_metazoa_v3.hmm	OK
2023-01-06 11:46:26,581 -INFO- Start classifying pipeline
2023-01-06 11:46:27,027 -INFO- total 30594 sequences
2023-01-06 11:46:27,028 -INFO- translating `Vigna_unguiculata_assembly.all.cds.fasta.code` in six frames
/home/wuzhikun/anaconda3/envs/EDTA/lib/python3.6/site-packages/Bio/Seq.py:2338: BiopythonWarning: Partial codon, len(sequenceof three. Explicitly trim the sequence or add trailing N before translation. This may become an error in future.
  BiopythonWarning,
2023-01-06 11:47:58,704 -INFO- HMM scanning against `/home/wuzhikun/anaconda3/envs/EDTA/lib/python3.6/site-packages/TEsorter/otein_database_viridiplantae_v3.0_plus_metazoa_v3.hmm`
2023-01-06 11:48:02,081 -INFO- Creating server instance (pp-1.6.4.4)
2023-01-06 11:48:02,081 -INFO- Running on Python 3.6.12 linux
2023-01-06 11:48:05,223 -INFO- pp local server started with 30 workers
2023-01-06 11:48:05,255 -INFO- Task 0 started
2023-01-06 11:48:05,256 -INFO- Task 1 started
2023-01-06 11:48:05,256 -INFO- Task 2 started
2023-01-06 11:48:05,257 -INFO- Task 3 started
2023-01-06 11:48:05,258 -INFO- Task 4 started
2023-01-06 11:48:05,259 -INFO- Task 5 started
2023-01-06 11:48:05,259 -INFO- Task 6 started
2023-01-06 11:48:05,260 -INFO- Task 7 started
2023-01-06 11:48:05,261 -INFO- Task 8 started
2023-01-06 11:48:05,262 -INFO- Task 9 started
2023-01-06 11:48:05,262 -INFO- Task 10 started
2023-01-06 11:48:05,263 -INFO- Task 11 started
2023-01-06 11:48:05,264 -INFO- Task 12 started
2023-01-06 11:48:05,264 -INFO- Task 13 started
2023-01-06 11:48:05,265 -INFO- Task 14 started
2023-01-06 11:48:05,266 -INFO- Task 15 started
2023-01-06 11:48:05,267 -INFO- Task 16 started
2023-01-06 11:48:05,267 -INFO- Task 17 started
2023-01-06 11:48:05,268 -INFO- Task 18 started
2023-01-06 11:48:05,269 -INFO- Task 19 started
2023-01-06 11:48:05,270 -INFO- Task 20 started
2023-01-06 11:48:05,271 -INFO- Task 21 started
2023-01-06 11:48:05,272 -INFO- Task 22 started
2023-01-06 11:48:05,272 -INFO- Task 23 started
2023-01-06 11:48:05,273 -INFO- Task 24 started
2023-01-06 11:48:05,274 -INFO- Task 25 started
2023-01-06 11:48:05,274 -INFO- Task 26 started
2023-01-06 11:48:05,275 -INFO- Task 27 started
2023-01-06 11:48:05,276 -INFO- Task 28 started
2023-01-06 11:48:05,276 -INFO- Task 29 started
2023-01-06 11:50:18,241 -INFO- generating gene anntations
2023-01-06 11:50:21,721 -INFO- 611 sequences classified by HMM
2023-01-06 11:50:21,722 -INFO- see protein domain sequences in `Vigna_unguiculata_assembly.all.cds.fasta.code.rexdb.dom.faa` f3 file in `Vigna_unguiculata_assembly.all.cds.fasta.code.rexdb.dom.gff3`
2023-01-06 11:50:21,722 -INFO- classifying the unclassified sequences by searching against the classified ones
2023-01-06 11:50:23,430 -INFO- using the 80-80-80 rule
2023-01-06 11:50:23,430 -INFO- run CMD: `makeblastdb -in ./tmp/pass1_classified.fa -dbtype nucl`
2023-01-06 11:50:23,590 -INFO- run CMD: `blastn -query ./tmp/pass1_unclassified.fa -db ./tmp/pass1_classified.fa -out ./tmp/p.fa.blastout -outfmt '6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen qcovs -num_threads 30`
2023-01-06 11:50:44,131 -INFO- 176 sequences classified in pass 2
2023-01-06 11:50:44,132 -INFO- total 787 sequences classified.
2023-01-06 11:50:44,132 -INFO- see classified sequences in `Vigna_unguiculata_assembly.all.cds.fasta.code.rexdb.cls.tsv`
2023-01-06 11:50:44,132 -INFO- writing library for RepeatMasker in `Vigna_unguiculata_assembly.all.cds.fasta.code.rexdb.cls.l
2023-01-06 11:50:45,511 -INFO- writing classified protein domains in `Vigna_unguiculata_assembly.all.cds.fasta.code.rexdb.cls
2023-01-06 11:50:45,548 -INFO- Summary of classifications:
Order           Superfamily      # of Sequences# of Clade Sequences    # of Clades# of full Domains
LTR             Bel-Pao                       9              0              0              0
LTR             Copia                       210            163             15              0
LTR             Gypsy                       237            190             14              0
LTR             Retrovirus                    2              0              0              0
pararetrovirus  unknown                       9              0              0              0
DIRS            unknown                       2              0              0              0
Penelope        unknown                       2              0              0              0
LINE            unknown                      27              0              0              0
TIR             MuDR_Mutator                 89              0              0              0
TIR             PIF_Harbinger                 7              0              0              0
TIR             hAT                          57              0              0              0
Helitron        unknown                      31              0              0              0
Maverick        unknown                     100              0              0              0
mixture         mixture                       5              0              0              0
2023-01-06 11:50:45,550 -INFO- Pipeline done.
2023-01-06 11:50:45,550 -INFO- cleaning the temporary directory ./tmp
				Remove CDS-related sequences in the EDTA library:

Fri Jan  6 12:01:42 CST 2023	EDTA final stage finished! You may check out:
				The final EDTA TE library: Vigna_unguiculata_assembly.fasta.mod.EDTA.TElib.fa
Fri Jan  6 12:01:42 CST 2023	Perform post-EDTA analysis for whole-genome annotation:

Fri Jan  6 12:01:42 CST 2023	Homology-based annotation of TEs using Vigna_unguiculata_assembly.fasta.mod.EDTA.TElib.fa fro

Fri Jan  6 13:15:54 CST 2023	TE annotation using the EDTA library has finished! Check out:
				Whole-genome TE annotation (total TE: 45.06%): Vigna_unguiculata_assembly.fasta.mod.EDTA.TEan
				Whole-genome TE annotation summary: Vigna_unguiculata_assembly.fasta.mod.EDTA.TEanno.sum
				Low-threshold TE masking for MAKER gene annotation (masked: 19.40%): Vigna_unguiculata_assembR.masked

Fri Jan  6 13:15:55 CST 2023	Evaluate the level of inconsistency for whole-genome TE annotation (slow step):

Mon Jan  9 21:03:15 CST 2023	Evaluation of TE annotation finished! Check out these files:

				Overall: Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.stat.all.sum
				Nested: Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.stat.nested.sum
				Non-nested: Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.stat.redun.sum

				If you want to learn more about the formatting and information of these files, please visit:
					https://github.com/oushujun/EDTA/wiki/Making-sense-of-EDTA-usage-and-outputs---Q&A



```

It needs five days


output files:

```
lrwxrwxrwx 1 wuzhikun wuzhikun   77 Jan  4 17:03 Vigna_unguiculata_assembly.fasta -> /home/wuzhikun/Project/Vigna/scaffold/final3/Vigna_unguiculata_assembly.fasta
-rw-rw-r-- 1 wuzhikun wuzhikun 506M Jan  4 17:03 Vigna_unguiculata_assembly.fasta.mod
lrwxrwxrwx 1 wuzhikun wuzhikun   89 Jan  4 17:03 Vigna_unguiculata_assembly.all.cds.fasta -> /home/wuzhikun/Project/Vigna/GenePrediction/gene/Vigna_unguiculata_assembly.all.cds.fasta
drwxrwxr-x 1 wuzhikun wuzhikun 4.0K Jan  5 03:19 Vigna_unguiculata_assembly.fasta.mod.EDTA.raw
drwxrwxr-x 1 wuzhikun wuzhikun 4.0K Jan  5 03:33 Vigna_unguiculata_assembly.fasta.mod.EDTA.combine
-rw-rw-r-- 1 wuzhikun wuzhikun 4.7M Jan  6 11:55 Vigna_unguiculata_assembly.fasta.mod.EDTA.TElib.fa
-rw-rw-r-- 1 wuzhikun wuzhikun 3.9M Jan  6 11:57 Vigna_unguiculata_assembly.fasta.mod.EDTA.intact.gff3
drwxrwxr-x 1 wuzhikun wuzhikun 4.0K Jan  6 11:57 Vigna_unguiculata_assembly.fasta.mod.EDTA.final
-rw-rw-r-- 1 wuzhikun wuzhikun 498M Jan  6 13:11 Vigna_unguiculata_assembly.fasta.mod.MAKER.masked
-rw-rw-r-- 1 wuzhikun wuzhikun  96M Jan  6 13:11 Vigna_unguiculata_assembly.fasta.mod.EDTA.TEanno.gff3
-rw-rw-r-- 1 wuzhikun wuzhikun 294K Jan  6 13:11 Vigna_unguiculata_assembly.fasta.mod.EDTA.TEanno.sum
drwxrwxr-x 1 wuzhikun wuzhikun 4.0K Jan  9 20:59 Vigna_unguiculata_assembly.fasta.mod.EDTA.anno

```

```
-rw-rw-r-- 1 wuzhikun wuzhikun  54M Jan  9 12:28 Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.iter31
-rw-rw-r-- 1 wuzhikun wuzhikun  54M Jan  9 14:34 Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.iter32
-rw-rw-r-- 1 wuzhikun wuzhikun  54M Jan  9 16:40 Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.iter33
-rw-rw-r-- 1 wuzhikun wuzhikun  64M Jan  9 17:11 Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.stat
-rw-rw-r-- 1 wuzhikun wuzhikun  54M Jan  9 18:49 Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.iter34
-rw-rw-r-- 1 wuzhikun wuzhikun  54M Jan  9 20:59 Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.cln
-rw-rw-r-- 1 wuzhikun wuzhikun 1.1K Jan  9 20:59 Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.stat.nested.sum
-rw-rw-r-- 1 wuzhikun wuzhikun 1.1K Jan  9 20:59 Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.stat.all.sum
-rw-rw-r-- 1 wuzhikun wuzhikun 1.1K Jan  9 20:59 Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.stat.redun.sum


$ grep -c '^>' Vigna_unguiculata_assembly.fasta.mod.EDTA.TE.fa.cln
145850


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


