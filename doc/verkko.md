
## [verkko](https://github.com/marbl/verkko)

telomere-to-telomere assembly of PacBio HiFi and Oxford Nanopore reads

### install verkko

```
conda install -c bioconda verkko
```

### parameters
```
(Assembly3) wuzhikun@mu02 10:10:22 ^_^ /home/wuzhikun/Project/Vigna/clean 
$ verkko --help
usage: /home/wuzhikun/anaconda3/envs/Assembly3/bin/verkko -d <output-directory> -hifi <hifi-reads ...> -nano <nanopore-reads ...>
  MANDATORY PARAMETERS:
    -d <output-directory>    Directory to use for verkko intermediate and final results.
                             Will be created if needed.
    --hifi <files ...>       List of files containing PacBio HiFi reads.
    --nano <files ...>       List of files containing Oxford Nanopore reads.

                             Input reads can be any combination of FASTA/FASTQ,
                             uncompressed or gzip/bzip2/xz compressed.  Any
                             number of files can be supplied; *.gz works.

  ALGORITHM PARAMETERS:
    --no-correction


    --base-k
    --max-k
    --window
    --threads
    
    --split-bases
    --split-reads
    --min-ont-length
    
    --correct-k-mer-size
    --correct-mer-threshold
    --correct-min-read-length
    --correct-min-overlap-length
    
    --seed-min-length
    --seed-max-length
    --align-bandwidth
    --score-fraction
    --min-identity
    --min-score
    --end-clipping
    --incompatible-cutoff
    --max-trace

  COMPUTATIONAL PARAMETERS:
    --python <interpreter>   Path or name of a python interpreter.  Default: 'python'.
    --mbg <path>             Path to MBG.             Default for both is the
    --graphaligner <path>    Path to GraphAligner.    one packaged with verkko.

    --local                  Run on the local machine (default).
    --sge                    Enable Sun Grid Engine support.
    --slurm                  Enable Slurm support.
    --lsf                    Enable IBM Spectrum LSF support.

    --snakeopts <string>     Append snakemake options in "string" to the
                             snakemake command.  Options MUST be quoted.

    --sto-run                Set resource limits for various stages.
    --mer-run                Format: number-of-cpus memory-in-gb time-in-hours
    --ovb-run                  --cns-run 8 32 2
    --red-run
    --mbg-run
    --utg-run
    --spl-run
    --ali-run
    --pop-run
    --utp-run
    --lay-run
    --sub-run
    --par-run
    --cns-run

  Verkko module path: /home/wuzhikun/anaconda3/envs/Assembly3/lib/verkko/

No output directory (-d) set.
No PacBio HiFi reads (-hifi) supplied.
No Oxford Nanopore reads (-nano) supplied.
```


### run verkko

```
source activate verkko &&  verkko -d /home/wuzhikun/Project/Vigna/verkko --hifi /home/wuzhikun/data/Vigna/HiFi/m64181_210607_093522.ccs.fasta  --nano /home/wuzhikun/data/bean/Nanopore/fastq/Vigna_unguiculata_ONT.fastq.gz --min-ont-length 10000  --threads 20   > /home/wuzhikun/Project/Vigna/verkko/verkko.log 2>&1

```




output files:

```
drwxrwxr-x 1 4.0K Nov  8 00:44 0-correction
drwxrwxr-x 1 4.0K Nov  8 07:12 1-buildGraph
drwxrwxr-x 1 4.0K Nov  8 20:00 2-processGraph
drwxrwxr-x 1 4.0K Nov  8 19:37 3-align
drwxrwxr-x 1 4.0K Nov  8 20:00 4-processONT
drwxrwxr-x 1 4.0K Nov  8 20:00 5-untip
drwxrwxr-x 1 4.0K Nov  8 20:05 6-layoutContigs
drwxrwxr-x 1 4.0K Nov  9 08:15 7-consensus
-rw-rw-r-- 1 581M Nov  8 23:41 assembly.fasta
-rw-rw-r-- 1  55K Nov  8 23:41 assembly.hifi-coverage.csv
-rw-rw-r-- 1 387M Nov  8 23:41 assembly.homopolymer-compressed.gfa
-rw-rw-r-- 1 108M Nov  8 23:41 assembly.homopolymer-compressed.layout
-rw-rw-r-- 1 328K Nov  8 23:41 assembly.homopolymer-compressed.noseq.gfa
-rw-rw-r-- 1  60K Nov  8 23:41 assembly.ont-coverage.csv
-rw-rw-r-- 1    0 Nov  7 10:37 emptyfile
-rw-rw-r-- 1 5.4G Nov  8 00:44 hifi-corrected.fasta.gz
-rwxrwxr-x 1  385 Nov  7 10:37 snakemake.sh
-rw------- 1  159 Nov  8 23:45 verkko.e391195
-rw-rw-r-- 1 155K Nov  8 23:41 verkko.log
-rw------- 1    0 Nov  7 10:41 verkko.o391195
-rwxr--r-- 1  386 Nov  7 10:37 verkko.pbs
-rw-rw-r-- 1  28K Feb 28 22:01 verkko_Vu02_Vu06.log
-rw-rw-r-- 1 3.0K Nov  7 10:37 verkko.yml


```

```
$ assembly-stats  unitig-popped.unassigned.fasta
stats for unitig-popped.unassigned.fasta
sum = 608417212, n = 3041, ave = 200071.43, largest = 54801929
N50 = 28801511, n = 8
N60 = 27670232, n = 10
N70 = 13878996, n = 14
N80 = 11014841, n = 19
N90 = 52932, n = 131
N100 = 3393, n = 3041
N_count = 0
Gaps = 0


$ seqkit seq --min-len 1000000  unitig-popped.unassigned.fasta > unitig-popped.unassigned_1Mb.fasta
[INFO] flag -g (--remove-gaps) is switched on when using -m (--min-len) or -M (--max-len)


$ assembly-stats  unitig-popped.unassigned_1Mb.fasta
stats for unitig-popped.unassigned_1Mb.fasta
sum = 534980042, n = 29, ave = 18447587.66, largest = 54801929
N50 = 31706152, n = 7
N60 = 28726364, n = 9
N70 = 18975275, n = 11
N80 = 13878996, n = 14
N90 = 11466408, n = 18
N100 = 1311137, n = 29
N_count = 0
Gaps = 0


```


```

$ seqkit split -i unitig-popped.unassigned_1Mb.fasta

$ ls -rhSs
total 519M
1.3M unitig-popped.unassigned_1Mb.id_unassigned-0002128.fasta
1.4M unitig-popped.unassigned_1Mb.id_unassigned-0002687.fasta
1.6M unitig-popped.unassigned_1Mb.id_unassigned-0000677.fasta
2.0M unitig-popped.unassigned_1Mb.id_unassigned-0003043.fasta
3.2M unitig-popped.unassigned_1Mb.id_unassigned-0000154.fasta
4.5M unitig-popped.unassigned_1Mb.id_unassigned-0002877.fasta
4.7M unitig-popped.unassigned_1Mb.id_unassigned-0002756.fasta
5.8M unitig-popped.unassigned_1Mb.id_unassigned-0001168.fasta
7.4M unitig-popped.unassigned_1Mb.id_unassigned-0000982.fasta
9.3M unitig-popped.unassigned_1Mb.id_unassigned-0002818.fasta
 11M unitig-popped.unassigned_1Mb.id_unassigned-0002691.fasta
 12M unitig-popped.unassigned_1Mb.id_unassigned-0002515.fasta
 12M unitig-popped.unassigned_1Mb.id_unassigned-0000958.fasta
 12M unitig-popped.unassigned_1Mb.id_unassigned-0000675.fasta
 13M unitig-popped.unassigned_1Mb.id_unassigned-0001683.fasta
 14M unitig-popped.unassigned_1Mb.id_unassigned-0002378.fasta
 15M unitig-popped.unassigned_1Mb.id_unassigned-0002472.fasta
 18M unitig-popped.unassigned_1Mb.id_unassigned-0001480.fasta
 19M unitig-popped.unassigned_1Mb.id_unassigned-0001016.fasta
 27M unitig-popped.unassigned_1Mb.id_unassigned-0001452.fasta
 28M unitig-popped.unassigned_1Mb.id_unassigned-0002089.fasta
 28M unitig-popped.unassigned_1Mb.id_unassigned-0002980.fasta
 31M unitig-popped.unassigned_1Mb.id_unassigned-0002677.fasta
 34M unitig-popped.unassigned_1Mb.id_unassigned-0001470.fasta
 35M unitig-popped.unassigned_1Mb.id_unassigned-0001268.fasta
 36M unitig-popped.unassigned_1Mb.id_unassigned-0000202.fasta
 36M unitig-popped.unassigned_1Mb.id_unassigned-0000087.fasta
 53M unitig-popped.unassigned_1Mb.id_unassigned-0000013.fasta
 54M unitig-popped.unassigned_1Mb.id_unassigned-0002033.fasta

```
