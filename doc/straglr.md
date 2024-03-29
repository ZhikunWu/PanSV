## [straglr](https://github.com/bcgsc/straglr)

Short-tandem repeat genotyping using long reads

### prameters
```
(HiC-Pro) wuzhikun@mu02 08:24:35 O_O /home/wuzhikun/software/straglr 
$ python /home/wuzhikun/software/straglr/straglr.py --help
usage: straglr.py [-h] [--reads_fasta READS_FASTA [READS_FASTA ...]] [--min_ins_size MIN_INS_SIZE]
                  [--exclude EXCLUDE] [--regions REGIONS] [--nprocs NPROCS] [--chroms CHROMS [CHROMS ...]]
                  [--loci LOCI] [--include_alt_chroms] [--min_support MIN_SUPPORT]
                  [--min_cluster_size MIN_CLUSTER_SIZE] [--genotype_in_size] [--max_str_len MAX_STR_LEN]
                  [--min_str_len MIN_STR_LEN] [--max_num_clusters MAX_NUM_CLUSTERS] [--max_cov MAX_COV]
                  [--use_unpaired_clips] [--include_partials]
                  [--trf_args Match Mismatch Delta PM PI Minscore MaxPeriod] [--tmpdir TMPDIR] [--debug]
                  [--version]
                  bam genome_fasta out_prefix

positional arguments:
  bam                   bam file
  genome_fasta          genome_fasta
  out_prefix            output prefix

optional arguments:
  -h, --help            show this help message and exit
  --reads_fasta READS_FASTA [READS_FASTA ...]
                        read indexed fasta file(s)
  --min_ins_size MIN_INS_SIZE
                        minimum insertion size. Default:100
  --exclude EXCLUDE     bed file to exclude regions
  --regions REGIONS     bed file for scanning only specific regions
  --nprocs NPROCS       number of processes
  --chroms CHROMS [CHROMS ...]
                        chromosomes
  --loci LOCI           bed file of loci for genotyping
  --include_alt_chroms  include alternate chromosomes. By default, only chroms 1-22,X,Y are considered in
                        genome scan
  --min_support MIN_SUPPORT
                        minimum number of supporting reads for detecting expansion. Default:2
  --min_cluster_size MIN_CLUSTER_SIZE
                        minimum number of supporting reads for allele clustering. Default:2
  --genotype_in_size    report genotype in size instead of copy numbers
  --max_str_len MAX_STR_LEN
                        maximum STR length. Default:50
  --min_str_len MIN_STR_LEN
                        minimum STR length. Default:2
  --max_num_clusters MAX_NUM_CLUSTERS
                        maximum number of clusters to try. Default:2
  --max_cov MAX_COV     maximum allowed coverage for ins inspection. Default:100
  --use_unpaired_clips  also examine unpaired clipped alignments in genome scan
  --include_partials    detect and report reads only capturing partial repeats when genotyping
  --trf_args Match Mismatch Delta PM PI Minscore MaxPeriod
                        tandem repeat finder arguments. Default:2 5 5 80 10 10 500
  --tmpdir TMPDIR       directory to use for generating tmp files instead of system TEMP
  --debug               debug mode i.e. keep trf output
  --version             show program's version number and exit


```

