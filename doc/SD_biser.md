
## [biser](https://github.com/0xTCG/biser/)



```
$ biser --help
usage: biser [-h] [--temp TEMP] [--threads THREADS] --output OUTPUT [--hard] [--keep-contigs] [--keep-temp]
             [--resume RESUME] [--no-decomposition] [--max-error MAX_ERROR] [--max-edit-error MAX_EDIT_ERROR]
             [--max-chromosome-size MAX_CHROMOSOME_SIZE] [--kmer-size KMER_SIZE] [--winnow-size WINNOW_SIZE]
             [--version] [--ld-path LD_PATH] [--gc-heap GC_HEAP]
             genomes [genomes ...]

Segmental duplication detection tool

positional arguments:
  genomes               Indexed genomes in FASTA format.

optional arguments:
  -h, --help            show this help message and exit
  --temp TEMP, -T TEMP  Temporary directory location
  --threads THREADS, -t THREADS
                        Number of threads
  --output OUTPUT, -o OUTPUT
                        Indexed genomes in FASTA format.
  --hard, -H            Are input genomes already hard-masked?
  --keep-contigs        Do not ignore contigs, unplaced sequences, alternate alleles, patch chromosomes and
                        mitochondrion sequences (i.e., chrM and chromosomes whose name contains underscore). Enable
                        this when running BISER on scaffolds and custom assemblies.
  --keep-temp, -k       Keep temporary directory after the execution. Useful for debugging.
  --resume RESUME       Resume the previously interrupted run (that was run with --keep-temp; needs the temp
                        directory for resume).
  --no-decomposition    Skip SD decomposition step.
  --max-error MAX_ERROR
                        Maximum SD error.
  --max-edit-error MAX_EDIT_ERROR
                        Maximum SD error (without large gaps).
  --max-chromosome-size MAX_CHROMOSOME_SIZE
                        Maximum chromosome size.
  --kmer-size KMER_SIZE
                        Search k-mer size.
  --winnow-size WINNOW_SIZE
                        Search winnow size.
  --version, -v         show program s version number and exit
  --ld-path LD_PATH     Override LD_LIBRARY_PATH (debug use only).
  --gc-heap GC_HEAP     Set GC_INITIAL_HEAP_SIZE.

```


### run

```
(Assembly3) wuzhikun@fat01 11:01:07 ^_^ /home/wuzhikun/Project/Centromere 
$ biser --threads 30 --output /home/wuzhikun/Project/Centromere/SD /home/wuzhikun/Project/Centromere/repeatMasker/Chromosome/Zygaena_filipendulae/Zygaena_filipendulae.fasta.masked


0. Hard-masking genomes
100%|████████████████████████████████████████| 1/1
  Hard-masking: 00:23s (single: 00:20s)

1. Putative SD detection
 98%|███████████████████████████████████████▎| 55/56
ERROR: BISER failed (-11):
[p] /home/wuzhikun/anaconda3/envs/Assembly3/lib/python3.8/site-packages/biser/exe/biser.exe search /tmp/biser.fwdxbvn4/genomes/Zygaena_filipendulae.fasta.fa -c OU015680.1 -o /tmp/biser.fwdxbvn4/search/Zygaena_filipendulae.fasta_OU015680.1_0.bed 0 --kmer-size=14 --winnow-size=16 --max-error=0.3 --max-edit-error=0.15 --max-chromosome-size=300000000
[o] 
[e] 
```


```
problem solved, just added the parameter (--gc-heap 10G)
```



```
(Assembly3) wuzhikun@fat01 11:48:38 O_O /home/wuzhikun/Project/Centromere 
$ biser --threads 30 --gc-heap 10G  --output /home/wuzhikun/Project/Centromere/SD /home/wuzhikun/Project/Centromere/repeatMasker/Chromosome/Zygaena_filipendulae/Zygaena_filipendulae.fasta.masked
Running BISER v1.4 on 1 genome(s): Zygaena_filipendulae.fasta

0. Hard-masking genomes
100%|████████████████████████████████████████| 1/1
  Hard-masking: 00:24s (single: 00:21s)

1. Putative SD detection
100%|████████████████████████████████████████| 56/56
  Search: 18:09s (single: 334:37s)

2. Putative SD alignment
Total alignments: 80,353
100%|████████████████████████████████████████| 60/60
  Align: 00:13s (single: 06:35s)

3. SD decomposition
100%|████████████████████████████████████████| 972/972
  Decomposition: 00:04s (single: 00:31s)
Done! Results are available in /home/wuzhikun/Project/Centromere/SD
BISER: 18:56s

```

output file:
```
-rw-rw-r-- 1 wuzhikun wuzhikun 2.7M Jun 13 12:04 /home/wuzhikun/Project/Centromere/SD
-rw-rw-r-- 1 wuzhikun wuzhikun 2.2M Jun 13 12:04 /home/wuzhikun/Project/Centromere/SD.elem.txt

```

--keep-contigs


