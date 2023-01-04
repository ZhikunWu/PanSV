## [RagTag](https://github.com/malonge/RagTag)

### install 


```
conda install -c bioconda ragtag
```


### parameters

```
(Assembly3) wuzhikun@mu03 10:46:45 ^_^ /home/wuzhikun/database/genome/Vigna_unguiculata/NJ 
$ ragtag.py

RagTag: Tools for fast and flexible genome assembly scaffolding and improvement.
Version: v2.0.1

usage: ragtag.py <command> [options]
    
    assembly improvement:
      correct         homology-based misassembly correction
      scaffold        homology-based assembly scaffolding
      patch           homology-based assembly patching
      merge           scaffold merging
      
    file utilities:
      agp2fa          build a FASTA file from an AGP file
      agpcheck        check for valid AGP file format
      asmstats        assembly statistics
      splitasm        split an assembly at gaps
      delta2paf       delta to PAF file conversion
      paf2delta       PAF to delta file conversion
      updategff       update gff intervals
      

    options:
      -c, --citation  
      -v, --version

```


```
$ ragtag.py scaffold --help
usage: ragtag.py scaffold <reference.fa> <query.fa>

Homology-based assembly scaffolding: Order and orient sequences in 'query.fa' by comparing them to sequences in 'reference.fa'>

positional arguments:
  <reference.fa>       reference fasta file (uncompressed or bgzipped)
  <query.fa>           query fasta file (uncompressed or bgzipped)

optional arguments:
  -h, --help           show this help message and exit

scaffolding options:
  -e <exclude.txt>     list of reference sequences to ignore [null]
  -j <skip.txt>        list of query sequences to leave unplaced [null]
  -J <hard-skip.txt>   list of query headers to leave unplaced and exclude from 'chr0' ('-C') [null]
  -f INT               minimum unique alignment length [1000]
  --remove-small       remove unique alignments shorter than '-f'
  -q INT               minimum mapq (NA for Nucmer alignments) [10]
  -d INT               maximum alignment merge distance [100000]
  -i FLOAT             minimum grouping confidence score [0.2]
  -a FLOAT             minimum location confidence score [0.0]
  -s FLOAT             minimum orientation confidence score [0.0]
  -C                   concatenate unplaced contigs and make 'chr0'
  -r                   infer gap sizes. if not, all gaps are 100 bp
  -g INT               minimum inferred gap size [100]
  -m INT               maximum inferred gap size [100000]

input/output options:
  -o PATH              output directory [./ragtag_output]
  -w                   overwrite intermediate files
  -u                   add suffix to unplaced sequence headers

mapping options:
  -t INT               number of minimap2/unimap threads [1]
  --aligner PATH       aligner executable ('nucmer', 'unimap' or 'minimap2') [minimap2]
  --mm2-params STR     space delimited minimap2 parameters ['-x asm5']
  --unimap-params STR  space delimited unimap parameters ['-x asm5']
  --nucmer-params STR  space delimted nucmer parameters ['--maxmatch -l 100 -c 500']

```


### run RagTag

```
(Assembly3) wuzhikun@fat01 15:58:03 ^_^ /home/wuzhikun/Project/Vigna/RagTag/A3-15 
$ ragtag.py scaffold /home/wuzhikun/Project/Vigna/scaffold/final3/Vigna_unguiculata_assembly.fasta  /home/wuzhikun/Project/Vigna/MultiGenome/A3-15/03.ctg_graph/nd.asm.fasta  -o /home/wuzhikun/Project/Vigna/RagTag/A3-15 -w -u -t 20 --aligner minimap2 


Wed Jan  4 15:59:33 2023 --- VERSION: RagTag v2.0.1
Wed Jan  4 15:59:33 2023 --- CMD: ragtag.py scaffold /home/wuzhikun/Project/Vigna/scaffold/final3/Vigna_unguiculata_assembly.fasta /home/wuzhikun/Project/Vigna/MultiGenome/A3-15/03.ctg_graph/nd.asm.fasta -o /home/wuzhikun/Project/Vigna/RagTag/A3-15 -w -u -t 20 --aligner minimap2
Wed Jan  4 15:59:33 2023 --- INFO: Mapping the query genome to the reference genome
Wed Jan  4 15:59:33 2023 --- INFO: Running: minimap2 -x asm5 -t 20 /home/wuzhikun/Project/Vigna/scaffold/final3/Vigna_unguiculata_assembly.fasta /home/wuzhikun/Project/Vigna/MultiGenome/A3-15/03.ctg_graph/nd.asm.fasta > /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.asm.paf 2> /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.asm.paf.log
Wed Jan  4 16:01:08 2023 --- INFO: Finished running : minimap2 -x asm5 -t 20 /home/wuzhikun/Project/Vigna/scaffold/final3/Vigna_unguiculata_assembly.fasta /home/wuzhikun/Project/Vigna/MultiGenome/A3-15/03.ctg_graph/nd.asm.fasta > /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.asm.paf 2> /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.asm.paf.log
Wed Jan  4 16:01:08 2023 --- INFO: Reading whole genome alignments
Wed Jan  4 16:01:08 2023 --- INFO: Filtering and merging alignments
Wed Jan  4 16:01:09 2023 --- INFO: Ordering and orienting query sequences
Wed Jan  4 16:01:09 2023 --- INFO: Writing scaffolds
Wed Jan  4 16:01:09 2023 --- INFO: Writing: /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.agp
Wed Jan  4 16:01:09 2023 --- INFO: Running: ragtag_agp2fa.py /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.agp /home/wuzhikun/Project/Vigna/MultiGenome/A3-15/03.ctg_graph/nd.asm.fasta > /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.fasta 2> /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.err
Wed Jan  4 16:01:15 2023 --- INFO: Finished running : ragtag_agp2fa.py /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.agp /home/wuzhikun/Project/Vigna/MultiGenome/A3-15/03.ctg_graph/nd.asm.fasta > /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.fasta 2> /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.err
Wed Jan  4 16:01:15 2023 --- INFO: Running: ragtag_stats.py /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.agp /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.confidence.txt > /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.stats 2> /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.err
Wed Jan  4 16:01:15 2023 --- INFO: Finished running : ragtag_stats.py /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.agp /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.confidence.txt > /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.stats 2> /home/wuzhikun/Project/Vigna/RagTag/A3-15/ragtag.scaffold.err
Wed Jan  4 16:01:15 2023 --- INFO: Goodbye

```

output files:
```
-rw-rw-r-- 1 wuzhikun wuzhikun    0 Jan  4 15:55 ragtag.scaffold.err
-rw-rw-r-- 1 wuzhikun wuzhikun 1.9M Jan  4 15:57 ragtag.scaffold.asm.paf
-rw-rw-r-- 1 wuzhikun wuzhikun  760 Jan  4 15:57 ragtag.scaffold.asm.paf.log
-rw-rw-r-- 1 wuzhikun wuzhikun  13K Jan  4 15:57 ragtag.scaffold.agp
-rw-rw-r-- 1 wuzhikun wuzhikun 3.9K Jan  4 15:57 ragtag.scaffold.confidence.txt
-rw-rw-r-- 1 wuzhikun wuzhikun 476M Jan  4 15:57 ragtag.scaffold.fasta
-rw-rw-r-- 1 wuzhikun wuzhikun  112 Jan  4 15:57 ragtag.scaffold.stats

```

```
$ cat ragtag.scaffold.stats
placed_sequences    placed_bp   unplaced_sequences  unplaced_bp gap_bp  gap_sequences
112 498620352   5   293695  10100   101

```


stats before and after RagTag
```
stats for nd.asm.fasta
sum = 498914047, n = 117, ave = 4264222.62, largest = 34984592
N50 = 15185130, n = 12
N60 = 11643077, n = 16
N70 = 9679230, n = 21
N80 = 6811185, n = 27
N90 = 4550420, n = 35
N100 = 33175, n = 117
N_count = 0
Gaps = 0

```

```
stats for ragtag.scaffold.fasta
sum = 498924147, n = 16, ave = 31182759.19, largest = 67088261
N50 = 44071752, n = 5
N60 = 43731116, n = 7
N70 = 42768203, n = 8
N80 = 41538818, n = 9
N90 = 39867062, n = 10
N100 = 33175, n = 16
N_count = 10100
Gaps = 101

```


