## [Polypolish](https://github.com/rrwick/Polypolish)


### [polypolish manual](https://github.com/rrwick/Polypolish/wiki/How-to-run-Polypolish)


### install

```
conda install -c bioconda polypolish
```


### parameters
```
$ polypolish_insert_filter.py
usage: polypolish_insert_filter.py --in1 IN1 --in2 IN2 --out1 OUT1 --out2 OUT2 [--orientation {fr,rf,ff,rr,auto}]
                                   [--low LOW] [--high HIGH] [-h] [--version]

  _____        _                       _  _       _     
 |  __ \      | |                     | |(_)     | |    
 | |__) |___  | | _   _  _ __    ___  | | _  ___ | |__  
 |  ___// _ \ | || | | || '_ \  / _ \ | || |/ __|| '_ \ 
 | |   | (_) || || |_| || |_) || (_) || || |\__ \| | | |
 |_|    \___/ |_| \__, || .__/  \___/ |_||_||___/|_| |_|
                   __/ || |                             
                  |___/ |_|                             

Polypolish insert size alignment filter v0.5.0
github.com/rrwick/Polypolish

Inputs:
  --in1 IN1                         Input SAM file (first read in pairs)
  --in2 IN2                         Input SAM file (first second in pairs)

Outputs:
  --out1 OUT1                       Output SAM file (first read in pairs)
  --out2 OUT2                       Output SAM file (first second in pairs)

Settings:
  --orientation {fr,rf,ff,rr,auto}  Expected pair orientation (default: determine automatically)
  --low LOW                         Low percentile threshold (default: 0.1)
  --high HIGH                       High percentile threshold (default: 99.9)

Help:
  -h, --help                        Show this help message and exit
  --version                         Show program's version number and exit

```


### quick start

```
bwa index draft.fasta
bwa mem -t 16 -a draft.fasta reads_1.fastq.gz > alignments_1.sam
bwa mem -t 16 -a draft.fasta reads_2.fastq.gz > alignments_2.sam
polypolish_insert_filter.py --in1 alignments_1.sam --in2 alignments_2.sam --out1 filtered_1.sam --out2 filtered_2.sam
polypolish draft.fasta filtered_1.sam filtered_2.sam > polished.fasta

```


```
(Assembly3) wuzhikun@mu02 09:34:28 O_O /home/wuzhikun/github/UsefulWeb/doc 
$ polypolish --help
  _____        _                       _  _       _
 |  __ \      | |                     | |(_)     | |
 | |__) |___  | | _   _  _ __    ___  | | _  ___ | |__
 |  ___// _ \ | || | | || '_ \  / _ \ | || |/ __|| '_ \
 | |   | (_) || || |_| || |_) || (_) || || |\__ \| | | |
 |_|    \___/ |_| \__, || .__/  \___/ |_||_||___/|_| |_|
                   __/ || |
                  |___/ |_|

Polypolish v0.5.0

short-read polishing of long-read assemblies
github.com/rrwick/Polypolish

USAGE:
    polypolish [OPTIONS] <ASSEMBLY> <SAM>...

ARGS:
    <ASSEMBLY>    Assembly to polish (one file in FASTA format)
    <SAM>...      Short read alignments (one or more files in SAM format)

OPTIONS:
    -d, --min_depth <MIN_DEPTH>                  A base must occur at least this many times in the pileup to be considered
                                                 valid [default: 5]
        --debug <DEBUG>                          Optional file to store per-base information for debugging purposes
    -h, --help                                   Print help information
    -i, --fraction_invalid <FRACTION_INVALID>    A base must make up less than this fraction of the read depth to be
                                                 considered invalid [default: 0.2]
    -m, --max_errors <MAX_ERRORS>                Ignore alignments with more than this many mismatches and indels [default:
                                                 10]
    -v, --fraction_valid <FRACTION_VALID>        A base must make up at least this fraction of the read depth to be
                                                 considered valid [default: 0.5]
    -V, --version                                Print version information

```











