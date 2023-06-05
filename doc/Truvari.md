## [Truvari](https://github.com/ACEnglish/truvari)

```
(nanovar2) wuzhikun@mu02 19:30:23 ^_^ /home/wuzhikun/github/SCZSV/script 
$ truvari --help
usage: truvari [-h] CMD ...

Truvari v4.0.0 Structural Variant Benchmarking and Annotation

Available commands:
    bench         Performance metrics from comparison of two VCFs
    collapse      Collapse possibly redundant VCF entries
    anno          Annotate a VCF
    consistency   Consistency report between multiple VCFs
    vcf2df        Turn a VCF into a pandas DataFrame
    segment       Normalization of SVs into disjointed genomic regions
    stratify      Count variants per-region in vcf
    divide        Divide a VCF into independent shards
    phab          Variant harmonization using MSA
    refine        Automated bench result refinement with phab
    version       Print the Truvari version and exit

positional arguments:
  CMD         Command to execute
  OPTIONS     Options to pass to the command

options:
  -h, --help  show this help message and exit
```

```
(nanovar2) wuzhikun@mu02 19:55:04 ^_^ /home/wuzhikun/github/SCZSV/script 
$ truvari collapse --help
usage: collapse [-h] -i INPUT [-o OUTPUT] [-c COLLAPSED_OUTPUT] [-f REFERENCE] [-k {first,maxqual,common}]
                [--debug] [-r REFDIST] [-p PCTSEQ] [-B MINHAPLEN] [-P PCTSIZE] [-O PCTOVL] [-t] [--use-lev] [--hap]
                [--chain] [--no-consolidate] [--null-consolidate NULL_CONSOLIDATE] [-s SIZEMIN] [-S SIZEMAX]
                [--passonly]

Structural variant collapser

Will collapse all variants within sizemin/max that match over thresholds

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        Input variants
  -o OUTPUT, --output OUTPUT
                        Output vcf (stdout)
  -c COLLAPSED_OUTPUT, --collapsed-output COLLAPSED_OUTPUT
                        Where collapsed variants are written (collapsed.vcf)
  -f REFERENCE, --reference REFERENCE
                        Indexed fasta used to call variants
  -k {first,maxqual,common}, --keep {first,maxqual,common}
                        When collapsing calls, which one to keep (first)
  --debug               Verbose logging
  --hap                 Collapsing a single individual's haplotype resolved calls (False)
  --chain               Chain comparisons to extend possible collapsing (False)
  --no-consolidate      Skip consolidation of sample genotype fields (True)
  --null-consolidate NULL_CONSOLIDATE
                        Comma separated list of FORMAT fields to consolidate into the kept entry by taking the
                        first non-null from all neighbors (None)

Comparison Threshold Arguments:
  -r REFDIST, --refdist REFDIST
                        Max reference location distance (500)
  -p PCTSEQ, --pctseq PCTSEQ
                        Min percent sequence similarity. Set to 0 to ignore. (0.95)
  -B MINHAPLEN, --minhaplen MINHAPLEN
                        Minimum haplotype sequence length to create (50)
  -P PCTSIZE, --pctsize PCTSIZE
                        Min pct allele size similarity (minvarsize/maxvarsize) (0.95)
  -O PCTOVL, --pctovl PCTOVL
                        Min pct reciprocal overlap (0.0) for DEL events
  -t, --typeignore      Variant types don't need to match to compare (False)
  --use-lev             Use the Levenshtein distance ratio instead of edlib editDistance ratio (False)

Filtering Arguments:
  -s SIZEMIN, --sizemin SIZEMIN
                        Minimum variant size to consider for comparison (50)
  -S SIZEMAX, --sizemax SIZEMAX
                        Maximum variant size to consider for comparison (50000)
  --passonly            Only consider calls with FILTER == PASS

```

