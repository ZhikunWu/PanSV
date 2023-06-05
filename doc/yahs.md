
## [yahs](https://github.com/sanger-tol/yahs)

yet another Hi-C scaffolding tool


### install
```
conda install -c bioconda yahs
```


```
$ yahs --help
Usage: yahs [options] <contigs.fa> <hic.bed>|<hic.bam>|<hic.bin>
Options:
    -a FILE           AGP file (for rescaffolding) [none]
    -r INT[,INT,...]  list of resolutions in ascending order [automate]
    -e STR            restriction enzyme cutting sites [none]
    -l INT            minimum length of a contig to scaffold [0]
    -q INT            minimum mapping quality [10]
    -o STR            prefix of output files [yahs.out]
    -v INT            verbose level [0]
    --version         show version number

```

