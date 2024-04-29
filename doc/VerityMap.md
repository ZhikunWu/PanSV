

## [VerityMap](https://github.com/ablab/VerityMap)


### parameters

```

(PanSV) wuzhikun@mu03 18:53:37 O_O /home/wuzhikun/software/VerityMap/build 
$ veritymap --help
Usage: veritymap [OPTIONS] [ASSEMBLY_FNAMES]...

Options:
  --reads PATH                    File with ONT/PacBio reads
  -o PATH                         Output folder  [required]
  -t INTEGER                      Threads
  -d [hifi-haploid|hifi-haploid-complete|hifi-diploid|ont-haploid-complete]
                                  Sequencing platform, supported types are:
                                  "hifi" for PacBio HiFi reads and "ont" for
                                  ONT reads.Please note that "ont" mode is
                                  experimental and should be used with extra
                                  care
  -f, --no-reuse                  Do not reuse old files
  --careful                       Run mapper in a careful mode to better
                                  detect inconsistencies. Can be time- and
                                  memory-consuming. Not recommended to run on
                                  the whole genome.
  -l TEXT                         Comma separated list of assembly labels
  --help                          Show this message and exit.

```


```
veritymap --reads /home/wuzhikun/Project/BAssembly/clean/Boleracea.HiFi.1k.fastq.gz  /home/wuzhikun/Project/BAssembly/Assembly/hifiOnt100kHic/hifiOnt100kHic.hic.p_ctg_200k.fasta  -o /home/wuzhikun/Project/BAssembly/pipeline/VerityMap/Boleracea  -d hifi-haploid -t 20 > Boleracea.log 2>&1
```



