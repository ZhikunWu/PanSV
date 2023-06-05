

## [nanoq](https://github.com/esteinig/nanoq)

### install nanoq

```
conda install -c bioconda nanoq
```


### stats

```
$ nanoq -i CPG01_ONT_0.1.fastq.gz -s > CPG01_ONT_0.1.fastq.report.txt
```

output file

```
$ cat CPG01_ONT_0.1.fastq.report.txt
192971 4241548555 31022 150433 139 21980 19756 10.2 10.2

```


```
number of reads
number of base pairs
N50 read length
longest read
shorted reads
mean read length
median read length
mean read quality
median read quality
```

### filt 

```
(nanovar2) wuzhikun@fat02 13:47:53 ^_^ /home/wuzhikun/data/Vigna/Nanopore/CPG01 
$  nanoq -i CPG01_ONT_0.1.fastq.gz -S 70 -E 20  -O g -c 6 > CPG01_ONT_0.1_filt.fastq.gz
```


