
## [VCF2Dis](https://github.com/hewm2008/VCF2Dis?tab=readme-ov-file)

### install

```
wget https://github.com/hewm2008/VCF2Dis/archive/v1.50.tar.gz
tar -zxf VCF2Dis-1.50.tar.gz
cd VCF2Dis-1.50
sh make.sh
```

### parameters

```
$ /home/wuzhikun/software/VCF2Dis-1.50/bin/VCF2Dis

    Usage: VCF2Dis -InPut  <in.vcf>  -OutPut  <p_dis.mat>

        -InPut     <str>     Input one or muti GATK VCF genotype File
        -OutPut    <str>     OutPut Sample p-Distance matrix

        -InList    <str>     Input GATK muti-chr VCF Path List
        -SubPop    <str>     SubGroup SampleList of VCFFile [ALLsample]
        -Rand      <float>   Probability (0-1] for each site to join Calculation [1]
        -KeepMF              Keep the Middle File diff & Use matrix

        -help                Show more help [hewm2008 v1.50]

```

### run
```
$ /home/wuzhikun/software/VCF2Dis-1.50/bin/VCF2Dis -InPut Sample_jointcall.twoallel.maf005.filtmiss.recode.vcf -OutPut sample_dist.mat
Total Sample Number to construct p-distance matrix is [ 215 ]
Start To Cal ...
Start To Create P_distance ...
P_distance is created done ...
```



### cotruct tree

```
PHYLIPNEW-3.69.650/bin/fneighbor  -datafile p_dis.matrix  -outfile tree.out1.txt -matrixtype s -treetype n -outtreefile tree.out2.tre
```


