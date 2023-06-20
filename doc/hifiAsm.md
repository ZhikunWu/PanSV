
### Assembly

```
$ hifiasm -o /home/wuzhikun/PublicData/HiFi_information/HiFi3/PRJEB58125/Chionographis_japonica.asm -t 60 /home/wuzhikun/PublicData/HiFi_information/HiFi3/PRJEB58125/PRJEB58125_HiFi.fastq.gz 
```


out files:
```
-rw-rw-r-- 1 wuzhikun wuzhikun  19G Jun 18 05:41 Chionographis_japonica.asm.ec.bin
-rw-rw-r-- 1 wuzhikun wuzhikun 9.3G Jun 18 05:42 Chionographis_japonica.asm.ovlp.source.bin
-rw-rw-r-- 1 wuzhikun wuzhikun 4.5G Jun 18 05:44 Chionographis_japonica.asm.ovlp.reverse.bin
-rw-rw-r-- 1 wuzhikun wuzhikun 2.7G Jun 18 06:04 Chionographis_japonica.asm.bp.r_utg.gfa
-rw-rw-r-- 1 wuzhikun wuzhikun  40M Jun 18 06:04 Chionographis_japonica.asm.bp.r_utg.noseq.gfa
-rw-rw-r-- 1 wuzhikun wuzhikun 9.4M Jun 18 06:09 Chionographis_japonica.asm.bp.r_utg.lowQ.bed
-rw-rw-r-- 1 wuzhikun wuzhikun 2.6G Jun 18 06:10 Chionographis_japonica.asm.bp.p_utg.gfa
-rw-rw-r-- 1 wuzhikun wuzhikun  38M Jun 18 06:10 Chionographis_japonica.asm.bp.p_utg.noseq.gfa
-rw-rw-r-- 1 wuzhikun wuzhikun 8.5M Jun 18 06:16 Chionographis_japonica.asm.bp.p_utg.lowQ.bed
-rw-rw-r-- 1 wuzhikun wuzhikun 1.4G Jun 18 06:17 Chionographis_japonica.asm.bp.p_ctg.gfa
-rw-rw-r-- 1 wuzhikun wuzhikun  21M Jun 18 06:17 Chionographis_japonica.asm.bp.p_ctg.noseq.gfa
-rw-rw-r-- 1 wuzhikun wuzhikun 4.7M Jun 18 06:20 Chionographis_japonica.asm.bp.p_ctg.lowQ.bed
-rw-rw-r-- 1 wuzhikun wuzhikun 1.4G Jun 18 06:20 Chionographis_japonica.asm.bp.hap1.p_ctg.gfa
-rw-rw-r-- 1 wuzhikun wuzhikun  20M Jun 18 06:20 Chionographis_japonica.asm.bp.hap1.p_ctg.noseq.gfa
-rw-rw-r-- 1 wuzhikun wuzhikun 4.6M Jun 18 06:23 Chionographis_japonica.asm.bp.hap1.p_ctg.lowQ.bed
-rw-rw-r-- 1 wuzhikun wuzhikun 1.2G Jun 18 06:24 Chionographis_japonica.asm.bp.hap2.p_ctg.gfa
-rw-rw-r-- 1 wuzhikun wuzhikun  19M Jun 18 06:24 Chionographis_japonica.asm.bp.hap2.p_ctg.noseq.gfa
-rw-rw-r-- 1 wuzhikun wuzhikun 3.7M Jun 18 06:26 Chionographis_japonica.asm.bp.hap2.p_ctg.lowQ.bed

```


### gfa to fasta
```
$ gfatools gfa2fa Chionographis_japonica.asm.bp.hap1.p_ctg.gfa  >  Chionographis_japonica.asm.bp.hap1.p_ctg.fasta
```


### stats
```
$ assembly-stats  Chionographis_japonica.asm.bp.hap1.p_ctg.fasta
stats for Chionographis_japonica.asm.bp.hap1.p_ctg.fasta
sum = 1379121753, n = 5242, ave = 263090.76, largest = 10754215
N50 = 1169935, n = 224
N60 = 777064, n = 371
N70 = 495114, n = 593
N80 = 285509, n = 955
N90 = 67893, n = 1915
N100 = 8770, n = 5242
N_count = 
```



