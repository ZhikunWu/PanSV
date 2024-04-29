
## [TransDecoder预测转录本的开放阅读框](https://lxz9.com/2021/01/19/TransDecoder/)


```
(Anno) wuzhikun@cu09 09:32:54 ^_^ /home/wuzhikun/Project/PanSV/GenePrediction/RNA/Stringtie/CPG01 
$ /home/wuzhikun/anaconda3/envs/Anno/bin/gtf_genome_to_cdna_fasta.pl CPG01_RNA_stringtie.gtf  /home/wuzhikun/Project/PanSV/Scaffold/CPG01.genome.fasta > CPG01_RNA_stringtie.fasta

$ /home/wuzhikun/anaconda3/envs/Anno/bin/gtf_to_alignment_gff3.pl CPG01_RNA_stringtie.gtf  > CPG01_RNA_stringtie.gff3
```


```
$ /home/wuzhikun/anaconda3/envs/Anno/bin/TransDecoder.LongOrfs -t CPG01_RNA_stringtie.fasta --output_dir /home/wuzhikun/Project/PanSV/GenePrediction/RNA/Stringtie/CPG01
```

outdir:
```
drwxrwxr-x 1 wuzhikun wuzhikun 4.0K Apr 26 09:37 CPG01_RNA_stringtie.fasta.transdecoder_dir


$ lk CPG01_RNA_stringtie.fasta.transdecoder_dir
total 141M
-rw-rw-r-- 1 wuzhikun wuzhikun   68 Apr 26 09:37 base_freqs.dat
-rw-rw-r-- 1 wuzhikun wuzhikun  27M Apr 26 09:40 longest_orfs.pep
-rw-rw-r-- 1 wuzhikun wuzhikun  69M Apr 26 09:40 longest_orfs.cds
-rw-rw-r-- 1 wuzhikun wuzhikun  46M Apr 26 09:40 longest_orfs.gff3
drwxrwxr-x 1 wuzhikun wuzhikun 4.0K Apr 26 09:40 __checkpoints_longorfs

```


```
(Anno) wuzhikun@cu09 09:43:16 O_O /home/wuzhikun/Project/PanSV/GenePrediction/RNA/Stringtie/CPG01 
$ /home/wuzhikun/anaconda3/envs/Anno/bin/cdna_alignment_orf_to_genome_orf.pl  /home/wuzhikun/Project/PanSV/GenePrediction/RNA/Stringtie/CPG01/CPG01_RNA_stringtie.fasta.transdecoder_dir/longest_orfs.gff3  CPG01_RNA_stringtie.gff3 CPG01_RNA_stringtie.fasta > CPG01_RNA_stringtie.fasta.transdecoder.genome.gff3

```

