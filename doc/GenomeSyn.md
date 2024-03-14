## [GenomeSyn](https://github.com/JM-SONG/GenomeSyn/)

```

seqkit seq -m 2000000 Vigna_unguiculata_assembly.fasta > Vigna_unguiculata_assembly-1.fasta


seqkit seq -m 5000000 vigun.IT97K-499-35.gnm1.QnBw.genome_main.fna  > vigun.IT97K-499-35.gnm1.QnBw.genome_main-1.fna


nucmer -c 400 -l 150 -t 20  -p IT97K2Fengchan Vigna_unguiculata_assembly-1.fasta vigun.IT97K-499-35.gnm1.QnBw.genome_main-1.fna


show-coords -T -r -c -l -k IT97K2Fengchan.delta  > IT97K2Fengchan.coords


delta-filter -i 90 -l 1000 -1  IT97K2Fengchan.delta > IT97K2Fengchan_filt.delta
show-coords -T -r -c -l -k IT97K2Fengchan_filt.delta   > IT97K2Fengchan_filt.coords

sed '1,5d' IT97K2Fengchan_filt.coords > IT97K2Fengchan_filt_deheader.coords



nucmer -c 400 -l 150 -t 20  -p A1472Fengchan  ~/Project/Vigna/Final/Vigna_unguiculata_assembly.fasta /home/wuzhikun/database/genome/Vigna_unguiculata/A147/GCA_026401125.1_ASM2640112v1_genomic_5M.fna





nucmer -c 400 -l 150 -t 20  -p A1472Fengchan  ~/Project/Vigna/Final/Vigna_unguiculata_assembly.fasta /home/wuzhikun/database/genome/Vigna_unguiculata/A147/GCA_026401125.1_ASM2640112v1_genomic_5M.fna





/home/wuzhikun/software/GenomeSyn-main/GenomeSyn-1.2.7/bin/GenomeSyn -g1 ~/Project/Vigna/Final/Vigna_unguiculata_assembly.fasta -g2 /home/wuzhikun/database/genome/Vigna_unguiculata/A147/GCA_026401125.1_ASM2640112v1_genomic_5M.fna -cf1  /home/wuzhikun/database/genome/Vigna_unguiculata/A147/A1472Fengchan_filt_deheader.coords -ref Fengchan6 -q1 A147 -cen1 /home/wuzhikun/Project/Vigna/pipeline/GenomeSyn/Vu_fengchan_centromere-1.txt  -tel1  /home/wuzhikun/Project/Vigna/pipeline/GenomeSyn/Vu_fengchan_telomere-1.txt

```




```
nucmer -c 400 -l 150 -t 20 -p Bo_align /home/wuzhikun/PublicData/Brassica_oleracea/GCA_034621495.1_ASM3462149v1_genomic.fna.split/GCA_034621495.1_chromosome.fasta /home/wuzhikun/Project/BAssembly/Assembly/NextDenovo50k/03.ctg_graph/ragtag/ragtag.scaffold.fasta.split/NextDenovo50k_chr.fasta

show-coords -r -c -l -T  Bo_align.delta  > Bo_align.coords

export PERL5LIB=/home/wuzhikun/anaconda3/envs/TE3/lib

(TE3) wuzhikun@mu02 10:17:17 ^_^ /home/wuzhikun/Project/BAssembly/pipeline/GenomeSyn 
$ /home/wuzhikun/software/GenomeSyn-main/GenomeSyn-1.2.7/bin/GenomeSyn  -t 3 -g1 GCA_034621495.1_chromosome.fasta -g2 NextDenovo50k_chr.fasta -cf1 Bo_align.coords

```

