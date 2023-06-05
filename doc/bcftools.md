

### extract chromosome should index the vcf file

```
tabix -p vcf  variants_freeze4_indel_insdel_alt.vcf.gz
tabix -p vcf  variants_freeze4_snv_snv_alt.vcf.gz

(Assembly3) wuzhikun@mu03 10:44:51 ^_^ /home/wuzhikun/database/HGSVC2/freeze4 
$ bcftools view -s HG02818 --output-type v  -o variants_freeze4_indel_insdel_alt_HG0288.vcf --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX  variants_freeze4_indel_insdel_alt.vcf.gz


$ grep -cv '^#' variants_freeze4_indel_insdel_alt_HG0288.vcf
2614038
```



```
bcftools view -s HG02818 --output-type v  -o /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_snv.vcf --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX /home/wuzhikun/database/HGSVC2/freeze4/variants_freeze4_snv_snv_alt.vcf.gz >/home/wuzhikun/Project/Centromere/log/extractVariant_HG02818.log 2>&1
bcftools view --no-header -s HG02818 --output-type v  -o /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_indel.vcf --regions chr1,chr2,chr3,chr4,chr5,chr6,chr7,chr8,chr9,chr10,chr11,chr12,chr13,chr14,chr15,chr16,chr17,chr18,chr19,chr20,chr21,chr22,chrX /home/wuzhikun/database/HGSVC2/freeze4/variants_freeze4_indel_insdel_alt.vcf.gz >/home/wuzhikun/Project/Centromere/log/extractVariant_HG02818.log 2>&1
cat /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_snv.vcf /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_indel.vcf | grep -v '\.|\.' | grep -v '0|0'  > /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_snv_indel.vcf


$ grep -cv '^#' /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_snv_indel.vcf
5658153


```

```
$ cat /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_snv.vcf /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_indel.vcf | grep -v '\.|\.' | grep -v '0|0'  | bcftools sort --output-file /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_snv_indel.vcf.gz  --output-type z -
Writing to /tmp/bcftools-sort.vIPtO3
Merging 3 temporary files
Cleaning
Done

```





```
bcftools consensus --fasta-ref /home/wuzhikun/database/HGSVC2/genome/hg38.no_alt_chrom.fa --haplotype 1 /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_snv_indel.vcf.gz | sed '/^>/ s/$/_hap1/' | bgzip -c > /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_hap1.fasta.gz
```


```
The fasta sequence does not match the REF allele at chr1:1423775:
   .vcf: [CTT]
   .vcf: [C] <- (ALT)
   .fa:  [Ttt]ttttttttttttttttttttttggct

```


```
chr1    1423775 chr1-1423775-SNV-C-T    C       T       .       .       ID=chr1-1423775-SNV-C-T;VARTYPE=SNV;SVTYPE=SNV;SVLEN=1;SAMPLE=HG02818;TIG_REGION=cluster5_000002F:4314-4314;QUERY_STRAND=-;AC=2;AN=2    GT      1|1
chr1    1423775 chr1-1423776-DEL-2      CTT     C       .       .       ID=chr1-1423776-DEL-2;VARTYPE=INDEL;SVTYPE=DEL;SVLEN=-2;SAMPLE=NA19238;TIG_REGION=cluster12_000030F:557791-557791;QUERY_STRAND=+;HOMOP_LEN=24;AC=2;AN=2 GT      1|1
```



```
$ bcftools norm --rm-dup all  --output /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_snv_indel_norm.vcf.gz --output-type z /home/wuzhikun/Project/Centromere/PairtoolsPhase/HG02818/genome/HG02818_genome_snv_indel.vcf.gz

```

