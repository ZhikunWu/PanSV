

## [NanoCaller](https://github.com/WGLab/NanoCaller)

### [NanoCaller manual](https://github.com/WGLab/NanoCaller/blob/master/docs/ONT%20Case%20Study.md)


### install

v3.0.1

```
conda install -c bioconda nanocaller
```

### parameters


```
(HiC-Pro) wuzhikun@mu02 09:47:48 O_O /home/wuzhikun/github/LR_SNV_Caller/pipeline 
$ NanoCaller --help
usage: NanoCaller [-h] [-mode MODE] [-seq SEQUENCING] [-cpu CPU] [-mincov MINCOV] [-maxcov MAXCOV]
                  [-keep_bam] [-o OUTPUT] [-prefix PREFIX] [-sample SAMPLE] [-include_bed INCLUDE_BED]
                  [-exclude_bed EXCLUDE_BED] [-start START] [-end END] [-p PRESET] -bam BAM -ref REF -chrom
                  CHROM [-snp_model SNP_MODEL] [-min_allele_freq MIN_ALLELE_FREQ]
                  [-min_nbr_sites MIN_NBR_SITES] [-nbr_t NEIGHBOR_THRESHOLD] [-sup]
                  [-indel_model INDEL_MODEL] [-ins_t INS_THRESHOLD] [-del_t DEL_THRESHOLD]
                  [-win_size WIN_SIZE] [-small_win_size SMALL_WIN_SIZE] [-impute_indel_phase] [-phase_bam]
                  [-enable_whatshap]

optional arguments:
  -h, --help            show this help message and exit

Required Arguments:
  -bam BAM, --bam BAM   Bam file, should be phased if 'indel' mode is selected (default: None)
  -ref REF, --ref REF   Reference genome file with .fai index (default: None)
  -chrom CHROM, --chrom CHROM
                        Chromosome (default: None)

Preset:
  -p PRESET, --preset PRESET
                        Apply recommended preset values for SNP and Indel calling parameters, options are
                        'ont', 'ul_ont', 'ul_ont_extreme', 'ccs' and 'clr'. 'ont' works well for any type of
                        ONT sequencing datasets. However, use 'ul_ont' if you have several ultra-long ONT
                        reads up to 100kbp long, and 'ul_ont_extreme' if you have several ultra-long ONT
                        reads up to 300kbp long. For PacBio CCS (HiFi) and CLR reads, use 'ccs'and 'clr'
                        respectively. Presets are described in detail here:
                        github.com/WGLab/NanoCaller/blob/master/docs/Usage.md#preset-options. (default: None)

Configurations:
  -mode MODE, --mode MODE
                        NanoCaller mode to run, options are 'snps', 'snps_unphased', 'indels' and 'both'.
                        'snps_unphased' mode quits NanoCaller without using WhatsHap for phasing. (default:
                        both)
  -seq SEQUENCING, --sequencing SEQUENCING
                        Sequencing type, options are 'ont', 'ul_ont', 'ul_ont_extreme', and 'pacbio'. 'ont'
                        works well for any type of ONT sequencing datasets. However, use 'ul_ont' if you have
                        several ultra-long ONT reads up to 100kbp long, and 'ul_ont_extreme' if you have
                        several ultra-long ONT reads up to 300kbp long. For PacBio CCS (HiFi) and CLR reads,
                        use 'pacbio'. (default: ont)
  -cpu CPU, --cpu CPU   Number of CPUs to use (default: 1)
  -mincov MINCOV, --mincov MINCOV
                        Minimum coverage to call a variant (default: 8)
  -maxcov MAXCOV, --maxcov MAXCOV
                        Maximum coverage of reads to use. If sequencing depth at a candidate site exceeds
                        maxcov then reads are downsampled. (default: 160)

Variant Calling Regions:
  -include_bed INCLUDE_BED, --include_bed INCLUDE_BED
                        Only call variants inside the intervals specified in the bgzipped and tabix indexed
                        BED file. If any other flags are used to specify a region, intersect the region with
                        intervals in the BED file, e.g. if -chom chr1 -start 10000000 -end 20000000 flags are
                        set, call variants inside the intervals specified by the BED file that overlap with
                        chr1:10000000-20000000. Same goes for the case when whole genome variant calling flag
                        is set. (default: None)
  -exclude_bed EXCLUDE_BED, --exclude_bed EXCLUDE_BED
                        Path to bgzipped and tabix indexed BED file containing intervals to ignore for
                        variant calling. BED files of centromere and telomere regions for the following
                        genomes are included in NanoCaller: hg38, hg19, mm10 and mm39. To use these BED files
                        use one of the following options: ['hg38', 'hg19', 'mm10', 'mm39']. (default: None)
  -start START, --start START
                        start, default is 1 (default: None)
  -end END, --end END   end, default is the end of contig (default: None)

SNP Calling:
  -snp_model SNP_MODEL, --snp_model SNP_MODEL
                        NanoCaller SNP model to be used (default: ONT-HG002)
  -min_allele_freq MIN_ALLELE_FREQ, --min_allele_freq MIN_ALLELE_FREQ
                        minimum alternative allele frequency (default: 0.15)
  -min_nbr_sites MIN_NBR_SITES, --min_nbr_sites MIN_NBR_SITES
                        minimum number of nbr sites (default: 1)
  -nbr_t NEIGHBOR_THRESHOLD, --neighbor_threshold NEIGHBOR_THRESHOLD
                        SNP neighboring site thresholds with lower and upper bounds seperated by comma, for
                        Nanopore reads '0.4,0.6' is recommended, for PacBio CCS anc CLR reads '0.3,0.7' and
                        '0.3,0.6' are recommended respectively (default: 0.4,0.6)
  -sup, --supplementary
                        Use supplementary reads (default: False)

Indel Calling:
  -indel_model INDEL_MODEL, --indel_model INDEL_MODEL
                        NanoCaller indel model to be used (default: ONT-HG002)
  -ins_t INS_THRESHOLD, --ins_threshold INS_THRESHOLD
                        Insertion Threshold (default: 0.4)
  -del_t DEL_THRESHOLD, --del_threshold DEL_THRESHOLD
                        Deletion Threshold (default: 0.6)
  -win_size WIN_SIZE, --win_size WIN_SIZE
                        Size of the sliding window in which the number of indels is counted to determine
                        indel candidate site. Only indels longer than 2bp are counted in this window. Larger
                        window size can increase recall, but use a maximum of 50 only (default: 40)
  -small_win_size SMALL_WIN_SIZE, --small_win_size SMALL_WIN_SIZE
                        Size of the sliding window in which indel frequency is determined for small indels
                        (default: 4)
  -impute_indel_phase, --impute_indel_phase
                        Infer read phase by rudimentary allele clustering if the no or insufficient phasing
                        information is available, can be useful for datasets without SNPs or regions with
                        poor phasing quality (default: False)

Output Options:
  -keep_bam, --keep_bam
                        Keep phased bam files. (default: False)
  -o OUTPUT, --output OUTPUT
                        VCF output path, default is current working directory (default: None)
  -prefix PREFIX, --prefix PREFIX
                        VCF file prefix (default: variant_calls)
  -sample SAMPLE, --sample SAMPLE
                        VCF file sample name (default: SAMPLE)

Phasing:
  -phase_bam, --phase_bam
                        Phase bam files if snps mode is selected. This will phase bam file without indel
                        calling. (default: False)
  -enable_whatshap, --enable_whatshap
                        Allow WhatsHap to change SNP genotypes when phasing using --distrust-genotypes and
                        --include-homozygous flags (this is not the same as regenotyping), considerably
                        increasing the time needed for phasing. It has a negligible effect on SNP calling
                        accuracy for Nanopore reads, but may make a small improvement for PacBio reads. By
                        default WhatsHap will only phase SNP calls produced by NanoCaller, but not change
                        their genotypes. (default: False)

```


### run
```
$ NanoCaller -bam /home/wuzhikun/Project/HuaXi2021/mapping/minimap2/M692-0.bam  -p ont -o /home/wuzhikun/Project/HuaXi2021/mapping/minimap2   -ref /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa  -cpu 20 -chrom 22 -prefix M698-0_nanocaller



2023-03-29 10:47:20.431917: Variant calling output can be found here: /home/wuzhikun/Project/HuaXi2021/mapping/minimap2/M698-0_nanocaller.final.vcf.gz

2023-03-29 10:47:20.431999: Printing variant calling summary statistics
Location                     : /home/wuzhikun/Project/HuaXi2021/mapping/minimap2/M698-0_nanocaller.final.vcf.gz
Failed Filters               : 0
Passed Filters               : 318223
SNPs                         : 174210
MNPs                         : 0
Insertions                   : 12162
Deletions                    : 123202
Indels                       : 8649
Same as reference            : 0
Phased Genotypes             : 75.3% (239466/318223)
SNP Transitions/Transversions: 1.64 (136911/83350)
Total Het/Hom ratio          : 3.56 (248489/69734)
SNP Het/Hom ratio            : 2.82 (128555/45655)
MNP Het/Hom ratio            : - (0/0)
Insertion Het/Hom ratio      : 22.57 (11646/516)
Deletion Het/Hom ratio       : 4.37 (100259/22943)
Indel Het/Hom ratio          : 12.95 (8029/620)
Insertion/Deletion ratio     : 0.10 (12162/123202)
Indel/SNP+MNP ratio          : 0.83 (144013/174210)


2023-03-29 10:47:23.732759: Saving variant calling summary statistics to: /home/wuzhikun/Project/HuaXi2021/mapping/minimap2/M698-0_nanocaller.summary

2023-03-29 10:47:23.734940: Total Time Elapsed: 3036.63 seconds
l

```


out files:
```
-rw-rw-r-- 1 wuzhikun wuzhikun 2.5M Mar 29 10:43 M698-0_nanocaller.indels.vcf.gz
-rw-rw-r-- 1 wuzhikun wuzhikun  17K Mar 29 10:43 M698-0_nanocaller.indels.vcf.gz.tbi
-rw-rw-r-- 1 wuzhikun wuzhikun 4.2M Mar 29 10:43 M698-0_nanocaller.final.vcf.gz
-rw-rw-r-- 1 wuzhikun wuzhikun  22K Mar 29 10:43 M698-0_nanocaller.final.vcf.gz.tbi
-rw-rw-r-- 1 wuzhikun wuzhikun  894 Mar 29 10:43 M698-0_nanocaller.summary

```

