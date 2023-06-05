
## [vg](https://github.com/vgteam/vg)

```
vg autoindex
vg rna
vg mpmap
rpvg

```

### parameters

```
$ vg autoindex --help
usage: vg autoindex [options]
options:
  output:
    -p, --prefix PREFIX    prefix to use for all output (default: index)
    -w, --workflow NAME    workflow to produce indexes for, can be provided multiple
                           times. options: map, mpmap, giraffe (default: map)
  input data:
    -r, --ref-fasta FILE   FASTA file containing the reference sequence (may repeat)
    -v, --vcf FILE         VCF file with sequence names matching -r (may repeat)
    -i, --ins-fasta FILE   FASTA file with sequences of INS variants from -v
    -g, --gfa FILE         GFA file to make a graph from
    -x, --tx-gff FILE      GTF/GFF file with transcript annotations (may repeat)
  configuration:
    -f, --gff-feature STR  GTF/GFF feature type (col. 3) to add to graph (default: exon)
    -a, --gff-tx-tag STR   GTF/GFF tag (in col. 9) for transcript ID (default: transcript_id)
  logging and computation:
    -T, --tmp-dir DIR      temporary directory to use for intermediate files
    -M, --target-mem MEM   target max memory usage (not exact, formatted INT[kMG])
                           (default: 1/2 of available)
    -t, --threads NUM      number of threads (default: all available)
    -V, --verbosity NUM    log to stderr (0 = none, 1 = basic, 2 = debug; default 1)
    -h, --help             print this help message to stderr and exit

```


```
vg autoindex --workflow mpmap -t 4 --prefix vg_rna --ref-fasta small/x.fa --vcf small/x.vcf.gz --tx-gff small/x.gtf

```

output file:
```
-rw-rw-r-- 1 wuzhikun wuzhikun 2.3K Nov 22 16:45 vg_rna.haplotx.gbwt
-rw-rw-r-- 1 wuzhikun wuzhikun 8.7K Nov 22 16:45 vg_rna.spliced.xg
-rw-rw-r-- 1 wuzhikun wuzhikun  10K Nov 22 16:45 vg_rna.spliced.dist
-rw-rw-r-- 1 wuzhikun wuzhikun 7.9K Nov 22 16:45 vg_rna.spliced.gcsa
-rw-rw-r-- 1 wuzhikun wuzhikun 2.8K Nov 22 16:45 vg_rna.spliced.gcsa.lcp
-rw-rw-r-- 1 wuzhikun wuzhikun  344 Nov 22 16:45 vg_rna.txorigin.tsv

```

```
(RNA2) wuzhikun@cu06 11:11:44 ^_^ /home/wuzhikun/Project/RNA/Test2 
$ vg autoindex --workflow mpmap -t 20 --prefix GRCh38_Gencode  --ref-fasta /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa    --tx-gff /home/wuzhikun/database/GENCODE/gencode.v41.annotation_dechr.gtf
[vg autoindex] Executing command: vg autoindex --workflow mpmap -t 20 --prefix GRCh38_Gencode --ref-fasta /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --tx-gff /home/wuzhikun/database/GENCODE/gencode.v41.annotation_dechr.gtf
error:[vg autoindex] Input is not sufficient to create indexes
Inputs
    GTF/GFF
    Reference FASTA
are insufficient to create target index Haplotype-Transcript GBWT

```



```
grep -v  '^M' /home/wuzhikun/database/GENCODE/gencode.v41.annotation_dechr.gtf >   /home/wuzhikun/database/GENCODE/gencode.v41.annotation_dechr_deM.gtf
```


```
bgzip  /home/wuzhikun/database/1000G_2504_high_coverage/1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf

vg autoindex --workflow mpmap -t 20 --prefix GRCh38_Gencode --ref-fasta /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --tx-gff /home/wuzhikun/database/GENCODE/gencode.v41.annotation_dechr.gtf --vcf /home/wuzhikun/database/1000G_2504_high_coverage/1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf.gz
[vg autoindex] Executing command: vg autoindex --workflow mpmap -t 20 --prefix GRCh38_Gencode --ref-fasta /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --tx-gff /home/wuzhikun/database/GENCODE/gencode.v41.annotation_dechr.gtf --vcf /home/wuzhikun/database/1000G_2504_high_coverage/1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf.gz
```


```
$ vg autoindex --workflow mpmap -t 20 --prefix GRCh38_Gencode_all --ref-fasta /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --tx-gff /home/wuzhikun/database/GENCODE/gencode.v41.annotation_dechr_deM.gtf --vcf /home/wuzhikun/database/1000G_2504_high_coverage/1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_filt_dechr_all.vcf.gz
[vg autoindex] Executing command: vg autoindex --workflow mpmap -t 20 --prefix GRCh38_Gencode --ref-fasta /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --tx-gff /home/wuzhikun/database/GENCODE/gencode.v41.annotation_dechr_deM.gtf --vcf /home/wuzhikun/database/1000G_2504_high_coverage/1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_filt_dechr_all.vcf.gz

```



```

wuzhikun@mu03 07:40:02 ^_^ /home/wuzhikun/database/1000G_2504_high_coverage 
$ cat header 1kGP_high_coverage_Illumina.chr*.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf  > 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_all_chr.vcf
wuzhikun@mu03 07:40:07 ^_^ /home/wuzhikun/database/1000G_2504_high_coverage 
$ bgzip 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_all_chr.vcf
```




```
wuzhikun@mu03 17:52:34 O_O /home/wuzhikun/database/1000G_2504_high_coverage 
$ cat header  1kGP_high_coverage_Illumina.chr1.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr10.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr11.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr12.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr13.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr14.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr15.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr16.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr17.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr18.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr19.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr2.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr20.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr21.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr22.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr3.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr4.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr5.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr6.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr7.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr8.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chr9.filtered.SNV_INDEL_SV_phased_panel_filt_dechr.vcf 1kGP_high_coverage_Illumina.chrX.filtered.SNV_INDEL_SV_phased_panel.v2_filt_dechr.vcf > 1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_all_chr.vcf
```



```
$ vg autoindex --workflow mpmap -t 20 --prefix GRCh38_Gencode_all --ref-fasta /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --tx-gff /home/wuzhikun/database/GENCODE/gencode.v41.annotation_dechr_deM.gtf --vcf /home/wuzhikun/database/1000G_2504_high_coverage/1kGP_high_coverage_Illumina.filtered.SNV_INDEL_SV_phased_panel_all_chr.vcf.gz
```



```
[IndexRegistry]: Exceeded disk use limit while generating k-mers. Rewinding to pruning step with more aggressive pruning to simplify the graph.
[IndexRegistry]: Pruning complex regions of spliced VG to prepare for GCSA indexing with GBWT unfolding.
[IndexRegistry]: Constructing GCSA/LCP indexes.
error: [write_gcsa_kmers()] size limit exceeded
[IndexRegistry]: Exceeded disk use limit while generating k-mers. Rewinding to pruning step with more aggressive pruning to simplify the graph.

```




```
$ vg rna

usage: vg rna [options] graph.[vg|pg|hg|og] > splicing_graph.[vg|pg|hg|og]

General options:
    -t, --threads INT          number of compute threads to use [1]
    -p, --progress             show progress
    -h, --help                 print help message

Input options:
    -n, --transcripts FILE     transcript file(s) in gtf/gff format; may repeat
    -m, --introns FILE         intron file(s) in bed format; may repeat
    -y, --feature-type NAME    parse only this feature type in the gtf/gff (parses all if empty) [exon]
    -s, --transcript-tag NAME  use this attribute tag in the gtf/gff file(s) as id [transcript_id]
    -l, --haplotypes FILE      project transcripts onto haplotypes in GBWT index file

Construction options:
    -j, --use-hap-ref          use haplotype paths in GBWT index as reference sequences (disables projection)
    -e, --proj-embed-paths     project transcripts onto embedded haplotype paths
    -c, --do-not-collapse      do not collapse transcript paths that are identical across references or haplotypes
    -k, --max-node-length      chop nodes longer than maximum node length (0 disables chopping) [0]
    -d, --remove-non-gene      remove intergenic and intronic regions (deletes all paths in the graph)
    -o, --do-not-sort          do not topological sort and compact the graph
    -r, --add-ref-paths        add reference transcripts as embedded paths in the graph
    -a, --add-hap-paths        add haplotype-specific transcripts as embedded paths in the graph

Output options:
    -u, --out-ref-paths        add reference transcript paths to pantranscriptome output
    -b, --write-gbwt FILE      write pantranscriptome transcript paths as GBWT index file
    -g, --gbwt-bidirectional   use bidirectional paths in GBWT index construction
    -f, --write-fasta FILE     write pantranscriptome transcript sequences as fasta file
    -i, --write-info FILE      write pantranscriptome transcript info table as tsv file

```



```
$ vg mpmap
usage: vg mpmap [options] -x graph.xg -g index.gcsa [-f reads1.fq [-f reads2.fq] | -G reads.gam] > aln.gamp
Multipath align reads to a graph.

basic options:
graph/index:
  -x, --graph-name FILE     graph (required; XG format recommended but other formats are valid, see `vg convert`) 
  -g, --gcsa-name FILE      use this GCSA2/LCP index pair for MEMs (required; both FILE and FILE.lcp, see `vg index`)
  -d, --dist-name FILE      use this snarl distance index for clustering (recommended, see `vg index`)
  -s, --snarls FILE         align to alternate paths in these snarls (unnecessary if providing -d, see `vg snarls`)
input:
  -f, --fastq FILE          input FASTQ (possibly gzipped), can be given twice for paired ends (for stdin use -)
  -i, --interleaved         input contains interleaved paired ends
algorithm presets:
  -n, --nt-type TYPE        sequence type preset: 'DNA' for genomic data, 'RNA' for transcriptomic data [RNA]
  -l, --read-length TYPE    read length preset: 'very-short', 'short', or 'long' (approx. <50bp, 50-500bp, and >500bp) [short]
  -e, --error-rate TYPE     error rate preset: 'low' or 'high' (approx. PHRED >20 and <20) [low]
output:
  -F, --output-fmt TYPE     format to output alignments in: 'GAMP for' multipath alignments, 'GAM' or 'GAF' for single-path
                            alignments, 'SAM', 'BAM', or 'CRAM' for linear reference alignments (may also require -S) [GAMP]
  -S, --ref-paths FILE      paths in the graph either 1) one per line in a text file, or 2) in an HTSlib .dict, to treat as
                            reference sequences for HTSlib formats (see -F) [all paths]
  -N, --sample NAME         add this sample name to output
  -R, --read-group NAME     add this read group to output
  -p, --suppress-progress   do not report progress to stderr
computational parameters:
  -t, --threads INT         number of compute threads to use [all available]

advanced options:
algorithm:
  -X, --not-spliced         do not form spliced alignments, even if aligning with --nt-type 'rna'
  -M, --max-multimaps INT   report (up to) this many mappings per read [10 rna / 1 dna]
  -a, --agglomerate-alns    combine separate multipath alignments into one (possibly disconnected) alignment
  -r, --intron-distr FILE   intron length distribution (from scripts/intron_length_distribution.py)
  -Q, --mq-max INT          cap mapping quality estimates at this much [60]
  -b, --frag-sample INT     look for this many unambiguous mappings to estimate the fragment length distribution [1000]
  -I, --frag-mean FLOAT     mean for a pre-determined fragment length distribution (also requires -D)
  -D, --frag-stddev FLOAT   standard deviation for a pre-determined fragment length distribution (also requires -I)
  -G, --gam-input FILE      input GAM (for stdin, use -)
  -u, --map-attempts INT    perform (up to) this many mappings per read (0 for no limit) [24 paired / 64 unpaired]
  -c, --hit-max INT         use at most this many hits for any match seeds (0 for no limit) [1024 DNA / 100 RNA]
scoring:
  -A, --no-qual-adjust      do not perform base quality adjusted alignments even when base qualities are available
  -q, --match INT           use this match score [1]
  -z, --mismatch INT        use this mismatch penalty [4 low error, 1 high error]
  -o, --gap-open INT        use this gap open penalty [6 low error, 1 high error]
  -y, --gap-extend INT      use this gap extension penalty [1]
  -L, --full-l-bonus INT    add this score to alignments that align each end of the read [mismatch+1 short, 0 long]
  -w, --score-matrix FILE   read a 4x4 integer substitution scoring matrix from a file (in the order ACGT)
  -m, --remove-bonuses      remove full length alignment bonuses in reported scores

```


```
vg mpmap -n rna -t 4 -x vg_rna.spliced.xg -g vg_rna.spliced.gcsa -d vg_rna.spliced.dist -f small/x_rna_1.fq -f small/x_rna_2.fq > mpmap.gamp

```

output files:
```
-rw-rw-r-- 1 wuzhikun wuzhikun 584K Nov 22 16:48 mpmap.gamp
```



