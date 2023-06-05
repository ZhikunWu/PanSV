
### SQANTI3 parameters

```
(RNA2) wuzhikun@cu01 19:54:26 ^_^ /home/wuzhikun/data/bean/A3-15/Nanopore/A3-15-1 
$ python /home/wuzhikun/software/SQANTI3-4.3/sqanti3_qc.py --help
R scripting front-end version 4.1.3 (2022-03-10)
usage: sqanti3_qc.py [-h] [--min_ref_len MIN_REF_LEN] [--force_id_ignore]
                     [--aligner_choice {minimap2,deSALT,gmap,uLTRA}]
                     [--cage_peak CAGE_PEAK]
                     [--polyA_motif_list POLYA_MOTIF_LIST]
                     [--polyA_peak POLYA_PEAK] [--phyloP_bed PHYLOP_BED]
                     [--skipORF] [--is_fusion] [--orf_input ORF_INPUT]
                     [--fasta] [-e EXPRESSION] [-x GMAP_INDEX] [-t CPUS]
                     [-n CHUNKS] [-o OUTPUT] [-d DIR] [-c COVERAGE] [-s SITES]
                     [-w WINDOW] [--genename] [-fl FL_COUNT] [-v]
                     [--saturation] [--report {html,pdf,both,skip}]
                     [--isoAnnotLite] [--gff3 GFF3]
                     [--short_reads SHORT_READS] [--SR_bam SR_BAM]
                     isoforms annotation genome

Structural and Quality Annotation of Novel Transcript Isoforms

positional arguments:
  isoforms              Isoforms (FASTA/FASTQ) or GTF format. It is
                        recommended to provide them in GTF format, but if it
                        is needed to map the sequences to the genome use a
                        FASTA/FASTQ file with the --fasta option.
  annotation            Reference annotation file (GTF format)
  genome                Reference genome (Fasta format)

optional arguments:
  -h, --help            show this help message and exit
  --min_ref_len MIN_REF_LEN
                        Minimum reference transcript length (default: 200 bp)
  --force_id_ignore     Allow the usage of transcript IDs non related with
                        PacBio's nomenclature (PB.X.Y)
  --aligner_choice {minimap2,deSALT,gmap,uLTRA}
  --cage_peak CAGE_PEAK
                        FANTOM5 Cage Peak (BED format, optional)
  --polyA_motif_list POLYA_MOTIF_LIST
                        Ranked list of polyA motifs (text, optional)
  --polyA_peak POLYA_PEAK
                        PolyA Peak (BED format, optional)
  --phyloP_bed PHYLOP_BED
                        PhyloP BED for conservation score (BED, optional)
  --skipORF             Skip ORF prediction (to save time)
  --is_fusion           Input are fusion isoforms, must supply GTF as input
  --orf_input ORF_INPUT
                        Input fasta to run ORF on. By default, ORF is run on
                        genome-corrected fasta - this overrides it. If input
                        is fusion (--is_fusion), this must be provided for ORF
                        prediction.
  --fasta               Use when running SQANTI by using as input a
                        FASTA/FASTQ with the sequences of isoforms
  -e EXPRESSION, --expression EXPRESSION
                        Expression matrix (supported: Kallisto tsv)
  -x GMAP_INDEX, --gmap_index GMAP_INDEX
                        Path and prefix of the reference index created by
                        gmap_build. Mandatory if using GMAP unless -g option
                        is specified.
  -t CPUS, --cpus CPUS  Number of threads used during alignment by aligners.
                        (default: 10)
  -n CHUNKS, --chunks CHUNKS
                        Number of chunks to split SQANTI3 analysis in for
                        speed up (default: 1).
  -o OUTPUT, --output OUTPUT
                        Prefix for output files.
  -d DIR, --dir DIR     Directory for output files. Default: Directory where
                        the script was run.
  -c COVERAGE, --coverage COVERAGE
                        Junction coverage files (provide a single file, comma-
                        delmited filenames, or a file pattern, ex:
                        "mydir/*.junctions").
  -s SITES, --sites SITES
                        Set of splice sites to be considered as canonical
                        (comma-separated list of splice sites). Default:
                        GTAG,GCAG,ATAC.
  -w WINDOW, --window WINDOW
                        Size of the window in the genomic DNA screened for
                        Adenine content downstream of TTS
  --genename            Use gene_name tag from GTF to define genes. Default:
                        gene_id used to define genes
  -fl FL_COUNT, --fl_count FL_COUNT
                        Full-length PacBio abundance file
  -v, --version         Display program version number.
  --saturation          Include saturation curves into report
  --report {html,pdf,both,skip}
                        select report format --html --pdf --both --skip
  --isoAnnotLite        Run isoAnnot Lite to output a tappAS-compatible gff3
                        file
  --gff3 GFF3           Precomputed tappAS species specific GFF3 file. It will
                        serve as reference to transfer functional attributes
  --short_reads SHORT_READS
                        File Of File Names (fofn, space separated) with paths
                        to FASTA or FASTQ from Short-Read RNA-Seq. If
                        expression or coverage files are not provided,
                        Kallisto (just for pair-end data) and STAR,
                        respectively, will be run to calculate them.
  --SR_bam SR_BAM       Directory or fofn file with the sorted bam files of
                        Short Reads RNA-Seq mapped against the genome

```


```
python /home/wuzhikun/software/SQANTI3-4.3/sqanti3_qc.py /home/wuzhikun/software/SQANTI3-4.3/example/UHR_chr22.gtf /home/wuzhikun/database/GENCODE/gencode.v41.annotation.gtf  /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly_chr.fa
```


```
python /home/wuzhikun/software/SQANTI3-4.3/sqanti3_qc.py  --aligner_choice minimap2 --fasta /home/wuzhikun/Project/RNA/flair/Trans_collapse.isoforms.fa /home/wuzhikun/database/GENCODE/gencode.v41.annotation.gtf  /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly_chr.fa


Invalid input IDs! Expected PB.X.Y or PB.X.Y|xxxxx or PBfusion.X format but saw ERR2856519.1348315;16 instead. Abort!
```


