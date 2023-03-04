

## [svim-asm](https://github.com/eldariont/svim-asm)

Structural variant identification method (Assembly edition)


### install
```
conda install -c bioconda svim-asm
```


### parameters

```
(HiC-Pro) wuzhikun@mu03 08:44:05 ^_^ /home/wuzhikun/Project 
$ svim-asm --help
usage: svim-asm [-h] [--version] {haploid,diploid} ...

SVIM-asm (pronounced SWIM-assem) is a structural variant caller for genome-genome alignments. 
It discriminates five different variant classes: deletions, insertions, tandem and interspersed duplications and inversions.
SVIM-asm analyzes alignments between a haploid or diploid query assembly and a reference assembly in SAM/BAM format. 
We recommend to produce the alignments using minimap2.

SVIM-asm has an haploid and a diploid mode depending on the input assembly and performs the following steps:
- COLLECT detects SVs from genome-genome alignments in BAM format
- PAIR merges the SV calls from the two haplotypes of a diploid assembly (diploid mode only)
- OUTPUT prints the found SVs in VCF format

positional arguments:
  {haploid,diploid}  modes
    haploid          Detect SVs from the alignment of an haploid query assembly to a reference assembly
    diploid          Detect SVs from the alignment of a diploid query assembly to a reference assembly

optional arguments:
  -h, --help         show this help message and exit
  --version, -v      show program version number and exit

```




```

(HiC-Pro) wuzhikun@mu03 08:52:21 O_O /home/wuzhikun/Project 
$ svim-asm haploid --help
usage: svim-asm haploid [-h] [--verbose] [--min_mapq MIN_MAPQ] [--min_sv_size MIN_SV_SIZE]
                        [--max_sv_size MAX_SV_SIZE] [--query_gap_tolerance QUERY_GAP_TOLERANCE]
                        [--query_overlap_tolerance QUERY_OVERLAP_TOLERANCE]
                        [--reference_gap_tolerance REFERENCE_GAP_TOLERANCE]
                        [--reference_overlap_tolerance REFERENCE_OVERLAP_TOLERANCE] [--sample SAMPLE]
                        [--types TYPES] [--symbolic_alleles] [--tandem_duplications_as_insertions]
                        [--interspersed_duplications_as_insertions] [--query_names]
                        working_dir bam_file genome

positional arguments:
  working_dir           Working and output directory. Existing files in the directory are overwritten. If the
                        directory does not exist, it is created.
  bam_file              SAM/BAM file with alignment of query assembly to reference assembly (needs to be
                        coordinate-sorted and indexed)
  genome                Reference genome file that the assembly was aligned to (FASTA)

optional arguments:
  -h, --help            show this help message and exit
  --verbose             Enable more verbose logging (default: False)

COLLECT:
  --min_mapq MIN_MAPQ   Minimum mapping quality of alignments to consider (default: 20). Alignments with a
                        lower mapping quality are ignored.
  --min_sv_size MIN_SV_SIZE
                        Minimum SV size to detect (default: 40). SVIM can potentially detect events of any size
                        but is limited by the signal-to-noise ratio in the input alignments. That means that
                        more accurate assemblies and alignments enable the detection of smaller events.
  --max_sv_size MAX_SV_SIZE
                        Maximum SV size to detect (default: 100000). This parameter is used to distinguish long
                        deletions (and inversions) from translocations which cannot be distinguished from the
                        alignment alone. Split read segments mapping far apart on the reference could either
                        indicate a very long deletion (inversion) or a translocation breakpoint. SVIM calls a
                        translocation breakpoint if the mapping distance is larger than this parameter and a
                        deletion (or inversion) if it is smaller or equal.
  --query_gap_tolerance QUERY_GAP_TOLERANCE
                        Maximum tolerated gap between adjacent alignment segments on the query (default: 50).
                        Example: Deletions are detected from two subsequent segments of a split query sequence
                        that are mapped far apart from each other on the reference. The query gap tolerance
                        determines the maximum tolerated length of the query gap between both segments. If
                        there is an unaligned query segment larger than this value between the two segments, no
                        deletion is called.
  --query_overlap_tolerance QUERY_OVERLAP_TOLERANCE
                        Maximum tolerated overlap between adjacent alignment segments on the query (default:
                        50). Example: Deletions are detected from two subsequent segments of a split query
                        sequence that are mapped far apart from each other on the reference. The query overlap
                        tolerance determines the maximum tolerated length of an overlap between both segments
                        in the query. If the overlap between the two segments in the query is larger than this
                        value, no deletion is called.
  --reference_gap_tolerance REFERENCE_GAP_TOLERANCE
                        Maximum tolerated gap between adjacent alignment segments on the reference (default:
                        50). Example: Insertions are detected from two segments of a split query sequence that
                        are mapped right next to each other on the reference but with unaligned sequence
                        between them on the query. The reference gap tolerance determines the maximum tolerated
                        length of the reference gap between both segments. If there is a reference gap larger
                        than this value between the two segments, no insertion is called.
  --reference_overlap_tolerance REFERENCE_OVERLAP_TOLERANCE
                        Maximum tolerated overlap between adjacent alignment segments on the reference
                        (default: 50). Example: Insertions are detected from two segments of a split query
                        sequence that are mapped right next to each other on the reference but with unaligned
                        sequence between them on the query. The reference overlap tolerance determines the
                        maximum tolerated length of an overlap between both segments on the reference. If there
                        is a reference gap larger than this value between the two segments, no insertion is
                        called.

OUTPUT:
  --sample SAMPLE       Sample ID to include in output vcf file (default: Sample)
  --types TYPES         SV types to include in output VCF (default: DEL,INS,INV,DUP:TANDEM,DUP:INT,BND). Give a
                        comma-separated list of SV types. The possible SV types are: DEL (deletions), INS
                        (novel insertions), INV (inversions), DUP:TANDEM (tandem duplications), DUP:INT
                        (interspersed duplications), BND (breakends).
  --symbolic_alleles    Use symbolic alleles, such as <DEL> or <INV> in the VCF output (default: False). By
                        default, deletions, insertions, and inversions are represented by their nucleotide
                        sequence in the output VCF.
  --tandem_duplications_as_insertions
                        Represent tandem duplications as insertions in output VCF (default: False). By default,
                        tandem duplications are represented by the SVTYPE=DUP:TANDEM and the genomic source is
                        given by the POS and END tags. When enabling this option, duplications are instead
                        represented by the SVTYPE=INS and POS and END both give the insertion point of the
                        duplication.
  --interspersed_duplications_as_insertions
                        Represent interspersed duplications as insertions in output VCF (default: False). By
                        default, interspersed duplications are represented by the SVTYPE=DUP:INT and the
                        genomic source is given by the POS and END tags. When enabling this option,
                        duplications are instead represented by the SVTYPE=INS and POS and END both give the
                        insertion point of the duplication.
  --query_names         Output names of supporting query sequences in INFO tag of VCF (default: False). If
                        enabled, the INFO/READS tag contains the list of names of the supporting query
                        sequences.

```


### run svim-asm

```

svim-asm haploid  --sample HG02622_maternal --types DEL,INS,INV,DUP --interspersed_duplications_as_insertions --query_names  /home/wuzhikun/Project/Centromere/SVCall/SvimAsm/HG02622/HG02622_maternal /home/wuzhikun/Project/Centromere/SVCall/mapping/HG02622/HG02622.maternal.bam /home/wuzhikun/database/genome/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa
```


output files:

```
/home/wuzhikun/Project/Centromere/SVCall/SvimAsm/HG02622/HG02622_maternal
├── SVIM_230222_091223.log
├── sv-lengths.png
└── variants.vcf

```



