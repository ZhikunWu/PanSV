

### diamond

```
(Assembly3) wuzhikun@cu11 11:34:41 ^_^ /home/wuzhikun/database/NCBI/NR 
$ diamond --help
diamond v2.0.11.149 (C) Max Planck Society for the Advancement of Science
Documentation, support and updates available at http://www.diamondsearch.org
Please cite: http://dx.doi.org/10.1038/s41592-021-01101-x Nature Methods (2021)

Syntax: diamond COMMAND [OPTIONS]

Commands:
makedb  Build DIAMOND database from a FASTA file
blastp  Align amino acid query sequences against a protein reference database
blastx  Align DNA query sequences against a protein reference database
view    View DIAMOND alignment archive (DAA) formatted file
help    Produce help message
version Display version information
getseq  Retrieve sequences from a DIAMOND database file
dbinfo  Print information about a DIAMOND database file
test    Run regression tests
makeidx Make database index

General options:
--threads (-p)           number of CPU threads
--db (-d)                database file
--out (-o)               output file
--outfmt (-f)            output format
    0   = BLAST pairwise
    5   = BLAST XML
    6   = BLAST tabular
    100 = DIAMOND alignment archive (DAA)
    101 = SAM

    Value 6 may be followed by a space-separated list of these keywords:

    qseqid means Query Seq - id
    qlen means Query sequence length
    sseqid means Subject Seq - id
    sallseqid means All subject Seq - id(s), separated by a ';'
    slen means Subject sequence length
    qstart means Start of alignment in query
    qend means End of alignment in query
    sstart means Start of alignment in subject
    send means End of alignment in subject
    qseq means Aligned part of query sequence
    qseq_translated means Aligned part of query sequence (translated)
    full_qseq means Query sequence
    full_qseq_mate means Query sequence of the mate
    sseq means Aligned part of subject sequence
    full_sseq means Subject sequence
    evalue means Expect value
    bitscore means Bit score
    score means Raw score
    length means Alignment length
    pident means Percentage of identical matches
    nident means Number of identical matches
    mismatch means Number of mismatches
    positive means Number of positive - scoring matches
    gapopen means Number of gap openings
    gaps means Total number of gaps
    ppos means Percentage of positive - scoring matches
    qframe means Query frame
    btop means Blast traceback operations(BTOP)
    cigar means CIGAR string
    staxids means unique Subject Taxonomy ID(s), separated by a ';' (in numerical order)
    sscinames means unique Subject Scientific Name(s), separated by a ';'
    sskingdoms means unique Subject Super Kingdom(s), separated by a ';'
    skingdoms means unique Subject Kingdom(s), separated by a ';'
    sphylums means unique Subject Phylum(s), separated by a ';'
    stitle means Subject Title
    salltitles means All Subject Title(s), separated by a '<>'
    qcovhsp means Query Coverage Per HSP
    scovhsp means Subject Coverage Per HSP
    qtitle means Query title
    qqual means Query quality values for the aligned part of the query
    full_qqual means Query quality values
    qstrand means Query strand

    Default: qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore
--verbose (-v)           verbose console output
--log                    enable debug log
--quiet                  disable console output
--header                 Write header lines to blast tabular format.

Makedb options:
--in                     input reference file in FASTA format
--taxonmap               protein accession to taxid mapping file
--taxonnodes             taxonomy nodes.dmp from NCBI
--taxonnames             taxonomy names.dmp from NCBI

Aligner options:
--query (-q)             input query file
--strand                 query strands to search (both/minus/plus)
--un                     file for unaligned queries
--al                     file or aligned queries
--unfmt                  format of unaligned query file (fasta/fastq)
--alfmt                  format of aligned query file (fasta/fastq)
--unal                   report unaligned queries (0=no, 1=yes)
--max-target-seqs (-k)   maximum number of target sequences to report alignments for (default=25)
--top                    report alignments within this percentage range of top alignment score (overrides --max-target-seqs)
--max-hsps               maximum number of HSPs per target sequence to report for each query (default=1)
--range-culling          restrict hit culling to overlapping query ranges
--compress               compression for output files (0=none, 1=gzip, zstd)
--evalue (-e)            maximum e-value to report alignments (default=0.001)
--min-score              minimum bit score to report alignments (overrides e-value setting)
--id                     minimum identity% to report an alignment
--query-cover            minimum query cover% to report an alignment
--subject-cover          minimum subject cover% to report an alignment
--fast                   enable fast mode
--mid-sensitive          enable mid-sensitive mode
--sensitive              enable sensitive mode)
--more-sensitive         enable more sensitive mode
--very-sensitive         enable very sensitive mode
--ultra-sensitive        enable ultra sensitive mode
--iterate                iterated search with increasing sensitivity
--global-ranking (-g)    number of targets for global ranking
--block-size (-b)        sequence block size in billions of letters (default=2.0)
--index-chunks (-c)      number of chunks for index processing (default=4)
--tmpdir (-t)            directory for temporary files
--parallel-tmpdir        directory for temporary files used by multiprocessing
--gapopen                gap open penalty
--gapextend              gap extension penalty
--frameshift (-F)        frame shift penalty (default=disabled)
--long-reads             short for --range-culling --top 10 -F 15
--matrix                 score matrix for protein alignment (default=BLOSUM62)
--custom-matrix          file containing custom scoring matrix
--comp-based-stats       composition based statistics mode (0-4)
--masking                enable tantan masking of repeat regions (0/1=default)
--query-gencode          genetic code to use to translate query (see user manual)
--salltitles             include full subject titles in DAA file
--sallseqid              include all subject ids in DAA file
--no-self-hits           suppress reporting of identical self hits
--taxonlist              restrict search to list of taxon ids (comma-separated)
--taxon-exclude          exclude list of taxon ids (comma-separated)
--seqidlist              filter the database by list of accessions
--skip-missing-seqids    ignore accessions missing in the database

Advanced options:
--algo                   Seed search algorithm (0=double-indexed/1=query-indexed/ctg=contiguous-seed)
--bin                    number of query bins for seed search
--min-orf (-l)           ignore translated sequences without an open reading frame of at least this length
--freq-sd                number of standard deviations for ignoring frequent seeds
--id2                    minimum number of identities for stage 1 hit
--xdrop (-x)             xdrop for ungapped alignment
--gapped-filter-evalue   E-value threshold for gapped filter (auto)
--band                   band for dynamic programming computation
--shapes (-s)            number of seed shapes (default=all available)
--shape-mask             seed shapes
--multiprocessing        enable distributed-memory parallel processing
--mp-init                initialize multiprocessing run
--mp-recover             enable continuation of interrupted multiprocessing run
--mp-query-chunk         process only a single query chunk as specified
--ext-chunk-size         chunk size for adaptive ranking (default=auto)
--no-ranking             disable ranking heuristic
--ext                    Extension mode (banded-fast/banded-slow/full)
--culling-overlap        minimum range overlap with higher scoring hit to delete a hit (default=50%)
--taxon-k                maximum number of targets to report per species
--range-cover            percentage of query range to be covered for range culling (default=50%)
--dbsize                 effective database size (in letters)
--no-auto-append         disable auto appending of DAA and DMND file extensions
--xml-blord-format       Use gnl|BL_ORD_ID| style format in XML output
--stop-match-score       Set the match score of stop codons against each other.
--tantan-minMaskProb     minimum repeat probability for masking (default=0.9)
--file-buffer-size       file buffer size in bytes (default=67108864)
--memory-limit (-M)      Memory limit for extension stage in GB
--no-unlink              Do not unlink temporary files.
--target-indexed         Enable target-indexed mode
--ignore-warnings        Ignore warnings

View options:
--daa (-a)               DIAMOND alignment archive (DAA) file
--forwardonly            only show alignments of forward strand

Getseq options:
--seq                    Sequence numbers to display.

Online documentation at http://www.diamondsearch.org

```



```
$ diamond blastx --db /home/wuzhikun/database/NCBI/NR/nr.dmnd --out CPG37_ONT_racon3.id_ctg000780.diamond.txt --outfmt 6 qseqid qlen qstart qend  sseqid slen sstart send evalue bitscore  stitle  --threads 20 --query /home/wuzhikun/Project/PanVigna/Assembly/Polish/CPG37/CPG37_ONT_racon3.fasta.split/CPG37_ONT_racon3.id_ctg000780.fasta

```


out file:
```
$ head -n 5   CPG37_ONT_racon3.id_ctg000780.diamond.txt
ctg000780   71362   62173   64086   QCE06174.1  613 1   613 8.53e-213   676 QCE06174.1 brassinosteroid insensitive 1-associated receptor kinase 1 [Vigna unguiculata]
ctg000780   71362   53732   50595   XP_027941542.1  698 1   559 4.14e-192   620 XP_027941542.1 LEAF RUST 10 DISEASE-RESISTANCE LOCUS RECEPTOR-LIKE PROTEIN KINASE-like 2.1 [Vigna unguiculata]
ctg000780   71362   30787   29912   XP_027941540.1  711 1   292 1.39e-190   617 XP_027941540.1 LEAF RUST 10 DISEASE-RESISTANCE LOCUS RECEPTOR-LIKE PROTEIN KINASE-like 2.1 isoform X1 [Vigna unguiculata]
ctg000780   71362   30787   29918   XP_027941541.1  700 1   290 1.78e-190   616 XP_027941541.1 LEAF RUST 10 DISEASE-RESISTANCE LOCUS RECEPTOR-LIKE PROTEIN KINASE-like 2.1 isoform X2 [Vigna unguiculata]
ctg000780   71362   62173   63972   QCE06515.1  582 1   575 4.35e-188   604 QCE06515.1 brassinosteroid insensitive 1-associated receptor kinase 1 [Vigna unguiculata]

```


