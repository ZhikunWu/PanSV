
## [pbmm2](https://github.com/PacificBiosciences/pbmm2)

### install
```
conda install -c bioconda pbmm2
```

### parameters
```
$ pbmm2 --help
pbmm2 - minimap2 with native PacBio BAM support

Usage:
  pbmm2 <tool>

  -h,--help    Show this help and exit.
  --version    Show application version and exit.

Tools:
  index      Index reference and store as .mmi file
  align      Align PacBio reads to reference sequences

Examples:
  pbmm2 index ref.referenceset.xml ref.mmi
  pbmm2 align ref.referenceset.xml movie.subreadset.xml ref.movie.alignmentset.xml

Typical workflows:
  A. Generate index file for reference and reuse it to align reads
     $ pbmm2 index ref.fasta ref.mmi
     $ pbmm2 align ref.mmi movie.subreads.bam ref.movie.bam

  B. Align reads and sort on-the-fly, with 4 alignment and 2 sort threads
     $ pbmm2 align ref.fasta movie.subreads.bam ref.movie.bam --sort -j 4 -J 2

  C. Align reads, sort on-the-fly, and create PBI
     $ pbmm2 align ref.fasta movie.subreadset.xml ref.movie.alignmentset.xml --sort

  D. Omit output file and stream BAM output to stdout
     $ pbmm2 align hg38.mmi movie1.subreadset.xml | samtools sort > hg38.movie1.sorted.bam

  E. Align CCS fastq input and sort on-the-fly
     $ pbmm2 align ref.fasta movie.Q20.fastq ref.movie.bam --preset CCS --sort --rg '@RG\tID:myid\tSM:mysample'

Copyright (C) 2004-2022     Pacific Biosciences of California, Inc.
This program comes with ABSOLUTELY NO WARRANTY; it is intended for
Research Use Only and not for use in diagnostic procedures.

```



