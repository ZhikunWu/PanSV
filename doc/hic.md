

```
hicSumMatrices \
--matrices hic/N_Ben_HiC2..fq.gz.bam-hic_matrix.h5  hic/N_Ben_HiC4..fq.gz.bam-hic_matrix.h5  \
--outFileName hic/replicateMerged.h5

hicMergeMatrixBins \
--matrix hic/replicateMerged.h5 --numBins 100 \
--outFileName hic/replicateMerged.100bins.h5

hicPlotMatrix \
--matrix hic/replicateMerged.100bins.h5 \
--log1p \
--dpi 300 \
--clearMaskedBins \
--chromosomeOrder chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
--colorMap jet \
--title "Hi-C matrix" \
--outFileName hic/plot_1Mb_matrix.png

hicMergeMatrixBins \
--matrix hic/replicateMerged.h5 --numBins 2 \
--outFileName hic/replicateMerged.matrix_20kb.h5

hicCorrectMatrix diagnostic_plot \
--chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
--matrix hic/replicateMerged.matrix_20kb.h5 --plotName hic/diagnostic_plot.png

hicCorrectMatrix correct \
--chromosomes chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
--matrix hic/replicateMerged.matrix_20kb.h5 \
--filterThreshold -2 3 --perchr --outFileName hic/replicateMerged.Corrected_20kb.h5

hicPlotMatrix \
--matrix hic/replicateMerged.Corrected_20kb.h5 \
--log1p \
--dpi 300 \
--clearMaskedBins \
--chromosomeOrder chr1 chr2 chr3 chr4 chr5 chr6 chr7 chr8 chr9 chr10 chr11 chr12 chr13 chr14 chr15 chr16 chr17 chr18 chr19 \
--colorMap jet \
--title "Corrected Hi-C matrix" \
--outFileName hic/plot_1Mb_corrected_matrix.png

hicFindTADs --matrix hic/replicateMerged.Corrected_20kb.h5 \
--minDepth 60000 --maxDepth 120000 --numberOfProcessors 8 --step 20000 \
--outPrefix hic/marks_et-al_TADs_20kb-Bins  --minBoundaryDistance 80000 \
--correctForMultipleTesting fdr --threshold 0.05
```

