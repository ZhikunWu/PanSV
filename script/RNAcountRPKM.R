#!/usr/bin/env Rscript
library(argparser)
library(edgeR)

arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--out', help='The output file.')
argv <- parse_args(arg)



countTable <- read.table(argv$input, sep = '\t', header=TRUE, quote="",check.names = FALSE)

rownames(countTable) <- countTable$Geneid
geneInfo <- countTable[,1:6]
counts <- countTable[, 7]
gene_length <- geneInfo$Length
rpkms <- rpkm(counts, gene_length)
rpkm_out <- cbind(geneInfo, rpkms)

write.table(rpkm_out, argv$out,  row.names=FALSE, sep="\t", quote=FALSE)