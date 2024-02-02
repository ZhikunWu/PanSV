#!/usr/bin/env Rscript
# library(CMplot)
### modified CMplot.r
source("/home/wuzhikun/github/SCZSV/script/CMplot.r") 
library(argparser)


### package
### https://cloud.r-project.org/src/contrib/CMplot_3.3.3.tar.gz

#usage: Rscript ~/github/TrioWGS/script/SVDistribution.R --input /home/wuzhikun/Project/NanoTrio/SVStats/Sniffles/minimap2/M628-2_position.xls --pdf /home/wuzhikun/Project/NanoTrio/SVStats/Sniffles/minimap2/M628-2_position.pdf


arg <- arg_parser('Density plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--maxNum', help='The max number for heatmap.')
argv <- parse_args(arg)



SNV_distribution <- function(in_file, out_file, maxNum){


    SNV_pos <- read.table(in_file, sep="\t", header=FALSE)
    colnames(SNV_pos) <- c("Chromosome", "Position" )
    SNV_pos$SNP <- paste(SNV_pos$Chromosome, SNV_pos$Position, sep="_")

    SNV_position <- data.frame(SNP=SNV_pos$SNP, Chromosome=SNV_pos$Chromosome, Position=SNV_pos$Position)

    #chroms <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22, "X")
    chroms <- c("Chr01", "Chr02", "Chr03", "Chr04", "Chr05", "Chr06", "Chr07", "Chr08", "Chr09", "Chr10", "Chr11")
    New_position <- SNV_position[SNV_position$Chromosome %in% chroms,]

    maxNum <- as.numeric(maxNum)
    CMplot(New_position, outName=out_file, maxNum=maxNum, plot.type="d",bin.size=1e5,col=c("darkgreen", "gold2", "red"),file="pdf",memo="",dpi=300, file.output=TRUE, verbose=TRUE, chr.labels=NULL)

}


SNV_distribution(argv$input, argv$pdf, argv$maxNum)
