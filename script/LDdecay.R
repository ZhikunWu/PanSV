#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)

#usage: Rscript ~/github/PanSV/script/ScaffoldRank.R  --input temp.txt --pdf rank.pdf --width 5 --height 4


arg <- arg_parser('Rank the chromosome length.')
arg <- add_argument(arg, '--input', help='The input file.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


LD_decay <- function(in_file, pdf_file, height, width) {
    height <- as.numeric(height)
    width <- as.numeric(width)
      
    LDValue <- read_tsv(in_file)
    LDValue$Distance <- LDValue$Distance / 1000
       
    colors <-c("limegreen", "royalblue1", "gold2", "tomato", "black")


    decay <-  ggplot(LDValue, aes(x = Distance, y = R2, colour=Group)) +
        geom_line() +
        # geom_point() +
        labs(x = "Distance (kb)", y = expression(italic(r)^2)) +
        xlim(0, 500) +
        theme(axis.text.x = element_text(angle=0,hjust=0.5)) + 
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.text.x=element_blank()) + #axis.ticks=element_blank()
        theme(axis.text.x = element_text(angle = 0, hjust = 0.5)) + #angle = 90, hjust = 1
        theme(plot.margin = margin(1,1,0.5,0, "cm")) +
        scale_fill_manual(values=colors)  + 
        theme(plot.margin = margin(1,1,1,1, "cm")) +
        theme(axis.text = element_text( size=rel(1.3 ))) +
        theme(axis.title = element_text( size=rel(1.4 )))  +
        theme(legend.position = "right")
        
    decay
    ggsave(pdf_file, width=width, height=height)
}

LD_decay(argv$input, argv$pdf, argv$height, argv$width)



