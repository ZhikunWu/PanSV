#!/usr/bin/env Rscript
library(ggplot2)
library(argparser)
library(readr)

arg <- arg_parser('hist of summary for quality control.')
arg <- add_argument(arg, '--input', help='In file with read quality.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


#usage: Rscript /home/wuzhikun/github/GenomeAssembly/script/LengthBar.R  --input /home/wuzhikun/Project/Vigna/Evaluation/SD_derepeat/final_iden0.9_len_category.txt --pdf /home/wuzhikun/Project/Vigna/Evaluation/SD_derepeat/final_iden0.9_len_category.pdf --width 5 --height 4



Quality_hist <- function(infile, pdf_file, width, height){

    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(infile)
    colnames(data) <- c("Repeat", "Sample", "Value")

    

    ggplot(data, aes(x=Value)) +  
      geom_histogram(aes(y=..density.., fill=Repeat), bins=30, position="identity",  alpha=0.6) +
      geom_density(aes(color=Repeat), alpha=0.7) +
      scale_fill_manual(values=c("tomato3", "forestgreen", "royalblue2")) +
      scale_color_manual(values=c("tomato3", "forestgreen", "royalblue2")) +
      #scale_fill_manual(values=c("#f57c6e", "#b8aeeb", "#84c3b7")) +
      #scale_color_manual(values=c("#f57c6e", "#b8aeeb", "#84c3b7")) +
      xlab("Glucose (mg/g)") + 
      #xlab("Length (kb)") + ylab("Number") +  
      theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
      theme(axis.text = element_text( size=rel(1.1 ))) + 
      theme(axis.title = element_text( size=rel(1.3 ))) +
      theme(legend.title = element_blank())

    ggsave(pdf_file, width=width, height=height)



}


Quality_hist(argv$input, argv$pdf, argv$width, argv$height)


