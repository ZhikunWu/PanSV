
#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(reshape2)
library(readr)


#usage: Rscript ~/github/NanoHub/script/SVTypeBarMultiple.R --input Samples_SV_summary.xls --hist hist.pdf --box box.pdf --width 8 --height 4

arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)



 
multipleStackBar <- function(input, pdf_file, height, width){

    in_file <- argv$input
    summary_data <- read_tsv(in_file)
    melt_dt <- reshape2::melt(summary_data)
    colnames(melt_dt) <- c("Sample", "Type", "Value")
    melt_dt$Type <- factor(melt_dt$Type, levels=c("Core", "Softcore", "Dispensable", "Private"))
    
    colors <-c("limegreen", "royalblue1", "gold2", "tomato",  "purple2",  "black")

    p <- ggplot(melt_dt, aes(x=reorder(Sample, -Value), y=Value, fill=Type)) + 
        # theme_bw() +
        geom_bar(stat="identity",width=0.7,) +
        xlab("Samples") + ylab("Number") +
        theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
        theme(axis.ticks = element_line(size=rel(0.5))) + 
        theme(axis.text.x = element_text(angle = 90, vjust =0.5, hjust = 0.5, size=rel(0.7))) +
        # theme(axis.ticks.x = element_blank()) + 
        # theme(axis.text.x = element_blank()) +
        theme(plot.margin = margin(0.3,0.5,0.5,0.3, "cm")) +
        scale_fill_manual(values=colors) +
        theme(axis.text.y = element_text( size=12 )) + #size=rel(1.5)
        theme(axis.title= element_text(size=18)) +
        theme(legend.box = "horizontal", legend.position = "top") + #legend.position = c(0.50, 0.85)
        # theme(legend.position = "right") + # bottom
        theme(legend.text=element_text(size=rel(0.7))) +
        theme(legend.key.width = unit(0.5, "cm"), legend.key.height = unit(0.5, "cm")) +
        theme(legend.title = element_blank()) 

    p

    height <- as.numeric(height)
    width <- as.numeric(width)
    ggsave(pdf_file, width=width, height=height)
}

 


multipleStackBar(argv$input, argv$pdf, argv$height, argv$width)

