#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(readr)
library(ggrepel)

#usage: 

arg <- arg_parser('Stack bar plot for total read bases.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output hist file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)


genotypePCA <- function(infile, pdf_file, width, height){
    width <- as.numeric(width)
    height <- as.numeric(height)

    data <- read_tsv(infile)
    print(data)


    # data$Group <- factor(data$Group, order=TRUE, levels=c("Neimenggu","Liaoning","Heilongjiang","Jilin","Hubei","Beijing","Tianjin","Xinjiang","Gansu","Ningxia","Qinghai","Shanxi","Shandong","Guizhou","Henan","Hebei","Hunan","Jiangsu","Guangdong","Guangxi","Unknown"))

    # colors <-c("cadetblue1", "cadetblue4","darkorchid1", "darkorchid3", "lightsteelblue1",  "lightsteelblue3", "darkolivegreen1",  "darkolivegreen4", "cornflowerblue", "blue", "orange", "orange3",  "plum1", "plum3", "springgreen", "springgreen3", "skyblue1", "skyblue4", "brown1", "firebrick4", "grey")

    # colors = c("dodgerblue3", "tomato2", "gray")

    # colors <- c("dodgerblue3", "tomato2", "mediumseagreen", "gold2")

    colors <- c("dodgerblue3", "tomato2", "mediumseagreen","gold2",  "gray")

    #####################################################
    #Target <- c("North", "South")
    #data <- data[which(data$Region %in% Target), ]
    ####################################################


    p<-ggplot(data,aes(x=PC3, y=PC4, label=IID)) + 
      geom_point(fill="blue", alpha=0.4) +
      geom_text(size=2, colour="tomato3", check_overlap=FALSE) +
      #geom_label_repel(point.padding=0.5, segment.color="tomato3") +
      theme( panel.background=element_blank(),panel.grid.major=element_blank(), panel.grid.minor=element_blank(), panel.border = element_blank(), plot.background=element_blank(), axis.line = element_line(colour = "black")) +
       #scale_colour_manual(values=colors) +
       #xlab("PC1 (16.1%)") + ylab("PC2 (8.1%)") + 
       xlab("PC3 (7.5%)") + ylab("PC4 (6.4%)") + 
        theme(axis.text = element_text( size=rel(1.1 ))) +
        theme(axis.title = element_text( size=rel(1.0 )))  +
       theme(legend.text=element_text(size=rel(0.5))) +
       theme(legend.key.width = unit(0.3, "cm"), legend.key.height = unit(0.3, "cm")) +
       theme(legend.title = element_blank()) 
    p

    ggsave(pdf_file, width=width, height=height)
}








genotypePCA(argv$input, argv$pdf, argv$width, argv$height)

