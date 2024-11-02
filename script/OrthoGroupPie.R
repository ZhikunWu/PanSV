
#!/usr/bin/Rscript
library(ggplot2)
library(argparser)
library(readr)
library(dplyr)

#usage: Rscript ~/github/NanoHub/script/SVTypeBarMultiple.R --input Samples_SV_summary.xls --hist hist.pdf --box box.pdf --width 8 --height 4

arg <- arg_parser('Stack bar plot for structural variant types.')
arg <- add_argument(arg, '--input', help='The file with type and number of SV.')
arg <- add_argument(arg, '--pdf', help='output file with pdf format.')
arg <- add_argument(arg, '--width', help='The width of picture.')
arg <- add_argument(arg, '--height', help='The height of picture.')
argv <- parse_args(arg)



pie_plot <- function(in_file, pdf_file, height, width){
    in_file <- argv$input
    data <- read_tsv(in_file)
    height <- as.numeric(height)
    width <- as.numeric(width)


    data <- data %>%
	    arrange(desc(Category)) %>%
	    mutate(prop = Number / sum(data$Number) * 100) %>%
	    mutate(ypos = cumsum(prop) - 0.5*prop)

    p <- ggplot(data, aes(x="", y=prop, fill=Category)) +
        geom_bar(stat="identity", width=1, color="white") +
        coord_polar("y", start=0) +
        theme_void() +
        theme(legend.position="none") +
        geom_text(aes(y=ypos, label=Category), color="white", size=3) +
        scale_fill_brewer(palette="Set1")
    p
    ggsave(pdf_file, width=width, height=height)
}



pie_plot(argv$input, argv$pdf, argv$height, argv$width)


