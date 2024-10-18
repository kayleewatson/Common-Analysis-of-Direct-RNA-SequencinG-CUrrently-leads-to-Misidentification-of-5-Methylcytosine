#!/usr/bin/env Rscript

library(ggplot2)
library(scales)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("Please supply a filename", call.=FALSE)
}

cat("Reading in the data...")
data1 = read.delim(args[1], header=FALSE, sep="\t")

cat("\nFormatting the data...\n\n")
colnames(data1) <- c("Value", "Identifier", "Sample")
data1$Identifier <- as.factor(data1$Identifier)
data1$Sample <- as.factor(data1$Sample)

data1_motif_count <- table(data1$Identifier)

print(data1_motif_count)

cat("\nPlotting...")

ggplot(aes(Value, group = Identifier), data = data1) +
    geom_histogram(binwidth=0.025, colour = 1, fill = "white") +
    geom_rect(aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf),
    colour = "black", fill = NA) +
    facet_wrap(~Identifier, scales="free") +
    scale_color_manual(name = "Sample:", values = c("black", "grey55")) +
    xlab("Methylated Fraction") + theme_linedraw() +
    theme(
    panel.grid.major = element_line(colour = "grey"), panel.grid.minor = element_blank(),
    axis.title=element_blank(),
    legend.text=element_text(size=11),
    legend.title=element_text(size=11),
    strip.text = element_text(size=10, face="bold"),
    axis.text=element_text(size=8),
    legend.position="right",
    strip.background = element_rect(color=NA),
    panel.spacing = unit(0.5, "lines")) +
    scale_y_continuous(labels = label_number(accuracy = 1)) +
    scale_x_continuous(breaks = c(0,0.25,0.5,0.75,1),
    labels = c(0,0.25,0.5,0.75,1))

sample_name = as.character(data1$Sample[1])
filename1 = paste(sample_name, "_histogram.svg", sep="")
ggsave(filename1, width = 7, height = 6, dpi = 600)

cat("\n")
