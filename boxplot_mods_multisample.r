#!/usr/bin/env Rscript

library(ggplot2)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("Please supply a filename", call.=FALSE)
}

cat("Gathering the data...\n")

df <- read.delim(args[1], header=FALSE, sep="\t")
df$V1 <- as.factor(df$V1)

df_nozero <- df[(df$V3 > 0),]
# removes all '0' modified fractions, so only nonzero fractions are plotted

df_nozero$V2 <- factor(df_nozero$V2, levels = c("All", "Non-GCU", "GCU"))

cat("Plotting...\n")

jpeg("plot.jpeg", width=2000, height=550, units="px")
# create a jpeg file with specified size and units

ggplot(df_nozero, aes(x=V2, y=V3, fill=V2)) +
geom_boxplot(outlier.size = 0.3) +
scale_fill_manual(values=c("white","grey70","grey30")) +
facet_grid(~factor(V1, levels=c("Bmalayi", "Dananassae",
    "Calbicans", "Ecoli", "JW18_SINV", "SINV_IVT"))) +
# create separate boxplot for each organism/sample in a specific order
ylab("Fraction of\nReads Modified") +
theme_bw() +
labs(fill="                                                            ") +
# use a large space for the legend title to shift the legend to the bottom right
    theme(
axis.title.y=element_text(size=45),
axis.title.x=element_blank(),
axis.text.x=element_blank(),
axis.ticks.x = element_blank(),
axis.text.y=element_text(size=35),
strip.text.x = element_text(size=25),
legend.text=element_text(size=50),
legend.title=element_text(size=50),
legend.position = "bottom",
legend.direction = "horizontal",
panel.border = element_rect(colour = "black", fill=NA, size=1.5)) +
scale_x_discrete(limits = c("All", "Non-GCU", "GCU"))

invisible(dev.off())
