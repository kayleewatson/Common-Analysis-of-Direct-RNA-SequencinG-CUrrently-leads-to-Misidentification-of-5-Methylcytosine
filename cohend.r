#!/usr/bin/env Rscript

library(psych)

args = commandArgs(trailingOnly=TRUE)
if (length(args)==0) {
    stop("Please supply a filename", call.=FALSE)
}

cat("Gathering the data...\n")

# import bed files
df1 <- read.delim(args[1], header=FALSE, sep="\t")
df2 <- read.delim(args[2], header=FALSE, sep="\t")

# option to perform a t-test (just uncomment)
# testing on the third column of the dataframe, this can be changed
#t.test(df1$V3, df2$V3, paired=FALSE)

# combine the two bed files
df_all <- rbind(df1, df2)
df_all$V1 = NULL

# perform cohen's d for effect size
cohen.d(df_all, "V2", alpha=.05)
