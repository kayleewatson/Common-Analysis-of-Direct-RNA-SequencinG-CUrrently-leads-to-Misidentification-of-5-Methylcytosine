library(dplyr)

args = commandArgs(trailingOnly=TRUE)

if (length(args)==0) {
    stop("Please supply a filename", call.=FALSE)
}
# import bed files
df1 <- read.delim(args[1],header=FALSE, sep="\t")
df1$Sample <- as.factor(df1$V1)
df1$Motif <- as.factor(df1$V2)
df1$Modified_Fraction <- df1$V3

df2 <- read.delim(args[2],header=FALSE, sep="\t")
df2$Sample <- as.factor(df2$V1)
df2$Motif <- as.factor(df2$V2)
df2$Modified_Fraction <- df2$V3

delta_0 <- 0
# calculate the variance for modified fractions
sigma_sq_1 <- var(df1$Modified_Fraction)
sigma_sq_2 <- var(df2$Modified_Fraction)

n_1 <- length(df1$Modified_Fraction)
n_2 <- length(df2$Modified_Fraction)

# perform a z-test
z_stat <- (mean(df1$Modified_Fraction) - mean(df2$Modified_Fraction) - delta_0) /
    sqrt(sigma_sq_1 / n_1 + sigma_sq_2 / n_2)

print(paste0("z stat: ", z_stat))

# calculate and print the p-value from the z-test
pv <- 2*pnorm(q=abs(z_stat),lower.tail=FALSE,log.p=TRUE)
print(paste0("Log p-value: ",pv))
