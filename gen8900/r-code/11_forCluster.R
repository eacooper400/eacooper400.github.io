### Set the Library Path to Install and
### Load any Packages; install devtools
.libPaths(c('~/compGenWS_112017/Rlibs', .libPaths()))
install.packages("devtools", repos='http://cran.us.r-project.org', lib='~/compGenWS_112017/Rlibs', dependencies=TRUE)

### Use command line arguments to specify input and output files,
### as well as the snp positions to look at
args = commandArgs(trailingOnly=TRUE)
input = args[1]
output = args[2]
snp1 = as.numeric(args[3])
snp2 = as.numeric(args[4])

### Read in my functions file,
### set the strings as factors = FALSE
library(devtools)
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)
options(stringsAsFactors=FALSE)

### Read in the sample data, then find the positions associated with the 2 SNPs
my.data=read.vcf(input, header=TRUE)
snp.of.interest = my.data[which(my.data$POS==snp1),]
control.snp = my.data[which(my.data$POS==snp2),]

### Calculate R-squared between the SNP of
### interest and ALL other SNPs in the data
my.results = data.frame(Pos=my.data$POS, Rsq=0)
my.results$Rsq = apply(my.data, 1, function(x,y)
    calc_r2(as.vector(x, mode="character"), as.vector(y, mode="character")), y=snp.of.interest)

### Get the average LD in 10Kb bins
BinStarts=seq(1, max(my.results$Pos), 10000)
smoothed.results = data.frame(BinStarts, BinEnds=BinStarts+10000, AvgRsq=0)
for (i in 1:nrow(smoothed.results)) {
    w = subset(my.results, (my.results$Pos>=smoothed.results$BinStarts[i]
        & my.results$Pos<smoothed.results$BinEnds[i]))
    smoothed.results$AvgRsq[i] = mean(w$Rsq, na.rm=TRUE)
}
### Now, pick a control SNP and do the same thing
my.results$Rsq2 = apply(my.data, 1, function(x,y)
    calc_r2(as.vector(x, mode="character"), as.vector(y, mode="character")), y=control.snp)
smoothed.results$AvgRsq2=0
for (i in 1:nrow(smoothed.results)) {
    w = subset(my.results, (my.results$Pos>=smoothed.results$BinStarts[i]
        & my.results$Pos<smoothed.results$BinEnds[i]))
    smoothed.results$AvgRsq2[i] = mean(w$Rsq2, na.rm=TRUE)
}

### Save the results to an output file
write.table(smoothed.results, paste0(output, "txt", collapse="."), quote=FALSE,
            row.names=FALSE, col.names=TRUE, sep="\t")

### Open a PDF file to plot the results in
pdf(paste0(output, "pdf", collapse="."))

### Plot the results,
### add a star where the SNP of interest is
plot(smoothed.results$BinStarts, smoothed.results$AvgRsq,
     xlab="Position (bp)", ylab="Linkage Disequilibrium (r2)",
     type="l", col="red", ylim=c(0,1))
points(snp.of.interest[,2], 1, pch=8, col="red", cex=2)
lines(smoothed.results$BinStarts, smoothed.results$AvgRsq2, col="blue")
points(control.snp[,2], 1, pch=8, col="blue", cex=1.5)

### Close the connection to the pdf file!
dev.off()
