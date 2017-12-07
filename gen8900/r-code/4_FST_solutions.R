### Set my working directory to the folder
### where I keep the sample data files for class
setwd("~/Desktop/Teaching/CompGenWS-GEN8900/")

### Get the most up to date version of
### the R functions script
library(devtools)
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)

### Read in the sample data file for this week
my.data = read.vcf("SampleData/sampleData_4_Fst.vcf", header=TRUE, stringsAsFactors=FALSE)

### Set up an empty results table with 9 columns
blank.col = rep(0, nrow(my.data))
my.results = data.frame(p=blank.col, n=blank.col, Hexp=blank.col, p1=blank.col, n1=blank.col, Hexp1=blank.col, p2=blank.col, n2=blank.col, Hexp2=blank.col)

### Determine which columns correspond to each sub-pop
spop1.columns = grep("^FU_", colnames(my.data))
spop2.columns = grep("^MA_", colnames(my.data))

### Since I need to get p, n, and Aa 3 separate times,
### I'm going to write a function for this to make it
### easier to repeat the task over and over
expected.het <- function(genotypes) {
    obs.counts = count.genotypes(genotypes)
    n = sum(obs.counts) - obs.counts["NN"]
    freqs = allele.freq(obs.counts)
    Hexp = 2 * freqs[1] * freqs[2]
    res = c(freqs["p"], n, Hexp)
    res=as.numeric(unname(res))
    return(res)
}

### Run the loop to calc. expAa for all pops
for (i in 1:nrow(my.data)) {
    my.row = as.vector(my.data[i,], mode="character")
    my.samples = my.row[10:length(my.row)]
    my.gens = get.field(my.samples, my.row[9], "GT")
    my.results[i,1:3] = expected.het(my.gens) # Since I made my fxn return all 3 values, I can fill in 3 results columns at once
    my.gens1 = my.gens[(spop1.columns-9)] # Notice that I need to subtract 9 from my original column index, to account for the first 9 columns of the VCF file that are NOT a part of my genotype vector
    my.results[i,4:6] = expected.het(my.gens1)
    my.gens2 = my.gens[(spop2.columns-9)]
    my.results[i,7:9] = expected.het(my.gens2)
}

### Calculate Fst
my.results$Hs = ((my.results$n1/my.results$n) * my.results$Hexp1) + ((my.results$n2/my.results$n) * my.results$Hexp2)
my.results$Fst = (my.results$Hexp - my.results$Hs)/my.results$Hexp

### Plot a Histogram
hist(my.results$Fst, xlim=c(0,1), breaks=50, xlab=expression(F[ST]), main="Histogram of FST", col="lightblue")

### Plot along the chromosome
plot(my.data$POS, my.results$Fst, xlab="Chromosomal Position (bp)", ylab=expression(F[ST]), ylim=c(0,1), pch=20, main="Spatial Distribution of Fst", col="hotpink3")

### Set up an empty matrix to hold numeric genotypes
numericGen = matrix(NA, nrow=nrow(my.data), ncol=(ncol(my.data) - 9))
colnames(numericGen) = colnames(my.data)[10:ncol(my.data)]

### Convert genotype strings to 1, 2, or 3 for PCA
for (i in 1:nrow(my.data)) {
    my.row = as.vector(my.data[i,], mode="character")
    my.samples = samples=my.row[10:length(my.row)]
    my.gens = get.field(samples=my.samples, format=my.row[9], fieldName="GT")
    my.gens = gsub("0/0", 1, gsub("0/1", 2, gsub("1/1", 3, gsub("\\./\\.", 0, my.gens))))
    numericGen[i,] = as.numeric(my.gens)
}

### Transpose the matrix
numGen.trans = t(numericGen)

### Run the pca function
pca.results = prcomp(numGen.trans, scale=FALSE)
summary(pca.results)
plot(pca.results)

### Get a vector of plot colors that correspond to species
plot.col = gsub("^FU_.*", "darkturquoise", gsub("^MA_.*", "darkviolet", rownames(numGen.trans)))

plot(pca.results$x, col=plot.col, pch=20, main="PCA of Finch Data")
legend(18, -1, c("Fuliginosa", "Magnirostris"), col=c("darkturquoise", "darkviolet"), pch=c(20,20), cex=0.75)

