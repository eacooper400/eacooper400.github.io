### Set my working directory
setwd("~/Desktop/Teaching/CompGenWS-GEN8900/")

### Read in the my functions file
library(devtools)
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)

### Read in the sample VCF file
my.data = read.vcf("SampleData/7_sampleData_AFS.vcf", header=TRUE, stringsAsFactors=FALSE)

### Calculate the derived allele count at each site
derived.counts = rep(0, nrow(my.data)) # First, I set up an vector of 0s the same length as my number of rows to hold results.

for (i in 1:nrow(my.data)) {
    my.row = as.vector(my.data[i,], mode="character")
    ## This next line combines the get.field function with the count.genotype function
    ## to return the genotype counts in 1 command (get.field parses the vcf row to get a vector of
    ## 0/0, 0/1, or 1/1 genotypes, the count.genotypes takes that vector and counts the number of each type
    gen.counts = count.genotypes(get.field(my.row[10:length(my.row)], my.row[9], "GT"))
    ## The vector from count.genotypes is "named," so I can extract the values I want with the notation below
    derived = (2*gen.counts["aa"]) + gen.counts["Aa"]
    derived.counts[i] = derived # lastly, I put my count into my results vector
}

### Create a histogram of the counts, and save the info.
N = ncol(my.data) - 9 # Count the total number of samples/individuals
n.chr = 2*N # multiply the number of samples by 2, to get the number of chromosomes
possible.counts = seq(1, n.chr, by=1) # Use seq to get a vector of every number between 1 and n.chr

histinfo = hist(derived.counts, breaks=possible.counts) # Save the histogram information in a variable
plot(histinfo, col=c("lightblue")) # Check the plot output

### Calculate Waterson's theta
Sn = nrow(my.data) # get the total number of SNPs (equivalent to the number of rows in the file)

### For the denominator, I know that my i values need to go from 1 to 2N-1.  1 to 2N is the range that I had to use for my histogram breakpoints, so I am going to re-use that possible.counts vector and subtract the very last value (48).
### To get 1 over every one of those values, I can rely on the fact that R knows how to do vectorized math, and just type 1/possible.counts.
### Then, I get the sum of all of those values, and that is my denominator.
theta.w = Sn/(sum(1/(possible.counts[-length(possible.counts)]))) 

### Now, get the expected distribution of counts based on theta
exp.counts = theta.w * (1/possible.counts[-length(possible.counts)])

compare.counts=matrix(c(histinfo$counts, exp.counts), ncol=2, byrow=FALSE) # get my observed counts (from histinfo) and expected counts into a 2 column matrix
colnames(compare.counts)=c("Obs", "Exp") # assign names to the columns
barplot(t(compare.counts), beside=T, col=c("lightblue", "black")) # use barplot, with the beside option set to "TRUE" to plot values next to each other
legend("topright", c("Obs", "Exp"), pch=c(15,15), col=c("lightblue", "black")) # add the legend
