### Set my working directory to the folder
### where I keep the sample data files for class
setwd("~/Desktop/Teaching/CompGenWS-GEN8900/")

### Get the most up to date version of
### the R functions script
library(devtools)
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)

### Read in the sample data file for this week
my.data = read.vcf("SampleData/sampleData_LDmini.vcf", header=TRUE, stringsAsFactors=FALSE)

### Extract the first 2 rows
row1=as.vector(my.data[1,], mode="character")
row2=as.vector(my.data[2,], mode="character")

### Calculate the distance
dist = as.numeric(row2[2]) - as.numeric(row1[2])

### Get pA and pB
g1 = get.field(row1[10:length(row1)], row1[9], "GT") # Using the fxn isn't necessary since I don't have other fields, but I'm going to have it in the code for the future when I might have other stuff there and still want to calc LD
g2 = get.field(row2[10:length(row2)], row2[9], "GT")

pA = unname(allele.freq(count.genotypes(g1))["p"])
pB = unname(allele.freq(count.genotypes(g2))["p"])

### Calculate pAB:
### First, get all of the HAPLOTYPES
get.haplotypes <- function(genotypes1, genotypes2) {
    a1 = gsub("\\|", "", genotypes1) # Get rid of the "|" symbol between alleles
    a2 = gsub("\\|", "", genotypes2)
    a1=unlist(strsplit(paste0(a1, collapse=""), split="")) # Collapse everything into 1 string, then break it apart again but now in a way such that each allele is a separate item in the list
    a2=unlist(strsplit(paste0(a2, collapse=""), split=""))
    haps = paste0(a1,a2)
    return(haps)
}

h=get.haplotypes(g1,g2)

### Now, just count the fraction of times "00" is the hap:
pAB = (length(h[h=="00"]))/length(h)

### Calculate D and r-squared
D = pAB - (pA * pB)
rsq = (D**2)/(pA*(1-pA)*pB*(1-pB))

### Check the other sites:
my.results = data.frame(matrix(0, nrow=10, ncol=5))
colnames(my.results)=c("First Site", "Second Site", "Distance (bp)", "D", "Rsq")
counter=1

for (i in 1:((nrow(my.data))-1)) {
    row1=as.vector(my.data[i,], mode="character")
    g1 = get.field(row1[10:length(row1)], row1[9], "GT")
    pA = unname(allele.freq(count.genotypes(g1))["p"])

    for (j in (i+1):(nrow(my.data))) {
        row2=as.vector(my.data[j,], mode="character")
        g2 = get.field(row2[10:length(row2)], row2[9], "GT")
        pB = unname(allele.freq(count.genotypes(g2))["p"])
        dist = as.numeric(row2[2]) - as.numeric(row1[2])
        h=get.haplotypes(g1,g2)
        pAB = (length(h[h=="00"]))/length(h)
        D = pAB - (pA * pB)
        rsq = (D**2)/(pA*(1-pA)*pB*(1-pB))
        my.results[counter,] = c(i,j,dist, D, rsq)
        counter = counter+1
    }
}
