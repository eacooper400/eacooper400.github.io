### Set the working directory
setwd("~/Desktop/Teaching/CompGenWS-GEN8900/")

### Load the Bioconductor scripts
source("https://bioconductor.org/biocLite.R")

### Install the packages for Genome Annotation
biocLite("VariantAnnotation")
biocLite("GenomicFeatures")
biocLite("Rsamtools")

### Load in the Reference Genome Data
library(Rsamtools)
fa = open(FaFile("SampleData/9_sampleData_MK/UrsMar_subset.fa"))

### Load in the Gene Table Information
library(GenomicFeatures)
txdb=makeTxDbFromGFF(file="SampleData/9_sampleData_MK/UrsMar_subset.gff3", format="gff3")

### Read in the VCF file with the
### VariantAnnotation package's function
library(VariantAnnotation)
vcf.object=readVcf(file="SampleData/9_sampleData_MK/bears.vcf")

### Use the predictCoding function to determine
### the effect of each SNP in the VCF file
### Save the GeneID that each SNP is found in as
### the ID column; save the effect as a CSQ in the INFO column
effects = predictCoding(vcf.object, txdb, fa)
id.from.effects = names(ranges(effects))
id.from.vcf = rownames(info(vcf.object))
m = match(id.from.vcf, id.from.effects)
info(vcf.object)$CSQ = effects$CONSEQUENCE[m]
names(vcf.object) = effects$GENEID[m]

### Read in the my functions file
library(devtools)
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)

### Simultaneously read out the annotated VCF file,
### and read it into a "my.data" table with our usual function
my.data=read.vcf(writeVcf(vcf.object, "SampleData/9_sampleData_MK/bears_wAnn.vcf"), header=TRUE, stringsAsFactors=FALSE)

### Remove SNPs in non-coding regions by getting rid of all
### rows with an "NA" value in the ID column
my.data=my.data[!(is.na(my.data$ID)),]

### Clean up the INFO column to only have the Effect (and not the other info)
my.data$INFO = gsub(".*;CSQ=", "", my.data$INFO)

### Set up a results table to categorize
### sites as Fixed or Polymorphic
new=my.data[,c(3,8)]

### Add a column for the allele frequency for each group
new$Brown_refAF = NA
new$Polar_refAF = NA

### Split the data into 2 data frames:
### one for brown bears and one for polar bears
brown.col=grep("^Ua", colnames(my.data))
polar.col=grep("^Umar", colnames(my.data))
brown.data=my.data[,-polar.col]
polar.data=my.data[,-brown.col]

### Calculate allele freq (p) in the Brown bears
for (i in 1:nrow(brown.data)) {
    my.row=as.vector(brown.data[i,], mode="character")
    genotypes=get.field(my.row[10:length(my.row)], my.row[9], "GT")
    all.counts=count.genotypes(genotypes)
    ref.freq=allele.freq(all.counts)["p"]
    new$Brown_refAF[i]=ref.freq
}

### Calculate allele freq(p) in the Polar bears
for (i in 1:nrow(polar.data)) {
     my.row=as.vector(polar.data[i,], mode="character")
     genotypes=get.field(my.row[10:length(my.row)], my.row[9], "GT")
     all.counts=count.genotypes(genotypes)
     ref.freq=allele.freq(all.counts)["p"]
     new$Polar_refAF[i]=ref.freq
}

### Calculate the allele frequency difference b/t the 2
### species, and assign each site a category of "Fixed" if the
### difference is 1, and "Polymorphic" otherwise
new$Diff = abs(new$Brown_refAF - new$Polar_refAF)
new$SiteClass = new$Diff
new$SiteClass[which(new$Diff==1)] <- "Fixed"
new$SiteClass[which(new$Diff<1)] <- "Polymorphic"

### Get the List of Genes
geneList=unique(my.data$ID)

### Set up an empty results table
my.results = data.frame(Genes=geneList, Fn=rep(0, length(geneList)), Fs=rep(0, length(geneList)), Pn=rep(0, length(geneList)), Ps=rep(0, length(geneList)))

### Loop through the gene list,
### get the subset of data for each gene,
### count the fixed and polymorphic non-syn and
### synonymous sites
for (i in 1:length(geneList)) {
    gene=new[which(new$ID==geneList[i]),]
    syn=which(gene$INFO=="synonymous")
    non=which(gene$INFO!="synonymous")
    my.results$Fn[i] = length(intersect(non, which(gene$SiteClass=="Fixed")))
    my.results$Fs[i] = length(which(syn %in% which(gene$SiteClass=="Fixed")))
    my.results$Pn[i] = nrow(subset(gene, (gene$SiteClass=="Polymorphic" & gene$INFO!="synonymous")))
    my.results$Ps[i] = nrow(gene) - sum(my.results[i,2:5])
}
    
### Calculate N:S for Between and Within Species
my.results$BetweenN_S = my.results$Fn/my.results$Fs
my.results$WithinN_S = my.results$Pn/my.results$Ps

### Calculate the Neutrality Index:
my.results$NI = (my.results$Pn/my.results$Ps)/(my.results$Fn/my.results$Fs)

### Get a p-value for each gene
my.results$pValues=apply(my.results[,2:5], 1, function(x) fisher.test(matrix(x, ncol=2))$p.value)

### Plot the results
for.plot=matrix(c(my.results$BetweenN_S, my.results$WithinN_S), ncol=2)
colnames(for.plot)=c("Between", "Within")
rownames(for.plot)=my.results$Genes

barplot(t(for.plot), beside=T, ylab="Ratio of Nonsynonymous to Synonymous Sites", xlab="Gene", col=c("darkgreen", "dodgerblue"), ylim=c(0,6))
legend("topleft", c("Between Species (F)", "Within Species (P)"), pch=15, col=c("darkgreen", "dodgerblue"), bty="n")

### Add notations to show p-values
### First "bracket":
lines(c(1.5,2.5), c(3,3))
lines(c(1.5,1.5), c(3,2.9))
lines(c(2.5,2.5), c(3,2.9))

### Bracket 2:
lines(c(4.5,5.5), c(4.5,4.5))
lines(c(4.5,4.5), c(4.5,4.4))
lines(c(5.5,5.5), c(4.5,4.4))

lines(c(4.5,5.5), c(4.5,4.5))
lines(c(4.5,4.5), c(4.5,4.4))
lines(c(5.5,5.5), c(4.5,4.4))

temp=rep("p =", 2)
p.labels=paste(temp, round(my.results$pValues, 2))
text(c(2,5), c(3.1,4.6), p.labels)





