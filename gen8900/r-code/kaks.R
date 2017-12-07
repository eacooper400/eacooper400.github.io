### Required libraries
library(GenomicFeatures)
library(Rsamtools)
library(VariantAnnotation)

###### EDIT THIS SECTION TO POINT TO YOUR OWN FILES!!!!! #############
### Perform VariantAnnotation pipeline
fa=open(FaFile("~/Desktop/Teaching/CompGenWS-GEN8900/SampleData/9_sampleData_MK/UrsMar_subset.fa"))
txdb=makeTxDbFromGFF(file="~/Desktop/Teaching/CompGenWS-GEN8900/SampleData/9_sampleData_MK/UrsMar_subset.gff3", format="gff3")
vcf=readVcf("~/Desktop/Teaching/CompGenWS-GEN8900/SampleData/9_sampleData_MK/bears.vcf")
effects = predictCoding(vcf, txdb, fa)
#####################################################################


### Find all of the genes present in this sample
my.genes = unique(effects$GENEID)

### Set up a table with 1 column of Gene Names
### and NINE empty columns:
### L0 = freq. nondegenerate sites
### L2 = freq. 2-fold degenerate
### L4 = freq 4-fold degenerate
### A0 = num. Transitions at nondegenerate sites
### A2 = num. Transitions at 2-fold degenerate
### A4 = num. Transitions at 4-fold degenerate
### B0 = num. Transversions at nondegenerate
### B2 = num. Transversions at 2-fold deg.
### B4 = num. Transversions at 4-fold deg.
my.results=data.frame(GeneID=my.genes,L0=0,L2=0,L4=0,A0=0,A2=0,A4=0,B0=0,B2=0,B4=0)

### Separate the effects table by gene, and save in a list
geneEffects=lapply(my.genes, function(x,effects) effects[effects$GENEID==x], effects=effects)

### Get the REFCODONS for each gene
### Then get the VARCODONS separately
ref.codons = lapply(geneEffects, function(x) unname(as.vector(x$REFCODON, mode="character")))
var.codons = lapply(geneEffects, function(x) unname(as.vector(x$VARCODON, mode="character")))

### Remove any codons longer than 3 bp
for (i in seq_along(ref.codons)) {
    to.remove = unique(c(which(nchar(ref.codons[[i]])>3), which(nchar(var.codons[[i]])>3)))
    if (length(to.remove)>0) {
        ref.codons[[i]]=ref.codons[[i]][-to.remove]
        var.codons[[i]]=var.codons[[i]][-to.remove]
    }
}
### Use the function below to classify each site in the gene
### as nondegenerate, 2-fold degen, or 4-fold degen
degen.sites <- function(codon.list) {
    degeneracy=list(T=list(T=list(T=c(3,3,2), C=c(3,3,2), A=c(2,3,2), G=c(2,3,2)),
                           C=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                           A=list(T=c(3,3,2), C=c(3,3,2), A=c(3,2,2), G=c(3,3,2)),
                           G=list(T=c(3,3,2), C=c(3,3,2), A=c(3,2,3), G=c(3,3,3))),
                    C=list(T=list(T=c(2,3,0), C=c(3,3,0), A=c(2,3,0), G=c(2,3,0)),
                           C=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                           A=list(T=c(3,2,2), C=c(3,2,2), A=c(3,2,2), G=c(3,2,2)),
                           G=list(T=c(3,3,0), C=c(3,3,0), A=c(2,3,0), G=c(2,3,0))),
                    A=list(T=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                           C=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                           A=list(T=c(3,2,2), C=c(3,2,2), A=c(3,2,2), G=c(3,2,2)),
                           G=list(T=c(3,3,2), C=c(3,3,2), A=c(2,3,2), G=c(2,3,2))),
                    G=list(T=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                           C=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0)),
                           A=list(T=c(3,2,2), C=c(3,2,2), A=c(3,2,2), G=c(3,2,2)),
                           G=list(T=c(3,3,0), C=c(3,3,0), A=c(3,3,0), G=c(3,3,0))))
    d=unlist(lapply(codon.list,
                    function(x) eval(parse(text=paste0(c("degeneracy",unlist(strsplit(x, split="")),
                                                         use.names=F), collapse="$")))))
    d=as.numeric(gsub("3", "0", gsub(0,4,d)))
    return(d)
}
### Convert each set of ref.codons to a vector of 0s, 2s, and 4s,
### to indicate level of degeneracy for every site
deg.codons=lapply(ref.codons, degen.sites)

### Calculate L0, L2, and L4
my.results$L0 = unlist(lapply(deg.codons, function(x) (length(which(x==0)))/(length(x))))
my.results$L2 = unlist(lapply(deg.codons, function(x) (length(which(x==2)))/(length(x))))
my.results$L4 = unlist(lapply(deg.codons, function(x) (length(which(x==4)))/(length(x))))

### Use the function below to find mutations
### and classify them as Transitins(Ts) or
### Transversions(Tv)
TsTv <- function(ref.codons, var.codons) {
    ref=paste0(ref.codons, collapse="")
    var=paste0(var.codons, collapse="")
    ref=unlist(strsplit(ref, split=""), use.names=F)
    var=unlist(strsplit(var, split=""), use.names=F)
    mut = which(ref != var)
    temp=ref[mut]
    temp[which(ref[mut]=="G")] = "A"
    temp[which(ref[mut]=="A")] = "G"
    temp[which(ref[mut]=="C")] = "T"
    temp[which(ref[mut]=="T")] = "C"
    transV = which(temp != var[mut])
    res = rep("N", length(ref))
    res[mut[transV]] = "Tv"
    res[mut[-transV]] = "Ts"
    return(res)
}

### Now run the function on all genes
tstv.by.gene = vector("list", length(ref.codons))
for (i in seq_along(ref.codons)) {
    tstv.by.gene[[i]]=TsTv(ref.codons[[i]], var.codons[[i]])
}

### Get the counts for A0,A2,A4,B0,B2,and B4 by gene
for (i in seq_along(tstv.by.gene)) {
    tmp=tstv.by.gene[[i]]
    x=deg.codons[[i]]
    my.results$A0[i]=length(intersect(which(tmp=="Ts"), which(x==0)))
    my.results$A2[i]=length(intersect(which(tmp=="Ts"), which(x==2)))
    my.results$A4[i]=length(intersect(which(tmp=="Ts"), which(x==4)))
    my.results$B0[i]=length(intersect(which(tmp=="Tv"), which(x==0)))
    my.results$B2[i]=length(intersect(which(tmp=="Tv"), which(x==2)))
    my.results$B4[i]=length(intersect(which(tmp=="Tv"), which(x==4)))
}

### Use Kimura's 2-parameter model equation to calculate Ka and Ks
my.results$Ka = my.results$A0 +
    (((my.results$L0 * my.results$B0) + (my.results$L2*my.results$B2))/(my.results$L0 + my.results$L2))
my.results$Ks = my.results$B4 +
    (((my.results$L2 * my.results$A2) + (my.results$L4*my.results$A4))/(my.results$L2 + my.results$L4))

### And finally, the KaKs ratio:
my.results$KaKs = my.results$Ka/my.results$Ks
