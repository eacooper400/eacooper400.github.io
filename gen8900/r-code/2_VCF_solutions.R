### Solutions for Week 2 Exercises

### Set the working directory to the folder where I keep my files for this class
### setwd("~/Desktop/Teaching/CompGenWS-GEN8900/SampleData")
setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Desktop/Teaching/CompGenWS-GEN8900/SampleData")

### First, I create a function that can read in VCF files
### (Discussed in class)
read.vcf <- function(file, special.char="##", ...) {
    my.search.term=paste0(special.char, ".*")
    all.lines=readLines(file)
    clean.lines=gsub(my.search.term, "",  all.lines)
    clean.lines=gsub("#CHROM", "CHROM", clean.lines)
    read.table(..., text=paste(clean.lines, collapse="\n"))
}

### Now, I can use the function I just created to read in
### the sample VCF file:
my.data=read.vcf(file="sampleData_wk2.vcf", header=TRUE, stringsAsFactors=FALSE)

#################################################################################
### Problem 1: Count the number of AA, Aa, aa, and N genotypes at each position
#################################################################################

### For my solution, I am going to start by creating a couple of functions
### that will help me extract JUST the genotype info. from each row

### This first function will parse the long string of sample
### information and return a specific field specified by fieldName
### The reason I created this function (instead of just one that gets or counts
### genotypes in one line) is because at some point in the future I might want to be
### able to extract a different piece from the sample info., so now I have a versatile fxn
get.field <- function(samples, format, fieldName) {
    x=strsplit(samples, split=":") # This breaks apart the string into a vector of things that were separated by ":"
    fields=unlist(strsplit(format, split=":")) # This breaks apart the FORMAT column (col. 9) so we can see what fields are there
    i=which(fields==fieldName) # This finds which item in the split apart list should correspond to the field name we want
    if (!(fieldName %in% fields)) stop('fieldName not found in format fields') # Check that the given field name is present
    return(sapply(x, `[[`, i)) # Finally, this uses that sapply trick to get the "ith" item from every vector in the list created by strsplit
}

### Now that I've got my function, my next step is to
### create a new table to hold the results:
AA=rep(0, nrow(my.data))
Aa=rep(0, nrow(my.data))
aa=rep(0, nrow(my.data))
NN=rep(0, nrow(my.data))
my.results=data.frame(AA, Aa, aa, NN)

### Next, I do a loop that goes through each line of the input vcf one at a time
for (i in 1:nrow(my.data)) {
    my.row = as.vector(my.data[i,], mode="character") # Get the ith row of the file, and make sure it is a vector object
    my.samples = samples=my.row[10:length(my.row)] # Get just the sample columns from the row (samples start at column 10)    
    my.gens = get.field(samples=my.samples, format=my.row[9], fieldName="GT") # Here is my get.field function in action (getting the GT field)
    my.results[i,] = table(factor(my.gens, levels=c("0/0","0/1","1/1","./."))) # Use the table command to automatically get the counts of each genotype
}

### Print the results
cat("Here are the Genotype Counts for Each Position:\n") # Cat is another way of printing something; it works well for text, but not very well for data objects
print(my.results)
cat("\n") # The \n symbol stands for "newline;" basically this prints out a blank line

#################################################################################
### Problem 2: Get the Percent Heterozygotes (Aa) for each INDIVIDUAL
#################################################################################

### To get the % heterozygotes for each sample, I am going to use my same get.field function,
### but this time run it on the columns instead of the rows

het.results = c() # I'm initializing an empty results vector that I will add to as I run my loop

for (i in 10:ncol(my.data)) { # Notice that here I start at 10 instead of at 1, because the first 9 columns of a VCF are NOT sample data
    my.col = as.vector(my.data[,i], mode="character") # Here I get the ith column (the same way I got the ith row before)
    format1 = my.data[1,9] # All of the formats should be the same for every site, so I just grab the first one
    my.gens = get.field(samples=my.col, format=format1, fieldName="GT") # Using my function again, but now I give it the column instead of the row
    hets = length(which(my.gens=="0/1")) # I only need to search for "0/1" to find heterozygotes
    total = (length(my.gens)) - (length(which(my.gens=="./.")))
    percent.het = (hets/total) * 100
    het.results=append(het.results, percent.het) # Add the result for the ith column to my results vector
}

names(het.results) = colnames(my.data)[10:ncol(my.data)] # This assigns the sample names to my results vector; I'll use these names to print the results

### Print the output:
cat("The parents samples (0% heterozygosity) are:\n\t")
print(het.results[which(het.results==0)]) # Here I select only the samples where the %het. is 0
cat("\nThe F1 samples (100% heterozygosity) are:\n\t")
print(het.results[which(het.results==100)])
cat("\nThe F2 samples (b/t 1% and 99% heterozygosity) are:\n\t")
print(het.results[which((het.results>0) & (het.results<100))])
cat("\n")

#################################################################################
### Problem 3: Convert the data to a table with nucleotides for each sample
#################################################################################

### Again I start by creating a new empty table to hold my results
num.samples=ncol(my.data) - 9 # I know my samples have to start at column 10, so I can figure out how many total samples there are like this
new.table = as.data.frame(matrix(NA, nrow=nrow(my.data), ncol=num.samples)) # Here is alternative way to make an empty table, starting with a matrix
colnames(new.table) = colnames(my.data)[10:ncol(my.data)] # I'm assigning the sample IDs as the column names of my results table
rownames(new.table) = my.data$ID # I'm assigning the SNP ID as the row name of my new table

homRef = paste0(my.data$REF, my.data$REF) # Use paste to create the diploid genotype strings for all homozygous reference cases
homAlt = paste0(my.data$ALT, my.data$ALT) # Do the same thing to make all of the possible homozygous ALT strings
het = paste0(my.data$REF, my.data$ALT) # And now make all of the heterozygous strings

for (i in 1:nrow(my.data)) {
    my.row = as.vector(my.data[i,], mode="character") # Same as above, get the row as a vector, then get the samples
    my.samples = samples=my.row[10:length(my.row)]
    new.samp = gsub("0/0.*", homRef[i], my.samples) # Use gsub to find a string starting with "0/0" followed by anything, and replace it with the homozygous Ref pattern
    new.samp = gsub("0/1.*", het[i], new.samp) # Do the same for the hets, BUT make sure to replace in new.samp (and not the old vector)
    new.samp = gsub("1/1.*", homAlt[i], new.samp)
    new.samp = gsub("\\./\\..*", "NN", new.samp) # Don't forget about escaping the special "." character for missing data!
    new.table[i,] = new.samp
}
cat("Here is the new table of diploid genotypes:\n")
print(new.table)
