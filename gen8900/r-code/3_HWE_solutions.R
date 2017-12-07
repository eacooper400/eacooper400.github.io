### Set my working directory to the folder
### where I keep the sample data files for class
setwd("~/Desktop/Teaching/CompGenWS-GEN8900/")

### Get the most up to date version of
### the R functions script
library(devtools)
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)

### Use the read.vcf function to read in
### the sample data file for week 3:
my.data=read.vcf(file="SampleData/sampleData_wk3.vcf", header=TRUE, stringsAsFactors=FALSE)

### EXERCISE 1: CALCULATE EXPECTED GENOTYPE COUNTS

### Since the first part of this exercise
### is the same as one of our exercises from
### last week, I'm mostly going to re-use code from
### the week 2 solutions
### But, instead of just copying and pasting,
### I'm going to use that code to make some new functions.

### I've already got a function to get the genotypes,
### now I'm adding a function to count them
count.genotypes <- function(genotypes) {
    genotypes = gsub("(\\||/)", "", genotypes) # First, I'm getting rid of the / or | symbol, so I can be sure the genotypes look the same no matter what
    gen.patterns = c("00", "01", "10", "11", "..") # Here are all the possible patterns to count
    my.counts=table(factor(genotypes, levels=gen.patterns)) # Then, I use the same code I had in my solution last week
    final.counts = c(my.counts[1], (my.counts[2]+my.counts[3]), my.counts[4:5]) # I need to sum together the 2 heterozygote options
    names(final.counts) = c("AA", "Aa", "aa", "NN") # Finally, I give some names to my vector to make it easy to find a specific count
    return(final.counts)
}

### The next function I'm going to create
### will calculate allele frequencies
### This function will take the genotype
### counts as input, and return p and q
allele.freq <- function(genotypeCounts) {
    n = sum(genotypeCounts) - genotypeCounts["NN"]
    p = ((2*genotypeCounts["AA"]) + genotypeCounts["Aa"])/(2*n)
    q = 1-p
    freqs = c(p,q)
    names(freqs) = c("p", "q")
    return(freqs)
}

### With my new functions all set up,
### I'm now going to do all of Exercise 1
### in one loop

### Set up an empty results table
blank.col = rep(0, nrow(my.data))
my.results=data.frame(obsAA=blank.col, obsAa=blank.col, obsaa=blank.col, p=blank.col, q=blank.col, N=blank.col, expAA=blank.col, expAa=blank.col, expaa=blank.col)


### Run the loop to do all calculations
### on each row
for (i in 1:nrow(my.data)) {
    my.row = as.vector(my.data[i,], mode="character")
    my.samples = samples=my.row[10:length(my.row)]
    my.gens = get.field(samples=my.samples, format=my.row[9], fieldName="GT")
    obs.counts = count.genotypes(my.gens) # Using my new function to get the counts
    my.results[i, 1:3] = c(obs.counts["AA"], obs.counts["Aa"], obs.counts["aa"]) # Fill in the observed counts (the first 3 columns of the my.results table) for the ith row
    obs.freq = allele.freq(obs.counts) # Use the other new function to get frequency
    my.results[i,4:5] = obs.freq # Fill in the next 2 columns (4 and 5) with the p and q values for the ith row
    my.results$N[i] = (sum(obs.counts)) - obs.counts["NN"] # Get the total of the Non-Missing values
    my.results$expAA[i] = (my.results$p[i]**2) * my.results$N[i] # Calculate expAA as p-squared times N (for the ith row)
    my.results$expAa[i] = 2 * my.results$p[i] * my.results$q[i] * my.results$N[i] # expAa = 2pq times N
    my.results$expaa[i] = (my.results$q[i]**2) * my.results$N[i] # expaa = q-squared times N    
}

### EXERCISE 2: Use Fisher's Exact Test to see which sites differ significantly from HWE

### First, I create an empty Pvalue column and add it to my results table:
Pvalue = blank.col
my.results = cbind(my.results, Pvalue)

### Next, I run a for() loop on my results table and use the
### code given along with the exercise to run the Fishers test
for (i in 1:nrow(my.results)) {
    a=as.integer(my.results$obsAA[i])
    b=as.integer(my.results$obsAa[i]/2)
    c=as.integer(my.results$obsAa[i] - b)
    d=as.integer(my.results$obsaa[i])
    x=matrix(c(a,b,c,d), nrow=2, ncol=2, byrow=TRUE)
    test=fisher.test(x)
    my.results$Pvalue[i]=test$p.value
}

### EXERCISE 3: Check Loci with Significant Deviation from HWE

### Part A: How many loci deviate from HWE?
sig.loci = which(my.results$Pvalue < 0.05) # This returns which row numbers have a pvalue less than 0.05
num.sig.loci = length(sig.loci) # To get the count I just take the length (you should get 44 loci)

### Part B:  What might cause these loci to deviate?

### As shown in the exercises, I can figure out how
### loci deviate by calculating the difference between
### observed and expected heterozygosity
my.results$Diff=my.results$obsAa-my.results$expAa

### To check if samples where obsAa is too low
### have low read depth, I use the get.fields
### function from week 2 solutions to get the read
### depth for every sample
my.results$rdDepth=rep(0,nrow(my.results))

for (i in 1:nrow(my.data)) {
    d=get.field(as.vector(my.data[i,10:ncol(my.data)], mode="character"), my.data[i,9], "DP")
    d[d=="."]=NA # Get rid of missing data characters so that we can treat the vector as numeric
    my.results$rdDepth[i]=mean(as.numeric(d), na.rm=TRUE)
}

### Now, I'll use subset() to look at loci
### where: 1) pvalue < 0.05, 2) Diff < 0 (obs Aa too low)
sig.results = my.results[sig.loci,] # this uses the results from which() that I got earlier to get the counts, and returns a table with only those rows
too.low = subset(sig.results, sig.results$Diff<0)

### When you print the "too.low" table, you should see 5 rows, and avg.
### read depth for all of them is around 1.5, so definitely too low!

### Next, let's look at samples where obs(Aa) is too high
### Notice that I can subset the original my.data table using
### a condition that I test in a different table (the results)
### this is because their rows line up (i.e. row 1 of my.results
### corresponds to a calculation done on row 1 of my.data, and so on...)
too.high = my.data[which((my.results$Pvalue<0.05) & (my.results$Diff>0)),] # I use which to test 2 conditions, and get only rows where both conditions are met

### To look at the results more easily, I'll only print the first
### 2 columns, since that is where the position info. is
too.high[,1:2]

### Notice that the rows at the end (row # 990-999) are
### 10 positions all right next to each other - this definitely
### jumps out as a possible paralog.  The rest of the sites are less
### obvious.  A few sections (15-18, 415-418) are a bit suspect, but
### the others don't jump out as mapping errors
