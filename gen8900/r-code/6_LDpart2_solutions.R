### Set my working directory
setwd("~/Desktop/Teaching/CompGenWS-GEN8900/")

### Read in the my functions file
library(devtools)
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)

### Read in the VCF file
my.data=read.vcf("SampleData/6_sampleData_LD.vcf", header=TRUE, stringsAsFactors=FALSE)

### Set up the empty results table

### To calculate the number of results
### I will have, I first use R's built-in
### combinatorics function
### the code below calculates (x choose 2),
### which is equivalent to x!/(2! * (x-2)!)
num.results = ncol(combn(x=nrow(my.data), m=2))

### Now, I can set up my empty table the same
### way I usually do:
my.results = data.frame(Distance=rep(0, num.results), Rsquared=rep(0, num.results))

### Set up a counter to help keep track
### of what row of the results I am on:
counter = 1

### Use a nested loop to calculate R-squared
### and distance for every SNP pair.
### Inside of this loop, you will see my same
### commands for getting pA, pB, pAB, D, and rsquared
### Then, whenever I get a result, I assign it to whatever
### row the "counter" is pointing to in my results table,
### and I advance the counter by 1
for (i in 1:(nrow(my.data)-1)) {
    row1=as.vector(my.data[i,], mode="character") # Instead of getting the 1st row, I get the ith row
    g1 = get.field(row1[10:length(row1)], row1[9], "GT") # This line and the next are copied exactly from my previous weeks solutions
    pA = unname(allele.freq(count.genotypes(g1))["p"])
    for (j in (i+1):nrow(my.data)) {
        row2=as.vector(my.data[j,], mode="character") # Here, instead of getting the 2nd row, I get the jth row
        my.results$Distance[counter] = as.numeric(row2[2]) - as.numeric(row1[2]) # Calculate the distance between the SNPs (copied exactly from my solution from last week)
        g2 = get.field(row2[10:length(row2)], row2[9], "GT") # Copied exactly from last week's solution
        pB = unname(allele.freq(count.genotypes(g2))["p"]) # Copied exactly from last week's solution
        h=get.haplotypes(g1,g2) # Copied exactly from my pAB method II solution
        pAB = (length(h[h=="00"]))/length(h) # Copied exactly from my pAB method II solution
        D = pAB - (pA * pB) # Copied exactly from last week
        my.results$Rsquared[counter] = (D**2)/(pA*(1-pA)*pB*(1-pB)) # The only change here is I assign my R-squared value to the results table; the equation for getting R-squared is copied from last week
        counter=counter+1 # Advance the counter by 1, so that on the next generation I'll put the results in the next row
    } # Close the j loop (the inner loop)
} # Close the i loop (the outer loop)

### Plot the Decay of LD with Distance
plot(my.results$Distance, y=my.results$Rsquared, xlab="Distance between SNPs (bp)", ylab=expression(R^2), pch=20, col=rgb(0,0,0,alpha=0.2), main="Decay of Linkage Disequilibrium")
        
### Add a smoothed moving average line
bins=seq(1,max(my.results$Distance), by=500) # Get all numbers between 1 and the maximum Distance value, in increments of 500
my.means=rep(0, length(bins)) # Set up an empty vector to hold the mean r-squared for each "bin"
LD.averages=data.frame(bins, my.means) # Create the results table (to hold calculations)

### Go through all of the "bin" values,
### For each one, find the subset of data that
### falls in that bin, and get the mean r-squared
for (i in 1:length(bins)) {
    data.interval=subset(my.results, (my.results$Distance >= bins[i] & my.results$Distance < (bins[i]+500))) # Get the subset of data that falls in the ith bin
    LD.averages$my.means[i]=mean(data.interval$Rsquared) # Calculate the average R-squared value for this data subset, and save in my "LD.averages" results table
}
### Add points and a line for my moving average on top of my exisiting plot
points(x=LD.averages$bins, y=LD.averages$my.means, col="red", pch=20) 
lines(x=LD.averages$bins, y=LD.averages$my.means, col="red", lwd=2)

### Find the point when LD drops below 0.1
### The which statement tells me which rows of the table have means
### greater than 0.1,
### the brackets are giving me the subset of the "bins" with averages
### that correspond to these rows, and the max statement tells me which bin has
### the greatest distance value
LD.drop = max(LD.averages$bins[which(LD.averages$my.means>0.1)])
abline(v=LD.drop, col="blue", lty=2) # Add a vertical line corresponding to the drop off point
text(25500, 0.8, "LD drops below 0.1", col="blue") # label that line

### Find the LD half life
LD.half = (max(LD.averages$my.means))/2 # Calculate half of the maximum average LD value
LD.half.point = max(LD.averages$bins[which(LD.averages$my.means>LD.half)]) # Same strategy as when I found the 0.1 drop off point
abline(v=LD.half.point, col="darkgreen", lty=2) # add a vertical line at the half life point
text(8500, 0.9, "LD half life", col="darkgreen") # label this line

### Calculate rho based on Hill Weir equation:
n=length(g1) # get the number of samples by finding the number of genotypes
hill.weir.eq = (Rsquared~(((10+(rho*Distance))/((2+(rho*Distance))*(11+(rho*Distance))))*(1+(((3+(rho*Distance))*(12+(12*(rho*Distance))+((rho*Distance)**2)))/(n*(2+(rho*Distance))*(11+(rho*Distance))))))) # Define the formula for the Hill-Weir equation

rho.start=0.1 # Set an arbitrary starting value for rho
m=nls(formula=hill.weir.eq, data=my.results, start=list(rho=rho.start)) # Run the function to test different rho values and find the best fit to the data
results.m=summary(m) # Get a summary of results from the regression model

### Plot Expected R-squared
rho.estimate=results.m$parameters[1] # extract the rho estimate that produced the best fit
Distance=sort(my.results$Distance) # sort the distance results from lowest to highest

exp.rsquared=(((10+(rho.estimate*Distance))/((2+(rho.estimate*Distance))*(11+(rho.estimate*Distance))))*(1+(((3+(rho.estimate*Distance))*(12+(12*(rho.estimate*Distance))+((rho.estimate*Distance)**2)))/(n*(2+(rho.estimate*Distance))*(11+(rho.estimate*Distance)))))) # Use the best rho estimate inside of the Hill-Weir formula to calculate the expected r-squared at each distance value

lines(Distance, exp.rsquared, col="purple", lwd=2) # plot a line for the expected values (the same as a best fit line)
legend(32000,0.99, c("Means", "Expected R2"), lty=c(1,1), col=c("red", "purple")) # add a legend
