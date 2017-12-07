### Read in the function source code
library(devtools)
source_gist("c954c757db02eb6c4ccfae1aa090658c", filename="compGen17_fxns.R", quiet=TRUE)

### Set the working directory,
### set stringsAsFactors=F
### Read in the sample data
options(stringsAsFactors=FALSE)
my.data=read.vcf("heliconius.vcf", header=TRUE)

### Find all of the window start positions
f=min(my.data$POS) # the "from" value will be the smallest position in the file
t=max(my.data$POS) # the "to" value will be the largest
win.start=seq(f,t,8000) # the "by" value = (window length - overlap) = (10000 - 2000 = 8000)
my.results=data.frame(Start=win.start, End=win.start+10000) # Now just add 10K to each start for the end pos.

### Find the number of base pairs in each window (amt. non-missing data)
my.results$NumBP = 0 # make a new column in results called NumBP
### In the next line, there are several things going on:
### the part that says c(x[1]:x[2]) is getting the full sequence of numbers between the window start (column 1) and the window end (column 2) position,
### y %in% then asks for every value of the y vector, which of them show up in that sequence of numbers I just made
### I define the y vector (values I want to test) as the POS column of my.data
### the y %in% x syntax by default just gets a True/False value for every item in y; which() forces it to only return the index of the ones where the they are "True", length tells me how many of those there are
### apply(my.results, 1 is tell R to automatically run this function on each row in my.results (1 is for row, if I used 2 it would operate on the columns)
### In a nutshell, this is a one line command that takes every start and end pos of each window, and gets the number of my.data positins that fall into it.
my.results$NumBP=apply(my.results, 1, function(x,y) length(which(y %in% c(x[1]:x[2]))), y=my.data$POS)

### Row 1775 has SNP, to see an example of variant vs. non-variant site
### Ways to filter: find sites where ALT is not "."
### Find sites where ALT count or maf is >0
### Find sites with more than one FORMAT tag (but won't always work-some SNP files will only have GT)
my.data[c(1,1775),1:13]

### I'm going to find rows where the ALT column just has "."
### I use the "-" notation inside of the brackets to get rid of them
my.data=my.data[-grep("^\\.$", my.data$ALT),]

### Now, I use the same fancy apply strategy as before to get the SNP counts
my.results$NumSNP=apply(my.results, 1, function(x,y) length(which(y %in% c(x[1]:x[2]))), y=my.data$POS)

### Set up 3 columns in results to hold TajimasD, Theta, and Pi
my.results$Tajima=0
my.results$Theta=0
my.results$Pi=0

### We don't have to worry about it in this data set,
### but I usually also make sure to set windows with 0 BP
### (i.e. no data) to NA, to distinguish them from actual 0 values
my.results$Tajima[which(my.results$NumBP==0)]=NA
my.results$Pi[which(my.results$NumBP==0)]=NA
my.results$Theta[which(my.results$NumBP==0)]=NA

### Loop through the windows, and for each one,
### calculate the 3 statistics
for (i in 1:nrow(my.results)) {
    w=my.data[which(my.data$POS %in% my.results$Start[i]:my.results$End[i]),] # get the data in the window   
    theta=waterson.theta(w, perBP=FALSE) # use my function for theta-w, set perBP = F (I'll scale by my own BP count, since i might have missing data)    
    pi=pi.diversity(w, perBP=FALSE) # same with pi as with theta-w
    n=2*(ncol(w)-9) # get the value for 2N
    s=nrow(w) # get the # snps (I could have also pulled this from my.results)
    var=variance.d(n,s) # get the variance of d
    my.results$Tajima[i]=(pi-theta)/sqrt(var) # calculate and store tajima's D
    my.results$Pi[i]=pi/my.results$NumBP[i]
    my.results$Theta[i]=theta/my.results$NumBP[i]
}

### Calculate Nucleotide Divergence (Dxy)
### First, split the data by population
aglaope=grep("^Ag", colnames(my.data))
amaryllis=grep("^Am", colnames(my.data))

### Set up an empty column for Dxy
my.results$Dxy=0
my.results$Dxy[which(my.results$NumBP==0)]=NA
for (i in 1:nrow(my.results)) {
    w=my.data[which(my.data$POS %in% my.results$Start[i]:my.results$End[i]),]
    vcf1=w[,-amaryllis]
    vcf2=w[,-aglaope]
    d = dxy(vcf1,vcf2, perBP=FALSE)
    my.results$Dxy[i]=d/my.results$NumBP[i]
}
###quantile(my.results$Dxy, na.rm=TRUE)

### Set the margins to allow room on the right side for another axis
par(mar=c(5,5,2,5))

### Create an empty plot
plot(my.results$Start, my.results$Theta, type="n",
     xlab="Position",
     ylab=expression(paste("Population Mutation Rate (", theta[W], " & ", pi, ")")), ylim=c(-0.054,0.07))

### Plot Theta as a filled in shape
polygon(c(my.results$Start[1], my.results$Start, my.results$Start[nrow(my.results)]),
        c(0, my.results$Theta, 0), col=rgb(0.545,0,0, alpha=0.99), border=NA)

### Plot Pi as another filled in shape,
### but make it semi-transparent so we can see the Theta
### values behind it
polygon(c(my.results$Start[1], my.results$Start, my.results$Start[nrow(my.results)]),
        c(0, my.results$Pi, 0), col=rgb(0.933,0.705,0.133, alpha=0.6), border=NA)

### Add a line for Dxy
lines(my.results$Start, my.results$Dxy, col="black")

### To add Tajima's D, which is scaled differently, I need to make a new plot on top of this one,
### and add a separate set of axes
par(new=T)
plot(my.results$Start, my.results$Tajima, type="l",
     col="black", ylim=c(-1.5,2), lty=3, axes=F, xlab=NA, ylab=NA)
axis(side=4)
mtext(side=4, line=3, "Tajima's D")

### Finally, let's make a legend
lgd <- legend("topright", legend=c(expression(theta[W]), expression(Pi), expression(D[XY]), "Tajima's D"),
              lty=c(NA,NA,1,3),
              col=c(NA,NA,"black","black"), bty="n")
legend(lgd$rect$left, lgd$rect$top, legend=c(NA,NA), pch=c(22,22),
       fill=c(rgb(0.545,0,0,alpha=0.99), rgb(0.933,0.705,0.133,alpha=0.6)), border=NA, bty="n", col=NA)
       
