# Computational Genomics Workshop:
## R coding for the analysis of SNP data and applied population genetics

### Description:
High-throughput sequencing technology has made it possible to obtain large scale genetic data sets for almost any organism, creating a need for computational tools and skill sets to process these data.  While the bioinformatics workflows for processing raw data into SNPs are typically well delineated, the path for analyzing and interpreting the resulting SNP data set can be less clear.  Classical population genetics theory is a rich field with many existing statistics and mathematical models, all created to draw biological inferences from nothing more than the patterns of observed gene and allele frequencies.  Since these frequencies can be calculated from any SNP data set, this makes them highly applicable to current research.  In this workshop, students will learn about statistics that test the neutral theory, and then apply the test to a sample SNP data set.  Throughout this process, students will learn to write their own R code for each analysis, and emphasis will be placed on algorithm design with scalability for large data sets in mind.  

### Goals:
1.	Achieve a thorough understanding of each statistic or model, especially in terms of neutral theory expectations, and interpretations of deviations from the neutral expectation.  Emphasis will also be placed on knowing when a certain model may or may not be appropriate for a particular data set or scenario.
2.	Make the connection between the table of genotypes obtained in a VCF file (standard next-gen pipeline output) and a theoretical equation that comes from a textbook or paper.  By writing their own code each week, students will get a feel for the best ways to set up algorithms by dividing tasks into re-usable functions, parsing difficult file structures, and creating useful output.  
3.	Learn (or improve) R programming skills with hands-on practice.  For each exercise, we will focus on writing code from scratch (rather than using pre-existing packages), so that students will become familiar with common coding structures such as functions and loops.  These skills are applicable to all coding languages, rather than being R specific.

### Prerequisites:
+	GEN 3000 or an equivalent undergrad genetics course from another institution.  Students must have a working knowledge of the basic concepts of genetics, including Mendelian inheritance and the molecular basis of heredity.  Prior knowledge of population genetics is NOT required.
+	ALL LEVELS OF PROGRAMMING EXPERIENCE ARE WELCOME! This course is designed for both biologists who want to learn bioinformatics, and for people with programming experience who want to learn more about biology or population genomics.
+	In previous semesters, some students who were completely new to R found it helpful to also go through an online tutorial of R to get more comfortable with it in the beginning.  They recommended ![this one](http://nathanieldphillips.com/thepiratesguidetor/) in particular.

### General Structure:
+  1-3 credits; Grades assigned based on participation (70%) and final project (30%)
+  Meet 1 time per week, for 2.5 hours
	++  First 30-45 minutes: Review results or questions from the previous week.  Then, discuss one key test statistic/concept; go over in detail the meaning, derivation, assumptions, and expectations of the central model or statistic for the week.
	++  1.5 hours: Briefly discuss algorithm for implementing the model in a VCF file, then everyone will work on their code while interacting with the instructor and the other students for assistance and feedback*. 

*The exact amount of time devoted to lecture vs. in-class programming will change from week to week, depending on what the topic is.  Additionally, some weeks the in-class exercise may not take the full amount of time, while other exercises might span into the next week of class.*

### Meeting Time:
TUESDAYS 2:00PM – 4:30PM

### Schedule
#### 1: Introduction to Course and R basics
+  Review the syllabus and course structure.
+  Brief introductions.
+  R-tutorial: work through creating and manipulating variables, common data structures, performing built-in functions, and writing loops.
+  Getting help: finding and reading R documentation.
+  *Exercises: Using R to perform some simple tasks, including creating a table and writing a small loop to perform some calculation on each row.*

### 2: The Variant Call Format (VCF)
+  Discuss the major concepts involved in next-gen sequencing, and review the standard raw data pipelines that are used to call SNPs and produce VCF files.
+  Look in detail at the VCF format, and learn how to read this format into R and manipulate the resulting table.
+  *Exercise: Read a sample VCF file into R, and perform some simple functions to extract, summarize, and re-code the SNP data.*

**Single-Site Tests**

### 3: Hardy-Weinberg Equilibrium
+  Review principle of HWE; its usefulness as a diagnostic test in next-gen data; calculating observed and expected frequencies; thinking about when you do and don’t expect deviations.
+  Use Fisher’s Exact test to find statistically significant deviation.
+  *Exercise: Write a script to calculate observed and expected genotype frequencies at each site in a sample VCF file, then extract sites that show too much deviation.*

### 4: Wright's FST 
+  Discussion of population structure, what causes it, and how FST is used to measure it.
+  Equations and assumptions for FST; effects of population size on genetic drift; relationship between FST and migration; review of alternative FST estimators (and when to use them).
+  Exercise: write code to calculate single site FST (using Wright’s original equation based on single site frequencies) between 2 populations in a sample VCF file; plot the distribution of FST.

**Pairwise Tests**

### 5: Linkage Disequilibrium pt. 1 
+  Definition of LD and underlying causes.  Estimators of LD, and using haplotypes versus inferred haplotypes.
+  *Exercise: write a script to find 2-locus haplotype frequencies (in addition to single site allele frequency), and calculate the LD estimators D and r2 between every PAIR of SNPs in a very small (test case) VCF file.  Use previous functions for genotype and allele frequencies, but create a new nested loop structure that will get site pairs.  Save functions and code to run on a larger file next week.*

### 6: Linkage Disequilibrium pt. 2
+  Discussion of final project guidelines (in time to review papers and pick data sets over fall break)
+  The interpretation and meaning of LD decay.  How to estimate the recombination rate (rho) based on the mathematical relationship between LD and recombination.
+  *Exercise: Using the code from last week, calculate r2 in a full 1000 row VCF file and plot the decay of LD with distance.  Fit a non-linear model to estimate recombination from the rate of decay.*


### Fall Break

**Whole Gene Tests**

### 7: McDonald-Kreitman Test
+  Review what positive natural selection is, and the theory underlying the MK test to find sites under selection.
+  Discussion of gene annotation: what this means, and what kind of programs give you this information.
+  *Exercise: using an already annotated sample VCF file, perform the MK test to find sites under selection.*

### 8: Allele Frequency Spectrum
+  A brief introduction to coalescent theory.
+  Neutral (coalescent) theory expectations of allele frequency distributions; selective and demographic forces causing deviations from neutral; interpretation of single gene vs. genome-wide spectra.
+  *Exercise: write a loop to calculate and plot the derived allele frequency spectrum for individuals in a VCF file.  Use the functions for allele frequencies from previous weeks.  Use R to simulate a neutral expectation density curve, and overlay on the observed spectra for comparison.*	

### 9: Population Mutation Rate (Theta)
+  Define Theta and its estimators; discuss Watterson’s theta vs. Pi; how to calculate each; the different expectations of each, and how this relates back to the Allele Frequency Spectrum.
+  *Exercise: For a VCF file with SNPs from a single locus; calculate both estimators of theta, and determine the per bp mutation rate based on the number of mapped sites.  Both estimators will require some new functions for the calculations; the Pi estimator in particular will require a new strategy of comparing columns in a VCF rather than rows; Watterson’s theta should be pretty straightforward.*

** Sliding Window Tests**

### 10: Tajima’s D 
+  Testing for selection using the difference in estimates of theta; how to calculate D; different types of selection predicted by negative vs. positive D.  Discuss extensions of D (e.g. Fay and Wu’s H).
+  *Exercise: Perform a sliding window analysis on a VCF file.  Within each window, calculate D (based on difference between theta-w and pi).  Make a genome-wide plot of the D values to find regions under selection.*
+  **Deadline to Select Topic for Final Project!**


### 11: Multi-site FST and Sliding Window Tests with Missing Data 
+  Finding, installing, and using R packages from external sources, R packages for VCF data, and dealing with large files.
+  Why and how do you get gaps or missing data in next-gen sequencing?  What are the implications of missing data in terms of 1) What tests you should perform on your data, and 2) what programs you should run to get the most information possible?
+  Combining what we learned about  and FST in previous weeks to calculate the multi-site estimator ST.  Discuss the limitations of this statistic, and go over a paper (Cruickshank and Hahn 2014) that demonstrates these limitations and an alternative DXY statistic.
+  *Exercise: perform a sliding window test to compare the FST and DXY statistics on a sample data set.*

** Projects**

### Week 12: Running R on the cluster 
+  Brief overview of the linux command line environment (navigating directories, creating, deleting, splitting, concatenating, and moving files).  
+  Cluster basics: logging in, transferring files, using interactive nodes and submitting job scripts.
+  Running R from the command line; writing output for data and figures to files; making code more flexible with command line arguments. Installing packages locally.
+	*In-class time to work on projects and get feedback as needed*

### Last Day of Class: Project Presentations
+  10 minute presentations by each person on the results of their project.
