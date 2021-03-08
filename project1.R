#####Install####
##if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install(version = "3.12")
###install packages######
#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install()
####

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install(version = "3.12")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("dada2")  

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ShortRead") 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ggplot2") 

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("phyloseq") 

#if (!requireNamespace("BiocManager", quietly = TRUE))
#install.packages("BiocManager")
#BiocManager::install("cutadapt") 

library(dada2); #packageVersion("dada2"); citation("dada2")
library(ShortRead); #packageVersion("ShortRead")
library(ggplot2); #packageVersion("ggplot2")
library(phyloseq); #packageVersion("phyloseq")

packageVersion("dada2")
packageVersion("phyloseq")

getwd()
#setwd("C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/fastqfiles")

#setting path to directory containing our fasstq files

path <- "/Users/gracebeery/Desktop/BI586/fishmicrobes_2021_assign1" 
path <- "C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/fastqfiles"
path <- "/usr4/bi594/vfrench3/assignment1/fishmicrobes_2021_assign1/fastqfiles"
path <- "/Users/Victoria1/Desktop/Grad School/Eco Gen./fishmicrobes_2021_assignment1/fastqfiles"
fns <- list.files(path)
fns


#read in names of fastq files, and divide them into matched lists of forward and reverse fastq files using _1 and _2 filename endings
fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures reads are in same order
fnFs <- fastqs[grepl("_1", fastqs)] #assigning forward reads to variable 
fnRs <- fastqs[grepl("_2", fastqs)] #assigning reverse reads to variable 


#removing .fastq on each on the file names
# these caused an error because we do not need specific sample.names for both forward and backward, because they're all the same set of samples
# sample.names.F <- sapply(strsplit(fnFs, ".fastq"), `[`, 1) #the last number will select the field for renaming
# sample.names.R <- sapply(strsplit(fnRs, ".fastq"), `[`, 1) #the last number will select the field for renaming
# sample.names.F
# sample.names.R
sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names



# Specify the full path to each fastq file
fnFs <- file.path(path, fnFs) 
fnRs <- file.path(path, fnRs)
fnFs
fnRs

#here we are visualizing the quality profiles of the forward and reverse reads
plotQualityProfile(fnFs[1:9])
plotQualityProfile(fnFs[10:18])

plotQualityProfile(fnRs[1:9])
plotQualityProfile(fnRs[10:18])

#placing trimmed and filtered files in a subdirectory
filt_path <- file.path(path, "trimmed") #assigning new trimmed folder to a variable 
if(!file_test("-d", filt_path)) dir.create(filt_path)  #creating trimmed directory 
filtFs <- file.path(path, "filtered", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "filtered", paste0(sample.names, "_R_filt.fastq.gz"))
names(filtFs) <- sample.names
names(filtRs) <- sample.names

#setting standard filtering parameters for maxN, truncQ=2, rm.phix=TRUE 
#we set trunclen to be 240 because no sample falls below a quality score of 20 after 250 bp. 
#After trial and error data manipulation, 240 maximizes read output during initial filtering and minimizes presence of bimeras. While still facilitating successful paired end merges. 
#Trimleft excluded as raw sequences from NCBI contain no primer sequences. Primer Sequences used were: 
#515F ('GTGYCAGCMGCCGCGGTAA') and 806R ('GGACTACNVGGGTWTCTAAT') for the 16S V4 region according to earth microbiome project. 

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=240,
                     maxN=0, 
                     maxEE=1, 
                     truncQ=2, 
                     rm.phix=TRUE,
                     compress=TRUE, multithread=FALSE) 

head(out)
tail(out)

#here DADA2 is using parametric error model to determine error rate of this amplicon dataset
setDadaOpt(MAX_CONSIST=30) #usually keep 30, increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

#estimated error rates (black lines) are good fit to the observed rates (points) and error rates drop with increased quality
plotErrors(errF, nominalQ=TRUE) #plotting error rates, want to see mostly linear relationships
plotErrors(errR, nominalQ=TRUE)

#dereplicating reads - collapsing all reads into unique sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <-sample.names
names(derepRs) <- sample.names
#names(derepFs) <- sample.names.F
#names(derepRs) <- sample.names.R


#infer sequence variants 
#Must change some of the DADA options b/c original program optomized for ribosomal data, not ITS - from github, "We currently recommend BAND_SIZE=32 for ITS data." leave as default for 16S/18S
setDadaOpt(BAND_SIZE=16) #Default is 16 for Illumina (low rates of indels) according to documentation 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)



#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[1]]
dadaFs[[5]]
dadaRs[[1]]

#merge paired ends - here we are merging the forward and reverse reads together to get full denoised sequences
#we merge by aligning denoised forward reads with reverse-complement of corresponding denoised reverse reads. we then make merged contig sequences
#Set minOverlap to ___ because of the length 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

#construct sequence table - now making amplicon sequence variant (ASV) table which is a higher res version of OTU table 
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
head(mergers)

#remove chimeras - chimeric sequences are identified and removed
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)
#identified __ bimeras out of __ input sequences 

sum(seqtab.nochim)/sum(seqtab)
#chimeras account for only about 1% of the merged sequence reads

write.csv(seqtab,file="fishmicrobes_seqtab.csv")
write.csv(seqtab.nochim,file="fishmicrobes_nochim.csv")

################################
##### Track Read Stats #######
################################

#here we are looking at number of reads that made it through each step in pipeline
getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names 
head(track)
tail(track)

write.csv(track,file="ReadFilterStats_AllData_final.csv",row.names=TRUE,quote=FALSE)
#this shows raw reads, how many filtered, denoised, etc.
#final number of reads under nonchim column 
#could include table of this in paper 

################################
##### Assign Taxonomy #######
################################


#assigning ASVs in seqtab.nochim to taxonomy via Silva v132 (most recent version); Silva primary database for 16S; arguments set at preset

#write taxa to csv file 
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/silva_species_assignment_v138.fa.gz")

taxa <- assignTaxonomy(seqtab.nochim, "/users/Victoria1/Desktop/Grad School/Eco Gen./fishmicrobes_2021_assignment1/fastqfiles/silva_nr_v132_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "/users/Victoria1/Desktop/Grad School/Eco Gen./fishmicrobes_2021_assignment1/fastqfiles/silva_species_assignment_v132.fa.gz")

taxa.print <- taxa
rownames(taxa.print) <- NULL #this gets rid of all the sequences
head(taxa) #prints full info, including sequence
head(taxa.print) #this just prints taxonomy


#Obtain a csv file for the taxonomy so that it's easier to map the sequences for the heatmap.

write.csv(taxa, file="taxa.csv",row.name=TRUE,quote=FALSE)
unname(head(taxa, 30))
unname(taxa)


#Now, save outputs so can come back to the analysis stage at a later point if desired
saveRDS(seqtab.nochim, file="final_seqtab_nochim.rds")
saveRDS(taxa, file="final_taxa_blastCorrected.rds")

#If you need to read in previously saved datafiles
seqtab.nochim <- readRDS("final_seqtab_nochim.rds")
taxa <- readRDS("final_taxa_blastCorrected.rds")
head(taxa)

################################
##### handoff 2 phyloseq #######
################################

#import dataframe holding sample information
#have your samples in the same order as the seqtab file in the rows, variables as columns
#before doing phyloseq, make sure sample names in the variable table have the same names as the sample names in the taxa and seqtab.nochim files
samdf<-read.csv("variabletable.csv")

head(samdf)
head(otu_table(seqtab.nochim, FALSE))
head(taxa)
rownames(samdf) <- samdf$Sample
head(samdf)


# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE),
               sample_data(samdf),
               tax_table(taxa))
ps


#replace sequences with shorter names (correspondence table output below)
ids<-taxa_names(ps)
ids <- paste0("sq",seq(1, length(colnames(seqtab.nochim))))
colnames(seqtab.nochim) <- ids

#Bar-plots
top90 <- names(sort(taxa_sums(ps), decreasing=TRUE))[1:90]
ps.top90 <- transform_sample_counts(ps, function(OTU) OTU/sum(OTU))
ps.top90 <- prune_taxa(top90, ps.top90)

plot_bar(ps.top90, x="Sample", fill="Class") 

#visusalize via counts rather than abundances:
plot_bar(ps, x = "sample", fill= "Class") #+ facet_wrap("tank")
#
#Obtain a csv file for the phyloseq data. Will give you the abundances for each sample and class. Useful for constructing the heatmap. Also, enables you to use ggplot, and construct more aesthetically pleasing plot.
psz <- psmelt(ps.top90)
write.csv(psz, file="Phyloseqoutputfinal.csv")
p <- ggplot(psz, aes(x=Sample, y=Abundance, fill=Class))
p + geom_bar(stat="identity", colour="black")
