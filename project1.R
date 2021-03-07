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

getwd()
#setwd("C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/fastqfiles")

#path <- "/Users/gracebeery/Desktop/BI586/fishmicrobes_2021_assign1" 
path <- "C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/fastqfiles"
path <- "/usr4/bi594/vfrench3/assignment1/fishmicrobes_2021_assign1/fastqfiles"
fns <- list.files(path)
#Let's make sure that all of our files are there
fns

fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures reads are in same order
fnFs <- fastqs[grepl("_1", fastqs)]
fnRs <- fastqs[grepl("_2", fastqs)]

sample.names.F <- sapply(strsplit(fnFs, ".fastq"), `[`, 1) #the last number will select the field for renaming
sample.names.R <- sapply(strsplit(fnRs, ".fastq"), `[`, 1) #the last number will select the field for renaming
sample.names.F
sample.names.R

# Specify the full path to each fastq file
fnFs <- file.path(path, fnFs)
fnRs <- file.path(path, fnRs)
fnFs
fnRs

plotQualityProfile(fnFs[1:9])
plotQualityProfile(fnFs[10:18])
plotQualityProfile(fnRs[1:9])
plotQualityProfile(fnRs[10:18])

filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, "filteredF", paste0(sample.names.F, "_F_filt.fastq.gz"))
filtRs <- file.path(filt_path, "filteredR", paste0(sample.names.R, "_R_filt.fastq.gz"))

out<- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=250, #set to 250 currently until Sarah responds. I think 250 works for most of the data but their are those weird outliers. 
                     maxN=0, #DADA does not allow Ns
                     maxEE=1, #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     trimLeft=19, #N nucleotides to remove from the start of each read: 16S (microbial community barcoding region) V4 region (used in paper even though V2-V3 higher resolution and species level identification; https://doi.org/10.1038/sdata.2019.7) forward primer 515 F = 19
                     rm.phix=TRUE, #remove reads matching phiX genome
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

out <- filterAndTrim(fnFs, filtFs, fnRs, filtRs, truncLen=250,
                     maxN=0, maxEE=1, truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)

head(out)
tail(out)

setDadaOpt(MAX_CONSIST=30) #usually keep 30, increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)
errR <- learnErrors(filtRs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE) #plotting error rates, want to see mostly linear relationships
plotErrors(errR, nominalQ=TRUE)

#dereplicating reads - collapsing all reads into unique sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
derepRs <- derepFastq(filtRs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names.F
names(derepRs) <- sample.names.R

#infer sequence variants 
#Must change some of the DADA options b/c original program optomized for ribosomal data, not ITS - from github, "We currently recommend BAND_SIZE=32 for ITS data." leave as default for 16S/18S
setDadaOpt(BAND_SIZE=16) #Default is 16 for Illumina (low rates of indels) according to documentation 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)
dadaRs <- dada(derepRs, err=errR, multithread=TRUE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[1]]
dadaRs[[1]]

#merge paired ends 
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose=TRUE)
head(mergers[[1]])

#construct sequence table
seqtab <- makeSequenceTable(mergers)
head(mergers)

#remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

write.csv(seqtab,file="fishmicrobes_seqtab.csv")
write.csv(seqtab.nochim,file="fishmicrobes_nochim.csv")

################################
##### Track Read Stats #######
################################

getN <- function(x) sum(getUniques(x))
track <- cbind(out, sapply(dadaFs, getN), sapply(dadaFs, getN), rowSums(seqtab), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoised", "merged", "tabled", "nonchim")
rownames(track) <- sample.names #getting error!!!!!!!!!!!
head(track)
tail(track)

write.csv(track,file="ReadFilterStats_AllData_final.csv",row.names=TRUE,quote=FALSE)
#this shows raw reads, how many filtered, denoised, etc.
#could include table of this in paper 

################################
##### Assign Taxonomy #######
################################

<<<<<<< HEAD
#assigning ASVs in seqtab.nochim to taxonomy via Silva v132 (most recent version); Silva primary database for 16S 
taxa <- assignTaxonomy(seqtab.nochim,'fastqfiles/silva_nr_v132_train_set.fa.gz',minBoot=50,multithread=TRUE,tryRC=TRUE,outputBootstraps=FALSE)
taxa <- addSpecies(taxa,'fastqfiles/silva_species_assignment_v132.fa.gz') #species level assignments  

#write taxa to csv file 
=======
taxa <- assignTaxonomy(seqtab.nochim, "C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/silva_nr99_v138_train_set.fa.gz", multithread=TRUE)
taxa <- addSpecies(taxa, "C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/silva_species_assignment_v138.fa.gz")
taxa.print <- taxa
rownames(taxa.print) <- NULL
head(taxa.print)
#should we use minboot, multithread, etc in line 153??


#Obtain a csv file for the taxonomy so that it's easier to map the sequences for the heatmap.
>>>>>>> f560edcf6102538485d487966a7bba04e2b1d590
write.csv(taxa, file="taxa.csv",row.name=TRUE,quote=FALSE)
unname(head(taxa, 30))
unname(taxa)

<<<<<<< HEAD
#save reads from sequence table for future analysis 
saveRDS(seqtab.nochim, file="final_seqtab_nochim.rds")
saveRDS(taxa, file="final_taxa_blastCorrected.rds")

#Visualization (phyloseq)

#Construct sample dataframe 

samples.out <- rownames(seqtab.nochim)
subject <- sapply(strsplit(samples.out, "D"), `[`, 1)
gender <- substr(subject,1,1)
subject <- substr(subject,2,999)
day <- as.integer(sapply(strsplit(samples.out, "D"), `[`, 2))
samdf <- data.frame(Subject=subject, Gender=gender, Day=day)
samdf$When <- "Early"
samdf$When[samdf$Day>100] <- "Late"
rownames(samdf) <- samples.out

# Construct phyloseq object (straightforward from dada2 outputs)
ps <- phyloseq(otu_table(seqtab.nochim, taxa_are_rows=FALSE), 
               sample_data(samdf), 
               tax_table(taxa))
ps
=======
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
samdf<-read.csv("variabletable.csv")

head(samdf)
head(seqtab.nochim)
head(taxa)
rownames(samdf) <- samdf$Sample

# Construct phyloseq object (straightforward from dada2 outputs)
otu <- otu_table(seqtab.nochim, taxa_are_rows=FALSE)
head(otu)
samp_data = sample_data(samdf)
head(samp_data)
tax_tab = tax_table(taxa)
head(tax_tab)
ps <- phyloseq(otu, 
               samp_data,
               tax_table(taxa))

head(tax_table(taxa))
View()
ps
sample_names()
>>>>>>> f560edcf6102538485d487966a7bba04e2b1d590

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
<<<<<<< HEAD

#Plot alpha-diversity 
plot_richness(ps, x="Day", measures=c("Shannon", "Simpson"), color="When")

#Bray Curtis 
#Transform to proportions 
ps.prop <- transform_sample_counts(ps, function(otu) otu/sum(otu))
ord.nmds.bray <- ordinate(ps.prop, method="NMDS", distance="bray")
#Plot 
plot_ordination(ps.prop, ord.nmds.bray, color="When", title="Bray NMDS")


                      
=======
>>>>>>> f560edcf6102538485d487966a7bba04e2b1d590
