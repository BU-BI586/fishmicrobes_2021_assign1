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
    
if (!requireNamespace("BiocManager", quietly = TRUE))
    install.packages("BiocManager")
    BiocManager::install("cutadapt") 
  
library(dada2); #packageVersion("dada2"); citation("dada2")
library(ShortRead); #packageVersion("ShortRead")
library(ggplot2); #packageVersion("ggplot2")
library(phyloseq); #packageVersion("phyloseq")

getwd()
#setwd("C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/fastqfiles")

#path <- "/Users/gracebeery/Desktop/BI586/fishmicrobes_2021_assign1" 
#path <- "C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/fastqfiles"
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

#edit these numbers
#change trim to match primer for our study
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

plotErrors(errF, nominalQ=TRUE) #plotting error rates, want to see mostly linear relationships

#dereplicating reads - collapsing all reads into unique sequences
derepFs <- derepFastq(filtFs, verbose=TRUE)
# Name the derep-class objects by the sample names
names(derepFs) <- sample.names

#infer sequence variants 
#Must change some of the DADA options b/c original program optomized for ribosomal data, not ITS - from github, "We currently recommend BAND_SIZE=32 for ITS data." leave as default for 16S/18S
setDadaOpt(BAND_SIZE=16) #Default is 16 for Illumina (low rates of indels) according to documentation 
dadaFs <- dada(derepFs, err=errF, multithread=TRUE)

#now, look at the dada class objects by sample
#will tell how many 'real' variants in unique input seqs
#By default, the dada function processes each sample independently, but pooled processing is available with pool=TRUE and that may give better results for low sampling depths at the cost of increased computation time. See our discussion about pooling samples for sample inference. 
dadaFs[[1]]
dadaFs[[25]]

#construct sequence table
seqtab <- makeSequenceTable(dadaFs)
head(seqtab)

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
rownames(track) <- sample.names
head(track)
tail(track)

write.csv(track,file="ReadFilterStats_AllData_final.csv",row.names=TRUE,quote=FALSE)
#this shows raw reads, how many filtered, denoised, etc.
#could include table of this in paper 

################################
##### Assign Taxonomy #######
################################



