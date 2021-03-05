install.packages("cutadapt")

library(dada2); #packageVersion("dada2"); citation("dada2")
library(ShortRead); #packageVersion("ShortRead")
library(ggplot2); #packageVersion("ggplot2")
library(phyloseq); #packageVersion("phyloseq")

getwd()
#setwd("C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/fastqfiles")

path <- "/Users/gracebeery/Desktop/BI586/fishmicrobes_2021_assign1" 
#path <- "C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1/fastqfiles"
fns <- list.files(path)
#Let's make sure that all of our files are there
fns

fastqs <- fns[grepl(".fastq$", fns)]
fastqs <- sort(fastqs) # Sort ensures reads are in same order

sample.names <- sapply(strsplit(fastqs, ".fastq"), `[`, 1) #the last number will select the field for renaming
sample.names
# Specify the full path to the fnFs
fnFs <- file.path(path, fastqs)
fnFs
plotQualityProfile(fnFs[c(1,2,3,4,5,6,7,8,9)])
plotQualityProfile(fnFs[c(10,11,12,13,14,15,16,17,18)])
plotQualityProfile(fnFs[c(19,20,21,22,23,24,25,26,27)])
plotQualityProfile(fnFs[c(28,29,30,31,32,33,34,35,36)])

filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))

#edit these numbers
#change trim to match primer for our study
out <- filterAndTrim(fnFs, filtFs, truncLen= 300, #change this, should be less than 300
                     maxN=0, #DADA does not allow Ns
                     maxEE=1, #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     trimLeft=20, #N nucleotides to remove from the start of each read: ITS2 primer = F 20bp
                     rm.phix=TRUE, #remove reads matching phiX genome
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

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
setDadaOpt(BAND_SIZE=32)
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



