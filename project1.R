install.packages("cutadapt")

library(dada2); #packageVersion("dada2"); citation("dada2")
library(ShortRead); #packageVersion("ShortRead")
library(ggplot2); #packageVersion("ggplot2")
library(phyloseq); #packageVersion("phyloseq")

getwd()
path <- "/Users/gracebeery/Desktop/BI586/fishmicrobes_2021_assign1" 
path <- "C:/Users/Maddy/Documents/BI586/fishmicrobes_2021_assign1"
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
plotQualityProfile(fnFs[c(1,2,3,4,5)])
plotQualityProfile(fnFs[c(6,7,8,9,10)])
plotQualityProfile(fnFs[c(11,12,13,14,15)])
plotQualityProfile(fnFs[c(16,17,18,19,20)])

filt_path <- file.path(path, "trimmed")
if(!file_test("-d", filt_path)) dir.create(filt_path)
filtFs <- file.path(filt_path, paste0(sample.names, "_F_filt.fastq.gz"))

#edit these numbers
out <- filterAndTrim(fnFs, filtFs, truncLen= 300, #end of single end reads = approx. 300 bp
                     maxN=0, #DADA does not allow Ns
                     maxEE=1, #allow 1 expected errors, where EE = sum(10^(-Q/10)); more conservative, model converges
                     truncQ=2, 
                     trimLeft=20, #N nucleotides to remove from the start of each read: ITS2 primer = F 20bp
                     rm.phix=TRUE, #remove reads matching phiX genome
                     compress=TRUE, multithread=FALSE) # On Windows set multithread=FALSE

head(out)
tail(out)

setDadaOpt(MAX_CONSIST=30) #increase number of cycles to allow convergence
errF <- learnErrors(filtFs, multithread=TRUE)

plotErrors(errF, nominalQ=TRUE) 


