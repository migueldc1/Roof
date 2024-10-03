# Author: M, de Celis Rodriguez
# Date: 19/09/2023
# Project: Marco - ASV sequence analysis (ITS2)

library(dada2)
library(ggplot2)
library(ShortRead)
library(Biostrings)


rm(list=ls()) #Clear R environment

# Set the project location as working directory
setwd("~/../OneDrive/Proyecto - MadreenRoof/")

#
#### CUSTOM FUNCTIONS ####
getN <- function(x) sum(getUniques(x))

allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

primerHits <- function(primer, fn) {
  # Counts number of reads in which the primer is found
  nhits <- vcountPattern(primer, sread(readFastq(fn)), fixed = FALSE)
  return(sum(nhits > 0))
}

#
#### GETTING READY ####
path <- "Reads/ITS" # CHANGE ME to the directory containing the fastq files after unzipping.

fnFs <- sort(list.files(path, pattern = ".fastq", full.names = TRUE))
fnFs <- fnFs[grepl("R1", fnFs)]
fnRs <- sort(list.files(path, pattern = ".fastq", full.names = TRUE))
fnRs <- fnRs[grepl("R2", fnRs)]

sample.names <- sapply(strsplit(basename(fnFs), "_"), `[`, 1)
sample.names <- gsub("NGS019-23-ITS2-", "", sample.names)

#Identify Primers
FWD <- "GTGARTCATCGAATCTTTG"
REV <- "TCCTCCGCTTATTGATATGC"

FWD.orients <- allOrients(FWD)
REV.orients <- allOrients(REV)

#
#### REMOVE Ns ####
fnFs.filtN <- file.path(path, "Remove.primers/filtN", basename(fnFs))
fnRs.filtN <- file.path(path, "Remove.primers/filtN", basename(fnRs))

filtN <- filterAndTrim(fnFs, fnFs.filtN, fnRs, fnRs.filtN, maxN = 0, verbose = TRUE)

save.image("Outputs/DADA2_fun.RData")

#
#### REMOVE PRIMERS ####
rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.filtN[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.filtN[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.filtN[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.filtN[[2]]))

#Set cutadapt path
cutadapt <- "D://Software/cutadapt.exe"
system2(cutadapt, args = "--version")

fnFs.cut <- file.path(path, "Remove.primers/cutadapt", basename(fnFs))
fnRs.cut <- file.path(path, "Remove.primers/cutadapt", basename(fnRs))

# Run Cutadapt
for(i in seq_along(fnFs)) {
  system2(cutadapt, args = c("-g", FWD, "-G", REV,
                             "-o", fnFs.cut[i], "-p", fnRs.cut[i],
                             fnFs.filtN[i], fnRs.filtN[i]))
}

rbind(FWD.ForwardReads = sapply(FWD.orients, primerHits, fn = fnFs.cut[[2]]), 
      FWD.ReverseReads = sapply(FWD.orients, primerHits, fn = fnRs.cut[[2]]), 
      REV.ForwardReads = sapply(REV.orients, primerHits, fn = fnFs.cut[[2]]), 
      REV.ReverseReads = sapply(REV.orients, primerHits, fn = fnRs.cut[[2]]))


#
#### INSPECT QUALITY PROFILES ####
qprFs <- plotQualityProfile(fnFs.cut[sample(1:length(sample.names), 3)])
qprRs <- plotQualityProfile(fnRs.cut[sample(1:length(sample.names), 3)])

dataF <- NULL
dataR <- NULL

for (cyc in 1:max(qprFs$data$Cycle)) {
  
  badF <- sum(qprFs$data[qprFs$data$Cycle == cyc,]$Count[qprFs$data[qprFs$data$Cycle == cyc,]$Score < 30])
  totF <- sum(qprFs$data[qprFs$data$Cycle == cyc,]$Count)
  
  dataF <- rbind(dataF, cbind.data.frame(Cycle = cyc, Count = badF, Proportion = badF*100/totF))
  
  badR <- sum(qprRs$data[qprRs$data$Cycle == cyc,]$Count[qprRs$data[qprRs$data$Cycle == cyc,]$Score < 30])
  totR <- sum(qprRs$data[qprRs$data$Cycle == cyc,]$Count)
  
  dataR <- rbind(dataR, cbind.data.frame(Cycle = cyc, Count = badR, Proportion = badR*100/totR))
  
}

ggplot() +
  geom_point(data = dataF, aes(x = Cycle, y = 100-Proportion), color = "blue") +
  geom_point(data = dataR, aes(x = Cycle, y = 100-Proportion), color = "red") +
  scale_x_continuous(breaks = seq(0, 300, 20))

#
#### FILTER AND TRIM ####
filtFs <- file.path(path, "Remove.primers/filtered2", paste0(sample.names, "_F_filt.fastq.gz"))
filtRs <- file.path(path, "Remove.primers/filtered2", paste0(sample.names, "_R_filt.fastq.gz"))

out <- filterAndTrim(fnFs.cut, filtFs, fnRs.cut, filtRs, truncLen = c(250, 185), maxEE = 6, truncQ = 2,
                     rm.phix = TRUE, compress = TRUE, verbose = TRUE, multithread = FALSE) # On Windows set multithread = FALSE

row.names(out) <- sapply(strsplit(row.names(out), "_"), `[`, 1)
row.names(out) <- gsub("NGS019-23-ITS2-", "", row.names(out))
head(out)


#
#### LEARN ERROR RATES ####
errF <- learnErrors(filtFs, multithread = FALSE)
errR <- learnErrors(filtRs, multithread = FALSE)
plotErrors(errF, nominalQ = TRUE)


#
#### SAMPLE INFERENCE ####
dadaFs <- dada(filtFs, err = errF, multithread = FALSE)
dadaRs <- dada(filtRs, err = errR, multithread = FALSE)


#
#### MERGE PAIRED READS ####
mergers <- mergePairs(dadaFs, filtFs, dadaRs, filtRs, verbose = TRUE, minOverlap = 8, maxMismatch = 2)
head(mergers[[1]])


#
#### CONSTRUCT SEQUENCE TABLE ####
seqtab <- makeSequenceTable(mergers)
dim(seqtab)
hist(nchar(getSequences(seqtab)))


#
## REMOVE CHIMERAS
seqtab.nochim <- removeBimeraDenovo(seqtab, method = "consensus", multithread = FALSE, verbose = TRUE)
dim(seqtab.nochim)
sum(seqtab.nochim)/sum(seqtab)


#
#### TRACK READS ####
track <- cbind.data.frame(out, sapply(dadaFs, getN), sapply(dadaRs, getN), sapply(mergers, getN), rowSums(seqtab.nochim))
colnames(track) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim")
rownames(track) <- sample.names
head(track)

track$Proportion <- track$nonchim*100/track$input


#
#### ASSIGN TAXONOMY ####
tax_fun <- assignTaxonomy(seqtab.nochim, "Database/sh_general_release_dynamic_25.07.2023_dev.fasta", outputBootstraps = FALSE, minBoot = 80)
tax_fun <- as.data.frame(tax_fun)

#
#### SAVE DATA ####
saveRDS(seqtab.nochim, "Data/Sequencing/Outputs/ASV.raw_fun.rds")
saveRDS(tax_fun, "Data/Sequencing/Outputs/tax.raw_fun.rds")
uniquesToFasta(seqtab.nochim, fout = "Data/Sequencing/Outputs/uniques_fun.fna", ids = colnames(seqtab.nochim))
write.table(track, "Data/Sequencing/Outputs/track_fun.txt", row.names = TRUE, sep = "\t")

save.image("Data/Sequencing/Outputs/DADA2_fun.RData")

#