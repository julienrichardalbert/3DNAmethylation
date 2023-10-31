library(diffloop)
library(diffloopdata)
library(ggplot2)
library(GenomicRanges)
library(ggrepel)
library(DESeq2)

# CODE KINDLY PROVIDED BY Ana Paula Azambuja AND ADAPTED BY JRA
#load bedpe directory and make loops object
beddir <- "/scratch/hicpro_all_samples_full_depth/D4_hichipper_prepeaks"
samplesWT <- c("D4_DnmtWT_rep1.filt.intra","D4_DnmtWT_rep2.filt.intra")
samplesTKO <- c("D4_DnmtTKO_rep1.filt.intra","D4_DnmtTKO_rep2.filt.intra")
WT <- loopsMake(beddir,samples=samplesWT)
dim(WT)
TKO <- loopsMake(beddir,samples=samplesTKO)
dim(TKO)

#run mango correction (FDR>0.05)
WT_mango <- mangoCorrection(WT, FDR = 0.05)
dim(WT_mango)
TKO_mango <- mangoCorrection(TKO, FDR = 0.05)
dim(TKO_mango)


#convert to and write table
WT_table <- summary(WT_mango)
TKO_table <- summary(TKO_mango)
setwd(dir = beddir)
write.table(WT_table, file = "D4_WT_FDR0.05_mango.csv", append = FALSE, quote = FALSE, row.names = FALSE, sep = ",")
write.table(TKO_table, file = "D4_TKO_FDR0.05_mango.csv", append = FALSE, quote = FALSE, row.names = FALSE, sep = ",")
