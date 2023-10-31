library(diffloop)
library(diffloopdata)
library(ggplot2)
library(GenomicRanges)
library(ggrepel)
library(DESeq2)

# CODE KINDLY PROVIDED BY Ana Paula Azambuja AND ADAPTED BY JRA
#load bedpe directory and make loops object
beddir <- "/scratch/hicpro_all_samples_full_depth/D0_hichipper_prepeaks"
samples <- c("D0_DnmtWT_rep1.filt.intra","D0_DnmtWT_rep2.filt.intra","D0_DnmtTKO_rep1.filt.intra","D0_DnmtTKO_rep2.filt.intra")
FullR1 <- loopsMake(beddir,samples=samples)
dim(FullR1)

#run mango correction (FDR>0.05)
FullR1_mango <- mangoCorrection(FullR1, FDR = 0.05)
dim(FullR1_mango)

#determine groups
groups <- c("WT","WT","TKO","TKO")
FullR1_mango_groups <- updateLDGroups(FullR1_mango, groups)
dim(FullR1_mango_groups)
res <- quickAssoc(FullR1_mango_groups)
#convert to and write table
res_table <- summary(res)
setwd(dir = beddir)
write.table(res_table, file = "D0_FDR0.05_quickAssoc.csv", append = FALSE, quote = FALSE, row.names = FALSE, sep = ",")
