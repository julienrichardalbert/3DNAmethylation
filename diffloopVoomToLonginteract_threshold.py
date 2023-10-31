#!/Users/jra/miniconda3/bin/python3
import sys
import pandas as pd
import numpy as np
import subprocess
pd.set_option('mode.chained_assignment', None)

userInputFile   = sys.argv[1]
inputFile       = pd.read_csv(userInputFile, sep=',')

# print(inputFile.head(1))
read_coverage_per_rep  = inputFile.iloc[:, 6:-10] # get 7th row and all coverage columns (should take all replicates)
#print(read_coverage_per_rep.head(5))
read_coverage          = read_coverage_per_rep.sum(axis=1)



df = inputFile[['chr_1','start_1','end_1','chr_2','start_2','end_2','logFC','adj.P.Val']]
df["score"] = read_coverage.clip(upper=1000, lower=0) # has to be in this range to appease UCSC
df.drop(df[df['score']<5].index, inplace=True) #remove rows with less than 5 total read alignments

cutoff = 0.05 # min adj.P.Val cutoff
#df["color"] = np.where( df['logFC'] < 0, "246,133,50", "38,103,153" ) # condition, true, false
df["color"] = np.where( df['logFC'] < 0, "1,160,115", "158,53,158" ) # for WT ESC vs EpiLC

df["color"] = np.where( df['adj.P.Val']   < cutoff, df["color"], "222,222,222" ) # set non-sigificant loops to grey

df["NA"] = "." # fill in unnecessary columns
df["colour_calc"] = df["score"] # was previously used for gradient colours

df = df.rename(columns={'chr_1': '#chr', 'start_1': 'start', 'end_1': 'end', 'start_1': 'start'})
df_to_print = df[['#chr', 'start', 'end_2', 'NA', 'score', 'colour_calc', 'NA', 'color', '#chr', 'start', 'end', 'NA', 'NA', 'chr_2', 'start_2', 'end_2', 'NA', 'NA']]
#print(df.head(5))
#print(df_to_print.head(2))
outTable = "%s_diffloops_x5_FDR_colors_thresh.bed" % (userInputFile)
df_to_print.to_csv(outTable, sep='\t', header=False, index=False)







'''
> library(diffloopdata)
> library(ggplot2)
> library(GenomicRanges)
> library(ggrepel)
> library(DESeq2)
> beddir <- "/scratch/hicpro_all_samples_full_depth/test6"
> samples <- c("D4_DnmtWT_rep1.filt.intra","D4_DnmtWT_rep2.filt.intra","D4_DnmtTKO_rep1.filt.intra","D4_DnmtTKO_rep2.filt.intra")
> FullR1 <- loopsMake(beddir,samples=samples)
dim(FullR1)

                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                    anchors interactions samples colData rowData
1   52310       636656       4       2       1
> FullR1_mango <- mangoCorrection(FullR1, FDR = 0.05)
dim(FullR1_mango)







^C
>
>
>
> FullR1_mango <- mangoCorrection(FullR1, FDR = 0.05)

>
> dim(FullR1_mango)
  anchors interactions samples colData rowData
1   50232       340992       4       2       3
> p1 <- loopDistancePlot(FullR1_mango)
> p1
>
>
>
> view(p1)
Error in view(p1) : could not find function "view"
> p1.show()
Error in p1.show() : could not find function "p1.show"
> show(p1)
>
>
> pdf("test_plot.pdf")
> plot(p1)
> dev.off()
pdf
  2
>
>
> loopMetrics(FullR1_mango)
       D4_DnmtWT_rep1.filt.intra D4_DnmtWT_rep2.filt.intra
unique                    268650                    276181
       D4_DnmtTKO_rep1.filt.intra D4_DnmtTKO_rep2.filt.intra
unique                     251928                     257501
> pcp1dat <-  FullR1_mango
> pcp1dat@colData$sizeFactor <- 1
> pcp1 <- pcaPlot(pcp1dat) + geom_text_repel(aes(label=samples)) +
    scale_x_continuous(limits = c(-140, 230)) + ggtitle("PC Plot with no Size Factor Correction") +
    theme(legend.position="none")
>
> pdf("test_plot2.pdf")
> plot(pcp1)
Warning messages:
1: Removed 1 rows containing missing values (geom_point).
2: Removed 1 rows containing missing values (geom_text_repel).
> dev.off()
pdf
  2
> samples
[1] "D4_DnmtWT_rep1.filt.intra"  "D4_DnmtWT_rep2.filt.intra"
[3] "D4_DnmtTKO_rep1.filt.intra" "D4_DnmtTKO_rep2.filt.intra"
> pcp1 <- pcaPlot(pcp1dat) + geom_text_repel(aes(label=samples)) +
    scale_x_continuous(ggtitle("PC Plot with no Size Factor Correction") +
    theme(legend.position="none")
+
+
+
+
> pcp1 <- pcaPlot(pcp1dat) + geom_text_repel(aes(label=samples)) +
     +
    theme(legend.position="none")
Error in `+.gg`:
! Cannot use `+.gg()` with a single argument. Did you accidentally put + on a new line?
Run `rlang::last_error()` to see where the error occurred.
>
> pcp1 <- pcaPlot(pcp1dat) + geom_text_repel(aes(label=samples)) +
    theme(legend.position="none")
> pdf("test_plot2.pdf")
> plot(pcp1)
> dev.off()
pdf
  2
> pcp2 <-pcaPlot(qc_filt) + geom_text_repel(aes(label=samples)) +
    theme(legend.position="none")
Error in h(simpleError(msg, call)) :
  error in evaluating the argument 'dlo' in selecting a method for function 'pcaPlot': object 'qc_filt' not found
> pcp2 <-pcaPlot(qc_filt) + geom_text_repel(aes(label=samples))
Error in h(simpleError(msg, call)) :
  error in evaluating the argument 'dlo' in selecting a method for function 'pcaPlot': object 'qc_filt' not found
> pcp2 <-pcaPlot(qc_filt) + geom_text_repel(aes(label=samples)))
Error: unexpected ')' in "pcp2 <-pcaPlot(qc_filt) + geom_text_repel(aes(label=samples)))"
> pcp2 <-pcaPlot(qc_filt) + geom_text_repel(aes(label=samples))))
Error: unexpected ')' in "pcp2 <-pcaPlot(qc_filt) + geom_text_repel(aes(label=samples)))"
> pcp2 <-pcaPlot(qc_filt) + geom_text_repel(aes(label=samples))
Error in h(simpleError(msg, call)) :
  error in evaluating the argument 'dlo' in selecting a method for function 'pcaPlot': object 'qc_filt' not found
> loopMetrics(FullR1_mango)
       D4_DnmtWT_rep1.filt.intra D4_DnmtWT_rep2.filt.intra
unique                    268650                    276181
       D4_DnmtTKO_rep1.filt.intra D4_DnmtTKO_rep2.filt.intra
unique                     251928                     257501
>
> client_loop: send disconnect: Broken pipe
(base) jra@greenberg14:~/Desktop$ ijm


ssh: connect to host epipax.ijm.univ-paris-diderot.priv port 22: Operation timed out
'''




#Information on file formats:
#https://genome.ucsc.edu/goldenPath/help/interact.html

#"interaction between two regions"
#  string chrom;        "Chromosome (or contig, scaffold, etc.). For interchromosomal, use 2 records"
#  uint chromStart;     "Start position of lower region. For interchromosomal, set to chromStart of this region"
#  uint chromEnd;       "End position of upper region. For interchromosomal, set to chromEnd of this region"
#  string name;         "Name of item, for display.  Usually 'sourceName/targetName/exp' or empty"
#  uint score;          "Score (0-1000)"
#  double value;        "Strength of interaction or other data value. Typically basis for score"
#  string exp;          "Experiment name (metadata for filtering). Use . if not applicable"
#  string color;        "Item color.  Specified as r,g,b or hexadecimal #RRGGBB or html color name, as in //www.w3.org/TR/css3-color/#html4. Use 0 and spectrum setting to shade by score"
#  string sourceChrom;  "Chromosome of source region (directional) or lower region. For non-directional interchromosomal, chrom of this region."
#  uint sourceStart;    "Start position in chromosome of source/lower/this region"
#  uint sourceEnd;      "End position in chromosome of source/lower/this region"
#  string sourceName;   "Identifier of source/lower/this region"
#  string sourceStrand; "Orientation of source/lower/this region: + or -.  Use . if not applicable"
#  string targetChrom;  "Chromosome of target region (directional) or upper region. For non-directional interchromosomal, chrom of other region"
#  uint targetStart;    "Start position in chromosome of target/upper/this region"
#  uint targetEnd;      "End position in chromosome of target/upper/this region"
#  string targetName;   "Identifier of target/upper/this region"
#  string targetStrand; "Orientation of target/upper/this region: + or -.  Use . if not applicable"
