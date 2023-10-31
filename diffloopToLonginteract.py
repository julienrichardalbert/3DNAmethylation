#!/Users/jra/miniconda3/bin/python3
import sys
import pandas as pd
import numpy as np
import subprocess
pd.set_option('mode.chained_assignment', None)

userInputFile   = sys.argv[1]
inputFile       = pd.read_csv(userInputFile, sep=',')

# print(inputFile.head(1))
read_coverage_per_rep  = inputFile.iloc[:, 6:-9] # get 7th row and all coverage columns (should take all replicates)
# print(read_coverage_per_rep.head(5))
read_coverage          = read_coverage_per_rep.sum(axis=1)

df = inputFile[['chr_1','start_1','end_1','chr_2','start_2','end_2','logFC','FDR']]
df["score"] = read_coverage
#print(df.head(5))

# min FDR in current run is E-7. min is 1.

# BUT I SHOULD NORMALIZE THE RANGE BETWEEN CELL TYPES, SO TAKE THE HIGHEST OF ALL 4 CONDITIONS NOT PAIRWISE!!!
# JUST SET A THRESHOLD ???

df["colour_calc"] = -np.log10(df["FDR"]+0.000000001).round().astype(int)+1 # must add 1 in the case where there are no statistically significant differential loops (as we find between WT and TKO ESCs)
max_value = df["colour_calc"].max()
#print("Maximum -log10 FDR value:", max_value)
df['colour_calc'] = np.where(df['logFC'] > 0, df['colour_calc'], -df['colour_calc'])
df['colour_calc'] = (df['colour_calc']/max_value).round().astype(int) # now the range is -1 to 1


# let's try to make a R-Y-B gradient. np.where(condition, if_true, if_false)
df["colourR"] = np.where( df['logFC'] < 0, 255, 255-255*df['colour_calc'] ) # enriched in TKO. we want it red
df["colourG"] = np.where( df['logFC'] < 0, 255-255*df['colour_calc'].abs(), 255-255*df['colour_calc'].abs() )
df["colourB"] = np.where( df['logFC'] < 0, 0, 255*df['colour_calc'] )

df["colourR"] = df["colourR"].clip(upper=255, lower=0)
df["colourG"] = df["colourG"].clip(upper=255, lower=0)
df["colourB"] = df["colourB"].clip(upper=255, lower=0)

df["color"] = df["colourR"].astype(str) + "," + df["colourG"].astype(str) + "," + df["colourB"].astype(str)
df["NA"] = "."
df["colour_calc"] = df["colour_calc"].clip(upper=1000, lower=0) # has to be in this range to appease UCSC
df["score"] = df["score"].clip(upper=1000, lower=0) # has to be in this range to appease UCSC
df = df.rename(columns={'chr_1': '#chr', 'start_1': 'start', 'end_1': 'end', 'start_1': 'start'})
df_to_print = df[['#chr', 'start', 'end_2', 'NA', 'score', 'colour_calc', 'NA', 'color', '#chr', 'start', 'end', 'NA', 'NA', 'chr_2', 'start_2', 'end_2', 'NA', 'NA']]
# print(df.head(5))
# print(df_to_print.head(2))
outTable = "%s_diffloops_FDR_colors.bed" % (userInputFile)
df_to_print.to_csv(outTable, sep='\t', header=False, index=False)





# range of logFC in current run is -6 to 6
max_value = df["logFC"].abs().max()
#print(max_value)
df['colour_calc'] = (df['logFC']/max_value).round().astype(int) # now the range is -1 to 1


# let's try to make a R-Y-B gradient. np.where(condition, if_true, if_false)
df["colourR"] = np.where( df['logFC'] < 0, 255, 255-255*df['colour_calc'] ) # enriched in TKO. we want it red
df["colourG"] = np.where( df['logFC'] < 0, 255-255*df['colour_calc'].abs(), 255-255*df['colour_calc'].abs() )
df["colourB"] = np.where( df['logFC'] < 0, 0, 255*df['colour_calc'] )

df["colourR"] = df["colourR"].clip(upper=255, lower=0)
df["colourG"] = df["colourG"].clip(upper=255, lower=0)
df["colourB"] = df["colourB"].clip(upper=255, lower=0)

df["color"] = df["colourR"].astype(str) + "," + df["colourG"].astype(str) + "," + df["colourB"].astype(str)
df["NA"] = "."
df["colour_calc"] = df["colour_calc"].clip(upper=1000, lower=0) # has to be in this range to appease UCSC
df["score"] = df["score"].clip(upper=1000, lower=0) # has to be in this range to appease UCSC
df = df.rename(columns={'chr_1': '#chr', 'start_1': 'start', 'end_1': 'end', 'start_1': 'start'})
df_to_print = df[['#chr', 'start', 'end_2', 'NA', 'score', 'colour_calc', 'NA', 'color', '#chr', 'start', 'end', 'NA', 'NA', 'chr_2', 'start_2', 'end_2', 'NA', 'NA']]
#print(df.head(5))
#print(df_to_print.head(2))
outTable = "%s_diffloops_logFC_colors.bed" % (userInputFile)
df_to_print.to_csv(outTable, sep='\t', header=False, index=False)












# set threshold to max adjpval = 0.0001 = 3
max_range = 4 # -log10(0.001) or log2(8)
#remove rows with less than 5 total read alignments
df.drop(df[df['score']<5].index, inplace=True)
# min FDR in current run is E-7. min is 1.
#df["colour_calc"] = -np.log10(df["FDR"]+0.000000001).round().astype(int)+1
df["colour_calc"] = -np.log10(df["FDR"]+0.000000001)
max_value = df["colour_calc"].max()
min_value = df["colour_calc"].min()
print('Initial max FDR value:')
print(max_value)
print('Initial min FDR value:')
print(min_value)
df['colour_calc'] = df['colour_calc'].clip(upper=max_range, lower=0)
max_value = df["colour_calc"].max()
minimum_val = df["colour_calc"].min()
print('Ceiling max FDR value:')
print(max_value)
print('Floor min FDR value:')
print(min_value)
df['colour_calc'] = np.where(df['logFC'] > 0, df['colour_calc'], -df['colour_calc']) # flip the value to reflect UP- or DOWN- in TKO cells
df['colour_calc'] = (df['colour_calc']/max_range) # now the range is -1 to 1
max_value = df["colour_calc"].max()
minimum_val = df["colour_calc"].min()
print('Normalized max FDR value:')
print(max_value)
print('Normalized min FDR value:')
print(min_value)
#print(df.head(5))
df["colourR"] = np.where( df['logFC'] < 0, 255, 255-255*df['colour_calc'] ) # enriched in TKO. we want it red
df["colourG"] = np.where( df['logFC'] < 0, 255-255*df['colour_calc'].abs(), 255-255*df['colour_calc'].abs() )
df["colourB"] = np.where( df['logFC'] < 0, 0, 255*df['colour_calc'] )
df["colourR"] = df["colourR"].clip(upper=255, lower=0).round().astype(int)
df["colourG"] = df["colourG"].clip(upper=255, lower=0).round().astype(int)
df["colourB"] = df["colourB"].clip(upper=255, lower=0).round().astype(int)
df["color"] = df["colourR"].astype(str) + "," + df["colourG"].astype(str) + "," + df["colourB"].astype(str)
df["NA"] = "."
df["colour_calc"] = df["colour_calc"].clip(upper=1000, lower=0) # has to be in this range to appease UCSC
df["score"] = df["score"].clip(upper=1000, lower=0) # has to be in this range to appease UCSC
df = df.rename(columns={'chr_1': '#chr', 'start_1': 'start', 'end_1': 'end', 'start_1': 'start'})
df_to_print = df[['#chr', 'start', 'end_2', 'NA', 'score', 'colour_calc', 'NA', 'color', '#chr', 'start', 'end', 'NA', 'NA', 'chr_2', 'start_2', 'end_2', 'NA', 'NA']]
#print(df.head(5))
#print(df_to_print.head(2))
outTable = "%s_diffloops_x5_FDR_colors.bed" % (userInputFile)
df_to_print.to_csv(outTable, sep='\t', header=False, index=False)



max_value = df["logFC"].abs().max()
print(max_value)
max_range = 3
df['colour_calc'] = df['logFC'].clip(upper=max_range, lower=-max_range)
max_value = df["colour_calc"].abs().max()
print(max_value)
df['colour_calc'] = (df['colour_calc']/max_range) # now the range is -1 to 1
# let's try to make a R-Y-B gradient. np.where(condition, if_true, if_false)
df["colourR"] = np.where( df['logFC'] < 0, 255, 255-255*df['colour_calc'] ) # enriched in TKO. we want it red
df["colourG"] = np.where( df['logFC'] < 0, 255-255*df['colour_calc'].abs(), 255-255*df['colour_calc'].abs() )
df["colourB"] = np.where( df['logFC'] < 0, 0, 255*df['colour_calc'] )
df["colourR"] = df["colourR"].clip(upper=255, lower=0)
df["colourG"] = df["colourG"].clip(upper=255, lower=0)
df["colourB"] = df["colourB"].clip(upper=255, lower=0)
df["color"] = df["colourR"].astype(str) + "," + df["colourG"].astype(str) + "," + df["colourB"].astype(str)
df["NA"] = "."
df["colour_calc"] = df["colour_calc"].clip(upper=1000, lower=0) # has to be in this range to appease UCSC
df["score"] = df["score"].clip(upper=1000, lower=0) # has to be in this range to appease UCSC
df = df.rename(columns={'chr_1': '#chr', 'start_1': 'start', 'end_1': 'end', 'start_1': 'start'})
df_to_print = df[['#chr', 'start', 'end_2', 'NA', 'score', 'colour_calc', 'NA', 'color', '#chr', 'start', 'end', 'NA', 'NA', 'chr_2', 'start_2', 'end_2', 'NA', 'NA']]
outTable = "%s_diffloops_x5_logFC_colors.bed" % (userInputFile)
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
