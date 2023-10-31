#!/Users/jra/miniconda3/bin/python3
import sys
import pandas as pd
import numpy as np
import subprocess
pd.set_option('mode.chained_assignment', None)

userInputFile   = sys.argv[1]
inputFile       = pd.read_csv(userInputFile, sep=',')

#print(inputFile.head(1))
read_coverage_per_rep = inputFile.iloc[:, 6:-4] # get 7th row and all coverage columns (should take all replicates)
read_coverage = read_coverage_per_rep.sum(axis=1)

df = inputFile.iloc[:, [0,1,2,3,4,5,-3]]



df = df.rename(columns={'chr_1': '#chr', 'start_1': 'start', 'end_1': 'end', 'start_1': 'start'})
df_to_print = df[['#chr', 'start', 'end', 'chr_2', 'start_2', 'end_2', 'mango.FDR']]
#print(df.head(5))
#print(df_to_print.head(2))

outTable = "%s_longRangeInteract.bedGraph" % (userInputFile)
df_to_print.to_csv(outTable, sep='\t', header=False, index=False)




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
