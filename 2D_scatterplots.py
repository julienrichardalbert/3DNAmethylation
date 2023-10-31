#!/usr/bin/env python
# coding: utf-8

# In[1]:


# I'm gonna make a 2D scatterplot with density contour here.
import sys
import seaborn as sb
import numpy as np
import pandas as pd
import matplotlib
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import matplotlib.font_manager as font_manager
pd.options.mode.chained_assignment = None # stop annoying warning when clipping columns. default='warn'


# In[2]:


#userInputFile = '/Users/jra/Desktop/CTCF_tmp/loops/CTCF_D4_pad1000_FDR0.05_quickAssoc_no11.csv'
#userInputFile = '/Users/jra/Desktop/CTCF_tmp/loops/H3K27ac/ESC_vs_EpiLC_FDR0.05_quickAssoc_no11.csv'
userInputFile = '/Users/jra/Library/CloudStorage/Dropbox/3dname/28_feb_2023/CnR_D0-4_WT-TKO_peaks_RPKM_DNAme_filtered_sort_limma_thresholded.txt'
inputFile = pd.read_csv(userInputFile, sep='\t')
print(inputFile.columns.tolist())


# In[3]:


#x1 = 'CTCF_D4_DnmtWT_rep1.filt.intra'
#x2 = 'CTCF_D4_DnmtWT_rep2.filt.intra'
#y1 = 'CTCF_D4_DnmtTKO_rep1.filt.intra'
#y2 = 'CTCF_D4_DnmtTKO_rep2.filt.intra'
x1 = 'D0_DnmtWT_rep1.filt.intra'
x2 = 'D0_DnmtWT_rep2.filt.intra'
y1 = 'D4_DnmtWT_rep1.filt.intra'
y2 = 'D4_DnmtWT_rep2.filt.intra'
#x1 = 'CTCF_D0_DnmtWT_rep1.filt.intra'
#x2 = 'CTCF_D0_DnmtWT_rep2.filt.intra'
#y1 = 'CTCF_D4_DnmtWT_rep1.filt.intra'
#y2 = 'CTCF_D4_DnmtWT_rep2.filt.intra'

x = 'E14_D4_WT_CTCF_CUTnRUN_AMS042021_rep1-2.bam_RPKM'
y = 'E14_D4_TKO_CTCF_CUTnRUN_AMS042021_rep1-2.bam_RPKM'
#x = 'E14_D0_WT_CTCF_CUTnRUN_AMS042021_rep1-2.bam_RPKM'
#y = 'E14_D0_TKO_CTCF_CUTnRUN_AMS042021_rep1-2.bam_RPKM'

logFC = 'logFC'
#FDR = 'FDR' # HiChIP
FDR = 'adj.P.Val' # CnR
# region = 'region'

#df = inputFile[[x1,x2,y1,y2,logFC,FDR,region]] # HiChIP
df = inputFile[[x,y,logFC,FDR]] # CnR
df


# In[4]:


# sum the replicates for HiChIP data
# df['x'] = df[x1] + df[x2]
# df['y'] = df[y1] + df[y2]

# reverse for non-diffloops comparisons
#df[logFC] = -df[logFC] # diffloops flipped on me
df[logFC] = df[logFC] # diffloops flipped on me

# set the filtering thresholds
FDR_cutoff = 0.05
logFC_cutoff = 1

# set the filtering conditions
conditions = [
    (df[FDR] < FDR_cutoff) & (df[logFC]  > logFC_cutoff),
    (df[FDR] < FDR_cutoff) & (df[logFC] < -logFC_cutoff)
]

# name the conditions
names = ['UP', 'DN']
#colours = {'NN':'0.333, 0.333, 0.333, 0.1', 'UP':'0.333, 0.333, 0.333, 1', 'DN':'0.878, 0.267, 0.278, 1'}
# add a new column with the names that correspond to the conditions met
df['thresh'] = np.select(conditions, names, default='NN')


colors_dict = {'UP':(0.878, 0.267, 0.278, 1.000), 
               'DN':(0.333, 0.333, 0.333, 0.050), 
               'NN':(0.333, 0.333, 0.333, 0.050),}
#colors_dict = {'UP':(0.004, 0.627, 0.451, 1.000), 
#               'DN':(0.62, 0.208, 0.62, 1.000), 
#               'NN':(0.333, 0.333, 0.333, 0.050),}



colors = [mcolors.to_rgba(colors_dict[c]) for c in df['thresh']]
# count the number of instances each condition is met
counts = df['thresh'].value_counts()
print(counts)


# In[7]:


# matplotlib.rcParams['svg.fonttype'] = 'none'  # saves fonts as a text object instead of a vector path
small_size = 12
medium_size = 16
plt.rcParams['font.size'] = medium_size
plt.rcParams["figure.figsize"] = (6,6)
font_path = '/opt/X11/share/system_fonts/HelveticaNeue.ttc'
font_name = 'Helvetica Neue'
prop = font_manager.FontProperties(fname=font_path)

# fig = df.plot(kind="scatter", x = 'x', y = 'y', c = colors, s=12) # HiChIP
fig = df.plot(kind="scatter", x = x, y = y, c = colors, s=12, marker='o',  edgecolors='none')

fig.set_xscale('log')
fig.set_yscale('log')

# this is really stupid, there must be a better way
for idx, count in enumerate(counts):
    fig.text(1.1, 0.95 - idx*0.1, f"{counts.index[idx]}: {count}", ha='left', va='top', transform=fig.transAxes)
 
fig.set_xlabel("ESC contacts", fontname=font_name, fontproperties=prop)
fig.set_ylabel("EpiLC contacts", fontname=font_name, fontproperties=prop)
fig.set_title("HiChIP loops", fontname=font_name, fontproperties=prop)

for tick in fig.get_xticklabels():
    tick.set_fontname(font_name)
    tick.set_fontsize(small_size)
for tick in fig.get_yticklabels():
    tick.set_fontname(font_name)
    tick.set_fontsize(small_size)
    
fig.spines['bottom'].set_linewidth(0.5)
fig.spines['left'].set_linewidth(0.5)
fig.spines['top'].set_linewidth(0)
fig.spines['right'].set_linewidth(0)
fig.tick_params(width=0.5)


# In[8]:


datapoint_count = len(df)
#outFigure = "%s_2D_scatter_%s_datapoints.svg" % (userInputFile, datapoint_count)
#fig.figure.savefig(outFigure)
outFigure = "%s_2D_scatter_%s_datapoints.png" % (userInputFile, datapoint_count)
fig.figure.savefig(outFigure, dpi=1200, bbox_inches='tight', transparent=True)


# In[ ]:




