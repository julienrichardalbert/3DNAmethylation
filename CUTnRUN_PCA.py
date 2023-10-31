#!/usr/bin/env python
# coding: utf-8

# In[101]:


import sys
import numpy as np
import pandas as pd
import seaborn as sb
import matplotlib
import matplotlib.pyplot as plt
pd.options.mode.chained_assignment = None # stop annoying warning when clipping columns. default='warn'
from combat.pycombat import pycombat


# In[97]:


# prepare data
# the indexes correspond to the gene names
# the column names correspond to the sample names

# CUT n RUN
input_df = pd.read_csv("/Users/jra/Desktop/CTCF_tmp/data_filtered_sort_limma.txt",
                            sep='\t', index_col=0)


# HiChIP 1D
#input_df = pd.read_csv("/Users/jra/Dropbox/3dname/28_feb_2023/E14_D0-4_WT-TKO_CTCF_HiChIP_rep1-2_R1-2_trimV1_mm10_macs2_broad_peaks_DNAme_filtered_1CpG_data_named.txt",
#                            sep='\t', index_col=0)
    




for column_name in input_df.columns:
    print(column_name)
print(input_df.shape)
#print(df_expression.columns)
#print(df_expression.head(2))


# In[98]:


# CUT n RUN
df = input_df[["E14_D0_WT_CTCF_CUTnRUN_AMS042021_rep1_crop36_trimV5_mm10.bam_RPKM",
               "E14_D0_WT_CTCF_CUTnRUN_AMS042021_rep2_crop36_trimV5_mm10.bam_RPKM",
               "E14_D0_TKO_CTCF_CUTnRUN_AMS042021_rep1_crop36_trimV5_mm10.bam_RPKM",
               "E14_D0_TKO_CTCF_CUTnRUN_AMS042021_rep2_crop36_trimV5_mm10.bam_RPKM",
               "E14_D4_WT_CTCF_CUTnRUN_AMS042021_rep1_crop36_trimV5_mm10.bam_RPKM",
               "E14_D4_WT_CTCF_CUTnRUN_AMS042021_rep2_crop36_trimV5_mm10.bam_RPKM",
               "E14_D4_TKO_CTCF_CUTnRUN_AMS042021_rep1_crop36_trimV5_mm10.bam_RPKM",
               "E14_D4_TKO_CTCF_CUTnRUN_AMS042021_rep2_crop36_trimV5_mm10.bam_RPKM"]]
 

    
# HiChIP 1D
#df = input_df[["E14_D4_WT_CTCF_HiChIP_deep_shallow_AMS112022_rep1_R1-2_trimV1_mm10.bam_RPKM",
#               "E14_D0_WT_CTCF_HiChIP_deep_shallow_AMS112022_rep1_R1-2_trimV1_mm10.bam_RPKM",
#               "E14_D4_TKO_CTCF_HiChIP_deep_shallow_AMS112022_rep1_R1-2_trimV1_mm10.bam_RPKM",
#               "E14_D0_TKO_CTCF_HiChIP_deep_shallow_AMS112022_rep1_R1-2_trimV1_mm10.bam_RPKM",
#               "E14_D4_TKO_CTCF_HiChIP_AMS012023_rep2_R1-2_trimV1_mm10.bam_RPKM",
#               "E14_D0_WT_CTCF_HiChIP_AMS012023_rep2_R1-2_trimV1_mm10.bam_RPKM",
#               "E14_D4_WT_CTCF_HiChIP_AMS012023_rep2_R1-2_trimV1_mm10.bam_RPKM",
#               "E14_D0_TKO_CTCF_HiChIP_AMS012023_rep2_R1-2_trimV1_mm10.bam_RPKM"]]

df
#df_expression
# SOMETHING IS WRONG WITH THE TABLE LETS FIND OUT WHAT
# df_expression.isnull().values.any()
# df_expression.isnull().sum()
# null_data = df_expression[df_expression.isnull().any(axis=1)]
# null_data
# should return 0 rows if the table is nice

df = df.loc[(df!=0).any(1)] # drop rows that have all 0 values
print("Dropping rows with only 0 values")
print(df.shape)
df


# In[99]:


matplotlib.rcParams['svg.fonttype'] = 'none'  # saves fonts as a text object instead of a vector path
plt.rcParams['font.size'] = 8
plt.rcParams["figure.figsize"] = (20,4)

# plot the distribution of log2(RPKM+1) gene expression values (uncorrected for batch effects)
boxPlot = plt.boxplot(np.log2(df.clip(upper=None, lower=0)+1), \
                    labels=df.columns, \
                    patch_artist=True, notch=True, \
                    medianprops={'color':'Black'}, \
                    flierprops={'markersize':2, 'marker':'o', 'markerfacecolor':'white', 'markeredgecolor':'black', 'markeredgewidth':0.5})
plt.xticks(rotation=90)
plt.ylim(-2.1, None)
plt.show() # showing the figure "moves" the figure to Shell, leaving nothing behind. Show after saving file.


# In[100]:


from sklearn.decomposition import PCA

samples = df.transpose()
#samples = df
pca = PCA(2) 
projected = pca.fit_transform(samples)
print(df.shape)
print(samples.shape)
print(projected.shape)

pca1variation = (pca.explained_variance_ratio_ * 100)[0].round(2)
pca2variation = (pca.explained_variance_ratio_ * 100)[1].round(2)

# use the explained variance to size the output graph
plt.rcParams["figure.figsize"] = ((pca1variation/4).astype(int), (pca2variation/4).astype(int))
plt.scatter(projected[:, 0], projected[:, 1],
            edgecolor='none', alpha=1, s=200)

x_title = 'PC1: ' + pca1variation.astype(str) + '% of variation explained'
y_title = 'PC2: ' + pca2variation.astype(str) + '% of variation explained'
plt.xlabel(x_title)
plt.ylabel(y_title)
plt.title('made by Julien using python in a jupyter notebook')

for i, label in enumerate(df.columns):
    plt.annotate(label, (projected[:,0][i], projected[:,1][i]))
outFigure = "/Users/jra/Desktop/PCA_plot_PCscaled.svg"
plt.savefig(outFigure)
plt.close()

plt.rcParams["figure.figsize"] = (4,4)
plt.scatter(projected[:, 0], projected[:, 1],
            edgecolor='none', alpha=1, s=200)

x_title = 'PC1: ' + pca1variation.astype(str) + '% of variation explained'
y_title = 'PC2: ' + pca2variation.astype(str) + '% of variation explained'
plt.xlabel(x_title)
plt.ylabel(y_title)
plt.title('made by Julien using python in a jupyter notebook')

for i, label in enumerate(df.columns):
    plt.annotate(label, (projected[:,0][i], projected[:,1][i]))
outFigure = "/Users/jra/Desktop/PCA_plot.svg"
plt.savefig(outFigure)
#plt.close()

