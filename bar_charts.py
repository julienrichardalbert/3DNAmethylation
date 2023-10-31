#!/usr/bin/env python
# coding: utf-8

# In[ ]:


# 23 Jan 2023
# JRA wrote this with the help of ChatGPT
# The prompt and example code are included at the end of this notebook. Thanks, ChatGPT bot!


# In[5]:


# Load libraries
import matplotlib.pyplot
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import ttest_ind

matplotlib.rcParams['svg.fonttype'] = 'none'  # saves fonts as a text object instead of a vector path
plt.rcParams['font.size'] = 8


# In[6]:


# Load data
#userInputFile = '/Users/jra/Dropbox/greenberg05_julien_share/epigenetic_switching/tables/rna/mm10.refGene.denovo.data.Visr0.15_join.txt'
#userInputFile = '/Users/jra/Dropbox/greenberg05_julien_share/epigenetic_switching/tables/rna/mm10.refGene.denovo.data.Visr0.15_join.txt'
userInputFile = '/Users/jra/Library/CloudStorage/Dropbox/greenberg05_julien_share/epigenetic_switching/tables/rna/mm10_refseq_genes_visr_RNA_data.txt'


inputFile = pd.read_csv(userInputFile, sep='\t')
inputFile

for column_name in inputFile.columns:
    print(column_name)


# In[ ]:


# List relevant groups
'''
E14_WT_D0_RNA_D129T01_rep1_trimV5_mm10.bam_RPKM
E14_WT_D0_RNA_D129T09_rep2_trimV5_mm10.bam_RPKM

E14_TKO_D0_RNA_D129T05_rep1_trimV5_mm10.bam_RPKM
E14_TKO_D0_RNA_D129T13_rep2_trimV5_mm10.bam_RPKM

E14_WT_D2_RNA_D129T02_rep1_trimV5_mm10.bam_RPKM
E14_WT_D2_RNA_D129T10_rep2_trimV5_mm10.bam_RPKM

E14_TKO_D2_RNA_D129T06_rep1_trimV5_mm10.bam_RPKM
E14_TKO_D2_RNA_D129T14_rep2_trimV5_mm10.bam_RPKM

E14_WT_D4_RNA_D129T03_rep1_trimV5_mm10.bam_RPKM
E14_WT_D4_RNA_D129T11_rep2_trimV5_mm10.bam_RPKM

E14_TKO_D4_RNA_D129T07_rep1_trimV5_mm10.bam_RPKM
E14_TKO_D4_RNA_D129T15_rep2_trimV5_mm10.bam_RPKM

E14_WT_D7_RNA_D129T04_rep1_trimV5_mm10.bam_RPKM
E14_WT_D7_RNA_D129T12_rep2_trimV5_mm10.bam_RPKM

E14_TKO_D7_RNA_D129T08_rep1_trimV5_mm10.bam_RPKM
E14_TKO_D7_RNA_D129T16_rep2_trimV5_mm10.bam_RPKM

E14_D7_DnmtTKO_EzhipWT_RNA_ams78_trimV5_mm10.bam_RPKM
E14_D7_DnmtTKO_EzhipWT_RNA_ams77_trimV5_mm10.bam_RPKM

E14_D7_DnmtTKO_EzhipKO_RNA_ams80_trimV5_mm10.bam_RPKM
E14_D7_DnmtTKO_EzhipKO_RNA_ams79_trimV5_mm10.bam_RPKM
'''


# In[7]:


# Assign relevant data to groups
gene_names = ['name']
D0_WT = ['E14_WT_D0_RNA_D129T01_rep1_trimV5_mm10.bam_RPKM', 'E14_WT_D0_RNA_D129T09_rep2_trimV5_mm10.bam_RPKM']
D0_TKO = ['E14_TKO_D0_RNA_D129T05_rep1_trimV5_mm10.bam_RPKM', 'E14_TKO_D0_RNA_D129T13_rep2_trimV5_mm10.bam_RPKM']
D2_WT = ['E14_WT_D2_RNA_D129T02_rep1_trimV5_mm10.bam_RPKM', 'E14_WT_D2_RNA_D129T10_rep2_trimV5_mm10.bam_RPKM']
D2_TKO = ['E14_TKO_D2_RNA_D129T06_rep1_trimV5_mm10.bam_RPKM', 'E14_TKO_D2_RNA_D129T14_rep2_trimV5_mm10.bam_RPKM']
D4_WT = ['E14_WT_D4_RNA_D129T03_rep1_trimV5_mm10.bam_RPKM', 'E14_WT_D4_RNA_D129T11_rep2_trimV5_mm10.bam_RPKM']
D4_TKO = ['E14_TKO_D4_RNA_D129T07_rep1_trimV5_mm10.bam_RPKM', 'E14_TKO_D4_RNA_D129T15_rep2_trimV5_mm10.bam_RPKM']
D7_WT = ['E14_WT_D7_RNA_D129T04_rep1_trimV5_mm10.bam_RPKM', 'E14_WT_D7_RNA_D129T12_rep2_trimV5_mm10.bam_RPKM']
D7_TKO = ['E14_TKO_D7_RNA_D129T08_rep1_trimV5_mm10.bam_RPKM', 'E14_TKO_D7_RNA_D129T16_rep2_trimV5_mm10.bam_RPKM']
D7_TKO_EzhipWT = ['E14_D7_DnmtTKO_EzhipWT_RNA_ams78_trimV5_mm10.bam_RPKM', 'E14_D7_DnmtTKO_EzhipWT_RNA_ams77_trimV5_mm10.bam_RPKM']
D7_TKO_EzhipKO = ['E14_D7_DnmtTKO_EzhipKO_RNA_ams80_trimV5_mm10.bam_RPKM', 'E14_D7_DnmtTKO_EzhipKO_RNA_ams79_trimV5_mm10.bam_RPKM']


# In[8]:


# Create dataframe with relevant groups
df = inputFile[gene_names + \
               D0_WT + D0_TKO + \
               D2_WT + D2_TKO + \
               D4_WT + D4_TKO + \
               D7_WT + D7_TKO + \
               D7_TKO_EzhipWT + D7_TKO_EzhipKO]
df



# In[21]:


# Filter the dataframe to only include the row corresponding to gene X
# gene_id = 'NM_022435'
gene_name = 'Mob3b'
gene = df.loc[df['name'] == gene_name ]
#gene

# Get the actual values, not some pandas series thing I don't understand
D0_WT_vals=gene[D0_WT].iloc[0].tolist()
D0_TKO_vals=gene[D0_TKO].iloc[0].tolist()
D2_WT_vals=gene[D2_WT].iloc[0].tolist()
D2_TKO_vals=gene[D2_TKO].iloc[0].tolist()
D4_WT_vals=gene[D4_WT].iloc[0].tolist()
D4_TKO_vals=gene[D4_TKO].iloc[0].tolist()
D7_WT_vals=gene[D7_WT].iloc[0].tolist()
D7_TKO_vals=gene[D7_TKO].iloc[0].tolist()
D7_TKO_EzhipWT_vals=gene[D7_TKO_EzhipWT].iloc[0].tolist()
D7_TKO_EzhipKO_vals=gene[D7_TKO_EzhipKO].iloc[0].tolist()

# Calculate the mean of each group
D0_WT_mean=np.mean(D0_WT_vals)
D0_TKO_mean=np.mean(D0_TKO_vals)
D2_WT_mean=np.mean(D2_WT_vals)
D2_TKO_mean=np.mean(D2_TKO_vals)
D4_WT_mean=np.mean(D4_WT_vals)
D4_TKO_mean=np.mean(D4_TKO_vals)
D7_WT_mean=np.mean(D7_WT_vals)
D7_TKO_mean=np.mean(D7_TKO_vals)
D7_TKO_EzhipWT_mean=np.mean(D7_TKO_EzhipWT_vals)
D7_TKO_EzhipKO_mean=np.mean(D7_TKO_EzhipKO_vals)
# print(D0_WT_mean, D0_TKO_mean, D2_WT_mean, D2_TKO_mean, D4_WT_mean, D4_TKO_mean, D7_WT_mean, D7_TKO_mean, D7_TKO_EzhipWT_mean, D7_TKO_EzhipTKO_mean)

# Create the base bar plot
x = [1,2,4,5,7,8,10,11]
plt.bar(x, [D0_WT_mean, D0_TKO_mean, D2_WT_mean, D2_TKO_mean, D4_WT_mean, D4_TKO_mean, D7_WT_mean, D7_TKO_mean], \
        color=['#3B90AF', '#F9C052', '#3584A8', '#F8A743', '#2D75A1', '#F7963C', '#266799', '#F68532', '#F68532', '#F5584C'],  edgecolor="black", linewidth=0.5, width=0.9)

# Remove the top and right 'axes'
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)

# Add jitter to the scatter plot of individual data points
jitter = np.random.uniform(-0.1,0.1, size=2)

plt.scatter([x[0]+jitter[0]], gene[D0_WT[0]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[0]+jitter[1]], gene[D0_WT[1]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[1]+jitter[0]], gene[D0_TKO[0]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[1]+jitter[1]], gene[D0_TKO[1]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)

plt.scatter([x[2]+jitter[0]], gene[D2_WT[0]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[2]+jitter[1]], gene[D2_WT[1]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[3]+jitter[0]], gene[D2_TKO[0]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[3]+jitter[1]], gene[D2_TKO[1]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)

plt.scatter([x[4]+jitter[0]], gene[D4_WT[0]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[4]+jitter[1]], gene[D4_WT[1]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[5]+jitter[0]], gene[D4_TKO[0]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[5]+jitter[1]], gene[D4_TKO[1]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)

plt.scatter([x[6]+jitter[0]], gene[D7_WT[0]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[6]+jitter[1]], gene[D7_WT[1]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[7]+jitter[0]], gene[D7_TKO[0]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)
plt.scatter([x[7]+jitter[1]], gene[D7_TKO[1]], color='white', edgecolors='black', marker='o', s=80, linewidth=0.5)


# Calculate the pvalues from each t-test between groups
# Two-tailed T-test assuming equal variance between groups
D0_pval = ttest_ind(D0_WT_vals, D0_TKO_vals).pvalue
D2_pval = ttest_ind(D2_WT_vals, D2_TKO_vals).pvalue
D4_pval = ttest_ind(D4_WT_vals, D4_TKO_vals).pvalue
D7_pval = ttest_ind(D7_WT_vals, D7_TKO_vals).pvalue



# ticks = [1,2,4,5,7,8,10,11,13,14]
margin_multiplyer = 1
# Add p-value annotation
if D0_pval < 0.05:
    line_margin = max(D0_WT_vals+D0_TKO_vals) + margin_multiplyer
    plt.annotate("*", xy=(1.5, line_margin), xytext=(1.5, line_margin))
    plt.plot([1,2], [line_margin,line_margin], color = 'black', linewidth=0.5)
if D2_pval < 0.05:
    line_margin = max(D2_WT_vals+D2_TKO_vals) + margin_multiplyer
    plt.annotate("*", xy=(4.5, line_margin), xytext=(4.5, line_margin))
    plt.plot([4,5], [line_margin,line_margin], color = 'black', linewidth=0.5)
if D4_pval < 0.05:
    line_margin = max(D4_WT_vals+D4_TKO_vals) + margin_multiplyer
    plt.annotate("*", xy=(7.5, line_margin), xytext=(7.5, line_margin))
    plt.plot([7,8], [line_margin,line_margin], color = 'black', linewidth=0.5)
if D7_pval < 0.05:
    line_margin = max(D7_WT_vals+D7_TKO_vals) + margin_multiplyer             
    plt.annotate("*", xy=(10.5, line_margin), xytext=(10.5, line_margin))
    plt.plot([10,11], [line_margin,line_margin], color = 'black', linewidth=0.5)
if Ezhip_pval < 0.05:
    line_margin = max(D7_TKO_EzhipWT_vals+D7_TKO_EzhipKO_vals) + margin_multiplyer                
    plt.annotate("*", xy=(13.5, line_margin), xytext=(13.5, line_margin))
    plt.plot([13,14], [line_margin,line_margin], color = 'black', linewidth=0.5)

#if Ezhip_D0_pval < 0.05:
#    line_margin = max(D0_WT_vals+D0_TKO_vals+D2_WT_vals+D2_TKO_vals+D4_WT_vals+D4_TKO_vals+D7_WT_vals+D7_TKO_vals+D7_TKO_EzhipWT_vals+D7_TKO_EzhipKO_vals) + margin_multiplyer*2                       
#    plt.annotate("*", xy=(7, line_margin), xytext=(7, line_margin))
#    plt.plot([1,14], [line_margin,line_margin], color = 'black', linewidth=0.5)

    
    
plt.xticks(x, ['D0_WT', 'D0_TKO', 'D2_WT', 'D2_TKO', 'D4_WT', 'D4_TKO', 'D7_WT', 'D7_TKO'])
plt.ylabel('RPKM')
plt.title(gene_name)

# Save p-values to text file
outText = "/Users/jra/Desktop/%s_pvalues.txt" % (gene_name)
with open(outText, "w") as file:
    file.write(gene_name + "\n\n")
    file.write("Two-sided T-test"+ "\n")
    file.write("D0 WT vs TKO: %s \n" % D0_pval)
    file.write("D2 WT vs TKO: %s \n" % D2_pval)
    file.write("D4 WT vs TKO: %s \n" % D4_pval)
    file.write("D7 WT vs TKO: %s \n" % D7_pval)
    file.write("D7 TKO EzhipWT vs EzhipKO: %s \n" % Ezhip_pval)
    file.write("D0 vs EzhipKO: %s \n\n" % Ezhip_D0_pval)
    file.write("Mean RPKM values"+ "\n")
    file.write("D0 mean: WT: %s \tTKO: %s\n" % (D0_WT_mean, D0_TKO_mean))
    file.write("D2 mean: WT: %s \tTKO: %s\n" % (D2_WT_mean, D2_TKO_mean))
    file.write("D4 mean: WT: %s \tTKO: %s\n" % (D4_WT_mean, D4_TKO_mean))
    file.write("D7 mean: WT: %s \tTKO: %s\n" % (D7_WT_mean, D7_TKO_mean))
    file.write("D7 mean: TKO_EzhipWT: %s \tTKO_EzhipKO: %s\n" % (D7_TKO_EzhipWT_mean, D7_TKO_EzhipKO_mean))

outFigure = "/Users/jra/Desktop/%s_bar.pdf" % (gene_name)
plt.savefig(outFigure, dpi=300)


# In[1]:


# Prompts:

# for a dataset with 4 samples and 3 replicates each, use python to plot a bar chart to plot the mean while keeping individual data points
# plot the individual data points in the previous bar chart
# also plot the statistical difference between samples
# add jitter to the scatterplot of individual data points
# redo the previous plot but do not include error bars and plot the p-values on top of the graph and instead of arrows indicate the test comparisons using a line
# instead of plotting the p value, indicate if the p value is below 0.05. Save all p values as a separate text file


# ChatGPT output:

'''
import matplotlib.pyplot as plt
import numpy as np
from scipy.stats import ttest_ind

# Create sample data
sample1 = [4, 5, 3]
sample2 = [5, 6, 2]
sample3 = [7, 6, 7]
sample4 = [2, 1, 3]

# Calculate mean of each sample
sample1_mean = np.mean(sample1)
sample2_mean = np.mean(sample2)
sample3_mean = np.mean(sample3)
sample4_mean = np.mean(sample4)

print(sample1_mean)


# Create a bar chart
x = np.arange(4)
plt.bar(x, [sample1_mean, sample2_mean, sample3_mean, sample4_mean], color='blue', width=0.4)


# Add jitter to the scatter plot of individual data points
jitter = np.random.uniform(-0.1,0.1,size=len(sample1))
for i in range(len(sample1)):
    plt.scatter([x[0]+jitter[i]], [sample1[i]], color='white', edgecolors='black', marker='o', s=80)
jitter = np.random.uniform(-0.1,0.1,size=len(sample2))
for i in range(len(sample2)):
    plt.scatter([x[1]+jitter[i]], [sample2[i]], color='white', edgecolors='black', marker='o', s=80)
jitter = np.random.uniform(-0.1,0.1,size=len(sample3))
for i in range(len(sample3)):
    plt.scatter([x[2]+jitter[i]], [sample3[i]], color='white', edgecolors='black', marker='o', s=80)
jitter = np.random.uniform(-0.1,0.1,size=len(sample4))
for i in range(len(sample4)):
    plt.scatter([x[3]+jitter[i]], [sample4[i]], color='white', edgecolors='black', marker='o', s=80)


# Perform t-test between samples
sample1_2_pvalue = ttest_ind(sample1,sample2).pvalue
sample1_3_pvalue = ttest_ind(sample1,sample3).pvalue
sample1_4_pvalue = ttest_ind(sample1,sample4).pvalue
sample2_3_pvalue = ttest_ind(sample2,sample3).pvalue
sample2_4_pvalue = ttest_ind(sample2,sample4).pvalue
sample3_4_pvalue = ttest_ind(sample3,sample4).pvalue





# Add p-value annotation
if sample1_2_pvalue < 0.05:
    plt.annotate("*", xy=(0.5, max(sample1_mean,sample2_mean)+1), xytext=(1, max(sample1_mean,sample2_mean)+1))
if sample1_3_pvalue < 0.05:
    plt.annotate("*", xy=(1.5, max(sample1_mean,sample3_mean)+1), xytext=(2, max(sample1_mean,sample3_mean)+1))
if sample1_4_pvalue < 0.05:
    plt.annotate("*", xy=(2.5, max(sample1_mean,sample4_mean)+1), xytext=(3, max(sample1_mean,sample4_mean)+1))
if sample2_3_pvalue < 0.05:
    plt.annotate("*", xy=(1.5, max(sample2_mean,sample3_mean)+1), xytext=(2, max(sample2_mean,sample3_mean)+1))
if sample2_4_pvalue < 0.05:
    plt.annotate("*", xy=(2.5, max(sample2_mean,sample4_mean)+1), xytext=(3, max(sample2_mean,sample4_mean)+1))
if sample3_4_pvalue < 0.05:
    plt.annotate("*", xy=(2.5, max(sample3_mean,sample4_mean)+1), xytext=(3, max(sample3_mean,sample4_mean)+1))

#Add lines to indicate test comparison
plt.plot([0.5,1.5], [max(sample1_mean,sample2_mean)+1,max(sample1_mean,sample2_mean)+1], color = 'black')
plt.plot([1.5,2.5], [max(sample1_mean,sample3_mean)+1,max(sample1_mean,sample3_mean)+1], color = 'black')
plt.plot([2.5,3.5], [max(sample1_mean,sample4_mean)+1,max(sample1_mean,sample4_mean)+1], color = 'black')
plt.plot([1.5,2.5], [max(sample2_mean,sample3_mean)+1,max(sample2_mean,sample3_mean)+1], color = 'black')
plt.plot([2.5,3.5], [max(sample3_mean,sample4_mean)+1,max(sample3_mean,sample4_mean)+1], color = 'black')

plt.xticks(x, ['Sample 1', 'Sample 2', 'Sample 3', 'Sample 4'])
plt.ylabel('Value')
plt.title('Individual Replicates')

# Save p-values to text file
with open("/Users/jra/Desktop/pvalues.txt", "w") as file:
    file.write("Sample1-Sample2: %s \n" % sample1_2_pvalue)
    file.write("Sample1-Sample3: {:.10f}\n".format(sample1_3_pvalue))
    file.write("Sample1-Sample4: {:.10f}\n".format(sample1_4_pvalue))
    file.write("Sample2-Sample3: {:.10f}\n".format(sample2_3_pvalue))
    file.write("Sample2-Sample4: {:.10f}\n".format(sample2_4_pvalue))
    file.write("Sample3-Sample4: {:.10f}\n".format(sample3_4_pvalue))

plt.show()
'''


# In[ ]:




