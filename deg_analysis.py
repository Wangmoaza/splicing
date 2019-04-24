import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

# mds Analysis
mds = pd.read_csv('/home/haeun/DATA/splicing/Analysis/DEG/MDSplot_coordinates.txt',
                  sep='\t', header=0, index_col=0)

sample_info = pd.read_csv('~/DATA/splicing/Analysis/HR_leafviz/BRCA_free/groups/BRCA_free.HRD_sig3_median_groups.txt',
                          sep='\t', header=None, index_col=0, squeeze=True)

clinical_info = pd.read_csv('/home/haeun/DATA/splicing/Analysis/DEG/BRCA_free.HRD_sig3_median_groups.clinical_info.txt',
                           sep='\t', header=0, index_col=0)

df = pd.concat([mds, sample_info, clinical_info], axis=1)
df.columns = ['mds1', 'mds2', 'HRD', 'ER', 'PR', 'HER2', 'stage', 'age']

fig, axes = plt.subplots(2, 3, figsize=(21, 14))
axes = axes.flatten()
sns.scatterplot(data=df, x="mds1", y="mds2", hue="HRD", ax=axes[0])
for i in range(1, 6):
    sns.scatterplot(data=df, x="mds1", y="mds2", style="HRD", hue=df.columns[i + 2], ax=axes[i])
    axes[i].set_title(df.columns[i + 2])
fig.tight_layout()
fig.savefig('/home/haeun/DATA/splicing/Analysis/DEG/MDS_plot.clustering.png')
fig.savefig('/home/haeun/DATA/splicing/Analysis/DEG/MDS_plot.clustering.pdf')
plt.close()
df.columns
