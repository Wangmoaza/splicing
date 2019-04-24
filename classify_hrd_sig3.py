import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import pearsonr

df = pd.read_csv("../data/TCGA/TCGA-BRCA.sig3_hrd.txt",
                 sep='\t', header=0, index_col=0)

# BRCA wildtype samples
free = []
with open('../data/TCGA/free_sample.txt', 'r') as f:
    for line in f.readlines():
        free.append(line.strip())

free_df = df.reindex(free).dropna()


sns.pairplot(data=free_df[['Sig3', 'NtAI', 'LST', 'HRD-LOH', 'HRDscore']])
sns.jointplot('Sig3', 'HRDscore', data=df, kind="reg")
fig, ax = plt.subplots(figsize=(5, 5))
sns.scatterplot('Sig3', 'HRDscore', data=df)
ax.axhline(df['HRDscore'].median(), c='black')
ax.axvline(df['Sig3'].median(), c='black')
sns.despine()
fig.savefig("../Analysis/HR_leafviz/HRDscore_sig3_median.png")


high_group = df[(df['HRDscore'] > df['HRDscore'].median())
                & (df['Sig3'] > df['Sig3'].median())].index
low_group = df[(df['HRDscore'] < df['HRDscore'].median())
               & (df['Sig3'] < df['Sig3'].median())].index
new_ser = pd.Series(index=list(high_group) + list(low_group))
new_ser[high_group] = 'high'
new_ser[low_group] = 'low'
new_ser.to_csv('../Analysis/HR_leafviz/HRD_sig3_median_groups.txt',
               sep='\t', header=False)

# median is based on all samples not just BRCA free samples
high_group = free_df[(free_df['HRDscore'] > df['HRDscore'].median()) & (
    free_df['Sig3'] > df['Sig3'].median())].index
low_group = free_df[(free_df['HRDscore'] < df['HRDscore'].median()) & (
    free_df['Sig3'] < df['Sig3'].median())].index
new_ser = pd.Series(index=list(high_group) + list(low_group))
new_ser[high_group] = 'high'
new_ser[low_group] = 'low'
new_ser.to_csv(
    '../Analysis/HR_leafviz/BRCA_free.HRD_sig3_median_groups.txt', sep='\t', header=False)

for i in ['NtAI', 'LST', 'HRD-LOH', 'HRDscore']:
    print pearsonr(free_df[i], free_df['Sig3'])


# 2019/02/27

import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

df = pd.read_csv('~/DATA/splicing/Analysis/Signature/TCGA-BRCA.mSig3_hrd.merged.txt',
                 sep='\t', header=0, index_col=0)

all_group = pd.read_csv('~/DATA/splicing/Analysis/HR_leafviz/All/groups/All.HRD_sig3_median_groups.txt',
                        sep='\t', header=None, index_col=0, squeeze=True)

free_group = pd.read_csv('~/DATA/splicing/Analysis/HR_leafviz/BRCA_free/groups/BRCA_free.HRD_sig3_median_groups.txt',
                         sep='\t', header=None, index_col=0, squeeze=True)

fig, ax = plt.subplots()
sns.scatterplot(data=df, x='sig3_abs', y='HRDscore', ax=ax)
sns.scatterplot(data=df.loc[free_group[free_group == 'high'].index, :], x='sig3_abs', y='HRDscore', ax=ax)
sns.scatterplot(data=df.loc[free_group[free_group == 'low'].index, :], x='sig3_abs', y='HRDscore', ax=ax)
