# Ha-Eun Hwangbo
# 2019/05/01

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import *
from statsmodels.stats.multitest import fdrcorrection
import os

################# Make dataframe #####################
dic = {}
path = '/home/haeun/DATA/splicing/data/PDX/dependency/'
for filename in os.listdir(path):
    sample = 'PDX_' + filename.rsplit('.', 1)[0].split('-')[1]
    tmp = pd.read_csv(path + filename, sep='\t',
                      header=None, index_col=0, usecols=[0, 2], squeeze=True)
    dic[sample] = tmp.to_dict()

df = pd.DataFrame.from_dict(dic)
df.index.name = 'gene'
df.round(4).to_csv(path + 'PDX_24.dependency.tsv', sep='\t')


#######################################################

prop = pd.read_csv('~/DATA/splicing/Analysis/transcript_usage/notappris_tier12.proportion.tsv',
                    sep='\t', index_col=0, header=0)
dep = pd.read_csv('/home/haeun/DATA/splicing/data/TCGA/TCGA-BRCA/dependency/TCGA-BRCA_1095.dependency.tsv',
                  sep='\t', header=0, index_col=0)
dep = dep.reindex(prop.columns, axis=1)


def get_corr(ser):
    r_list = []
    r, p = pearsonr(ser, )
## correlation
dep.columns
sns.distplot(dep.loc['RAN', :])
sns.distplot(dep.loc['XRCC6', :])
sns.jointplot(data=dep, x='TCGA-4H-AAAK', y='TCGA-Z7-A8R6')
dep.tail(10)

df = pd.DataFrame(np.nan, index=dep.index, columns=prop.index)
df.to_csv('~/DATA/splicing/data/TCGA/TCGA-BRCA/dependency/notappris_tier12.dependency.pearsonr.tsv', sep='\t')
for prop_idx in range(prop.shape[0]):
    print prop_idx
    for dep_idx in range(dep.shape[0]):
        r, _ = pearsonr(prop.iloc[prop_idx, :].values, dep.iloc[dep_idx, :].values)
        df.iat[dep_idx, prop_idx] = r

sns.kdeplot(df['XRCC6'].dropna())

sns.kdeplot(df.loc['BRIP1', prop.index])

sns.regplot(x=prop.loc['BRIP1', :], y=dep.loc['BRIP1', :])
sns.regplot(x=prop.loc['BRIP1', :], y=dep.loc['ATM', :])

sns.scatterplot(dep.loc['ZZZ3', :], dep.loc['XRCC6', :])

dep.corr(method='spearman')
dep.reindex(prop.index).dropna().T.corr()

df.dropna()
