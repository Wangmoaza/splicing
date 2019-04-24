import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns

sample_info = pd.read_csv('~/DATA/splicing/Analysis/DEG/sample_list/all_438_samples.txt',
                          sep='\t', header=None, index_col=0, names=['HRD', 'BRCAness'])

protein_coding = pd.read_csv('~/DATA/splicing/data/GENCODE/gencode.v22.protein_coding.id_name_mapping.txt',
                             sep='\t', header=None, index_col=0, squeeze=True)

tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                   sep='\t', header=0, index_col=0)

deg = pd.read_csv('~/DATA/splicing/Analysis/DEG/BRCA_free.HRD_sig3_median_groups.DEG_list.er_covariate.bonferroni.txt',
                  sep='\t', header=0, index_col=0)


df = pd.read_csv('~/DATA/splicing/Analysis/DEG/all_438.log_cpm.txt', sep='\t', header=0, index_col=0)

sample_info.head()
sample_info = sample_info[sample_info['BRCAness'] == "BRCA_free"]
df = df[sample_info.index]
df = df.reindex(deg[np.abs(deg['logFC']) > 0.5].index)
df.shape
# Draw the full plot
sns.clustermap(df, z_score=0, figsize=(13, 13), cmap='vlag')
plt.savefig('/home/haeun/DATA/splicing/Analysis/DEG/clustermap.png')
