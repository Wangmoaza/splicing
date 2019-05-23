# Ha-Eun Hwangbo
# 2019/04/24

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import *

group_info = pd.read_csv('~/DATA/splicing/Analysis/HR_leafviz/BRCA_free/groups/BRCA_free.HRD_sig3_median_groups.txt',
                         sep='\t', header=None, names=['sample', 'group'],
                         index_col=0, squeeze=True)
high_group = group_info[group_info == 'high'].index
low_group = group_info[group_info == 'low'].index
tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                   header=0, index_col=0, sep='\t')
tmap = pd.read_csv('~/DATA/splicing/Analysis/Quant/merged/stringtie.BRCA_751.merged.gtf.tmap_extended.tsv',
                   sep='\t', header=0, index_col=4)

hrd = pd.read_csv('~/DATA/splicing/Analysis/Signature/TCGA-BRCA.mSig3_hrd.merged.txt',
                  sep='\t', header=0, index_col=0)
hrd = hrd.reindex(group_info.index)
# change to exp values
# exclude tier 3 genes
tmap = tmap[tmap['ref_gene_id'].isin(tier[tier['Tier'] != 3].index)]

exp = pd.read_csv('~/DATA/splicing/Analysis/transcript_usage/novel_merged.txt',
                  sep='\t', header=0, index_col=None)
exp = exp.set_index(exp['gene_ENST'].str.split(
    '-', 1, expand=True)[0]).drop('gene_ENST', axis=1)
exp = exp.reindex(tmap.index)
merged = pd.concat([tmap['ref_gene_id'], exp], join='inner', axis=1)
gene_exp = merged.groupby('ref_gene_id').sum()

merged = merged.reset_index().set_index(['ref_gene_id', 'qry_id'])

prop = pd.DataFrame(np.nan, index=merged.index, columns=merged.columns)
for i in range(merged.shape[0]):
    prop.iloc[i, :] = merged.iloc[i, :] / gene_exp.loc[merged.index[i][0], :]

prop = prop.reindex(group_info.index, axis=1)

# not protein_coding
nonfunc = tmap[(tmap['ref_transcript_type'] != 'protein_coding') | (
    tmap['class_code'].isin(['i', 'y', 'p', 's']))].index

# not APPRIS (annotated transcirpt) or non-annotated transcript
nonfunc = tmap[(tmap['ref_transcript_code'].isna()) |
               (tmap['class_code'] != "=")].index

grouped = prop.loc[(slice(None), nonfunc), :].groupby('ref_gene_id')
result = grouped.sum()

df = pd.concat([result.T, hrd], axis=1)
df = df.sort_values('group')
sns.set()
for gene in result.index:
    r1, p1 = pearsonr(df[gene], df['HRDscore'])
    r2, p2 = pearsonr(df[gene], df['sig3_rel'])
    if p1 < 0.05 or p2 < 0.05:
        fig, axes = plt.subplots(1, 2, figsize=(12, 6))
        sns.regplot(data=df, x=gene, y='HRDscore', scatter_kws={"alpha": 0.1}, color='black', ax=axes[0])
        sns.regplot(data=df, x=gene, y='sig3_rel', scatter_kws={"alpha": 0.1}, color='black', ax=axes[1])
        sns.scatterplot(data=df, x=gene, y='HRDscore', hue='group', palette=sns.husl_palette(2), ax=axes[0])
        sns.scatterplot(data=df, x=gene, y='sig3_rel', hue='group', palette=sns.husl_palette(2), ax=axes[1])
        axes[0].set_title('r = {0:.5f}, p = {1:.5f}'.format(r1, p1))
        axes[1].set_title('r = {0:.5f}, p = {1:.5f}'.format(r2, p2))
        fig.suptitle(gene + ' proportion of notappris')
        plt.tight_layout()
        fig.savefig(
            '/home/haeun/DATA/splicing/Analysis/transcript_usage/{0}_notappris.hrd_correlation.png'.format(gene))
        plt.close()
