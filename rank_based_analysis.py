# Ha-Eun Hwangbo
# 2019/04/03

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import *
from statsmodels.stats.multitest import fdrcorrection

PATH = "/home/haeun/DATA/splicing/Analysis/heatmap_input/"


SAMPLE_INFO = pd.read_csv(PATH + '285sample/shared_sample_with_HRD.txt',
                          header=None, index_col=0, squeeze=True, sep='\t',
                          names=['sample', 'group'])

TIER_GENES = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                         sep='\t', header=0, index_col=0)

DF = pd.read_csv(PATH + '285sample/me_shared_sample_tier1-3.txt',
                 sep='\t', index_col=0, header=0)
DF = DF[~DF.index.duplicated(keep='first')]


def draw_clustermap(hm_df):
    sns.set()
    col_lut = dict(zip(SAMPLE_INFO.unique(), sns.husl_palette(2)))
    row_lut = dict(zip(TIER_GENES['Tier'].unique(), 'rbg'))
    col_colors = SAMPLE_INFO[hm_df.columns].map(col_lut)
    row_colors = TIER_GENES.loc[hm_df.index, 'Tier'].map(row_lut)
    ax = sns.clustermap(hm_df, figsize=(15, 15), row_cluster=False, col_cluster=False,
                        row_colors=row_colors, col_colors=col_colors)
    return ax


def get_postfix(i):
    if i == 0:
        return "1"
    elif i == 1:
        return "2"
    elif i == 2:
        return "3"
    elif i == 3:
        return "12"
    elif i == 4:
        return "123"


def main():
    tier_1 = TIER_GENES[TIER_GENES['Tier'] == 1].index
    tier_2 = TIER_GENES[TIER_GENES['Tier'] == 2].index
    tier_3 = TIER_GENES[TIER_GENES['Tier'] == 3].index
    tier_12 = TIER_GENES[TIER_GENES['Tier'] != 3].index
    tier_123 = TIER_GENES.index
    gene_subgroup_list = [tier_1, tier_2, tier_3, tier_12, tier_123]
    flag = 'raw'

    for i in range(len(gene_subgroup_list)):
        new_df = DF.reindex(gene_subgroup_list[i])[SAMPLE_INFO.index]
        new_df = new_df.dropna(axis=0)

        # normalize by genes
        if flag == 'zscore':
            normalized_df = pd.DataFrame(
                zscore(new_df, axis=1), index=new_df.index, columns=new_df.columns)
        elif flag == 'rank':
            normalized_df = new_df.rank(axis=1, method='min')
        else:
            normalized_df = new_df

        # get lowest rank / value genes per sample
        result = pd.concat([SAMPLE_INFO, normalized_df.max(axis=0)], axis=1)
        result.columns = ['group', 'lowest_rank']
        result['lowest_genes'] = None
        for sample in result.index:
            new = normalized_df.loc[:, sample]
            result.at[sample, 'lowest_genes'] = new[new
                                                    == result.at[sample, 'lowest_rank']].index.values

        # draw boxplot
        fig, ax = plt.subplots(figsize=(4, 6))
        sns.boxplot(data=result, x='group', y='lowest_rank', ax=ax)
        sns.stripplot(data=result, x='group',
                      y='lowest_rank', alpha=0.5, ax=ax)
        ax.set_title('Tier ' + get_postfix(i))
        plt.tight_layout()
        fig.savefig(
            '/home/haeun/DATA/splicing/Analysis/rank_based/methyl_tier{0}_{1}_box.png'.format(get_postfix(i), flag))
        plt.close()

        # draw heatmap
        hm_df = pd.DataFrame(
            np.nan, index=gene_subgroup_list[i], columns=SAMPLE_INFO.index)
        for k in result.index:
            hm_df.loc[result.loc[k, 'lowest_genes'],
                      k] = result.loc[k, 'lowest_rank']

        draw_clustermap(hm_df).savefig(
            '/home/haeun/DATA/splicing/Analysis/rank_based/methyl_tier{0}_{1}_heatmap.png'.format(get_postfix(i), flag))

        # statistical test
        stat, p = mannwhitneyu(result[result['group'] == 'low']['lowest_rank'],
                               result[result['group'] == 'high']['lowest_rank'])
        # permutation
#        p_values = [p]
#        for _ in range(1000):
#            permuted = pd.Series(np.random.permutation(
#                SAMPLE_INFO.values), index=SAMPLE_INFO.index)
#            high_idx, low_idx = permuted[permuted == 'high'].index, permuted[permuted == 'low'].index
#            stat, p = mannwhitneyu(result.loc[high_idx, 'lowest_rank'],
#                                   result.loc[low_idx, 'lowest_rank'])
#            p_values.append(p)

        print get_postfix(i), p


main()


def sub(x):
    if x < 0 or np.isnan(x):
        return 0
    return x


df_cnv = pd.read_csv(
    PATH + '285sample/CNV_shared_sample_tier1-3.txt', sep='\t', index_col=0, header=0)
df_snv = DF.reindex(df_cnv.index).dropna()[df_cnv.columns]
#aa = (df_cnv - (df_snv == 4).astype(int) <= -2).astype(int)
aa = (df_snv == 4).astype(int) - df_cnv
aa = aa.applymap(sub)
tier_1 = TIER_GENES[TIER_GENES['Tier'] == 1].index
tier_2 = TIER_GENES[TIER_GENES['Tier'] == 2].index
tier_3 = TIER_GENES[TIER_GENES['Tier'] == 3].index
tier_12 = TIER_GENES[TIER_GENES['Tier'] != 3].index
tier_123 = TIER_GENES.index
gene_subgroup_list = [tier_1, tier_2, tier_3, tier_12, tier_123]
for i in range(len(gene_subgroup_list)):
    bb = aa.reindex(gene_subgroup_list[i]).dropna(how='all').astype(int)
    bb.to_csv(
        '/home/haeun/DATA/splicing/Analysis/two-hit/twohit_somatic_nonbinary.tier_{0}.txt'.format(get_postfix(i)), sep='\t')
    #hm_df = bb.replace(0, np.nan)[SAMPLE_INFO.index]
    # draw_clustermap(hm_df).savefig('/home/haeun/DATA/splicing/Analysis/two-hit/twohit_somatic.tier_{0}.heatmap.png'.format(get_postfix(i)))
