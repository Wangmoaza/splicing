# Ha-Eun Hwangbo
# 2019/05/10

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import *
from statsmodels.stats.multitest import fdrcorrection

PATH = "/home/haeun/DATA/splicing/Analysis/transcript_usage/PDC/"

tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                   header=0, index_col=0, sep='\t')
tmap = pd.read_csv('/home/haeun/DATA/splicing/Analysis/Quant/PDC/merged/stringtie.PDX_24.merged.tmap_extended.tsv',
                   sep='\t', header=0, index_col=4)


######################## make prop file ##########################

def calc_prop():
    exp = pd.read_csv(PATH + 'merged.txt',
                      sep='\t', header=0, index_col=None)

    # remove gene names from transcripts
    exp = exp.set_index(exp['gene_ENST'].str.split(
        '-', 1, expand=True)[0]).drop('gene_ENST', axis=1)

    # change name to CCLE convention
    # ccle_name = pd.Series(exp.columns).str.split(".", expand=True)[1].str.replace(
    #     '-', '').str.replace('_', '').str.upper() + '_BREAST'
    # exp.columns = ccle_name

    # keep transcripts in tmap
    exp = exp.reindex(tmap.index)

    merged = pd.concat([tmap['ref_gene_id'], exp], join='inner', axis=1)
    gene_exp = merged.groupby('ref_gene_id').sum()
    merged = merged.reset_index().set_index(['ref_gene_id', 'qry_id'])
    prop = pd.DataFrame(np.nan, index=merged.index, columns=merged.columns)
    for i in range(merged.shape[0]):
        prop.iloc[i, :] = merged.iloc[i, :] / gene_exp.loc[merged.index[i][0], :]

    prop = prop.sort_index(axis=1)
    prop.to_csv(PATH + 'allgene_transcript_proportion.tsv', sep='\t')

calc_prop()
prop = pd.read_csv(PATH + 'allgene_transcript_proportion.tsv', sep='\t', header=0, index_col=[0, 1])
#exclude tier 3 genes
tmap = tmap[tmap['ref_gene_id'].isin(tier[tier['Tier'] != 3].index)]
prop.loc[(slice(None), tmap.index), :].to_csv(PATH + 'tier12_transcript_proportion.tsv', sep='\t')
########## make group_info file
#group_info = pd.read_csv(PATH + 'CCLE_sig3.txt', sep='\t', header=None, names=['Cell', 'sig3_abs', 'sig3_rel'], index_col=0)
#group_info.index = group_info.index + '_OVARY'
#with open(PATH + 'LOH_mut_sample.txt', 'r') as f:
#    brca_list = [i.strip() + '_OVARY' for i in f.readlines()]
#brca_list = brca_list[:-1] # remove OVCAR4
#group_info['BRCA_status'] = "BRCA_free"
#group_info.loc[brca_list, 'BRCA_status'] = 'BRCA_event'
#group_info.index.name = 'Cell_id'
#group_info['group'] = 'low'
#group_info.loc[group_info['sig3_abs'] >= 80, 'group'] = 'high'
#group_info.to_csv(PATH + 'CCLE-OV.group_info.tsv', sep='\t')


prop = pd.read_csv(PATH + 'tier12_transcript_proportion.tsv',
                   sep='\t', header=0, index_col=[0, 1])
group_info = pd.read_csv(PATH + 'CCLE-BRCA.group_info.tsv',
                         sep='\t', header=0, index_col=0)

group_info = group_info[group_info['BRCA_status'] == 'BRCA_free']

group_info = group_info.loc[np.intersect1d(prop.columns, group_info.index), :]
high_group = group_info[group_info['group'] == 'high'].index
low_group = group_info[group_info['group'] == 'low'].index

tmap = tmap[tmap['ref_gene_id'].isin(tier[tier['Tier'] != 3].index)]

target_list = ['BRIP1', 'ATM', 'RAD51', 'UBE2V2']

################# diversity #####################
# number of transcripts with > 1% expression within each gene
grouped = prop[prop > 0.01].groupby('ref_gene_id')
diversity = grouped.count()
diversity.to_csv(PATH + 'tier12_diversity.tsv', sep='\t')
merged = pd.concat([diversity.T, group_info['group']], axis=1)

merged.head()
sns.set(font_scale=3)
sns.set_style('white')
for gene in target_list:
    stat, p = mannwhitneyu(merged[merged['group'] == 'low'][gene],
                           merged[merged['group'] == 'high'][gene],
                           alternative='less')
    fig, ax = plt.subplots(figsize=(10, 9))
    ax.set_title('{0} Diversity (p = {1:.5})'.format(gene, p))
    sns.boxplot(x='group', y=gene, data=merged,
                ax=ax, palette=sns.husl_palette(2),
                showfliers=False, width=.5, linewidth=2, notch=False)
    sns.stripplot(x="group", y=gene, data=merged, color=".2", alpha=".6")
    fig.savefig(PATH + 'figures/{0}_diversity.png'.format(gene))
    # plt.close()


################ nonfunctional transcripts #################

# not protein_coding
nonfunc = tmap[(tmap['ref_transcript_type'] != 'protein_coding') | (
    tmap['class_code'].isin(['i', 'y', 'p', 's']))].index

# not APPRIS (annotated transcirpt) or non-annotated transcript
nonfunc = tmap[(tmap['ref_transcript_code'].isna())
               | (tmap['class_code'] != "=")].index

grouped = prop.loc[(slice(None), nonfunc), :].groupby('ref_gene_id')
result = grouped.sum()
merged = pd.concat([result.T, group_info['group']], axis=1)

for gene in target_list:
    stat, p = mannwhitneyu(merged[merged['group'] == 'low'][gene],
                           merged[merged['group'] == 'high'][gene],
                           alternative='less')
    fig, ax = plt.subplots(figsize=(10, 9))
    ax.set_title('{0} Minor (notcoding) (p = {1:.5})'.format(gene, p))
    sns.boxplot(x='group', y=gene, data=merged,
                ax=ax, palette=sns.husl_palette(2),
                showfliers=False, width=.5, linewidth=2, notch=False)
    sns.stripplot(x="group", y=gene, data=merged, color=".2", alpha=".6")
    fig.savefig(PATH + 'figures/{0}_notcoding.png'.format(gene))
    # plt.close()



############### Do some checking ##################

# log2 gene expression
sns.set()
fig, ax = plt.subplots(figsize=(5, 40))
np.log2(gene_exp + 0.00001).T.plot.box(vert=False, ax=ax)
plt.tight_layout()
fig.savefig(
    '/home/haeun/DATA/splicing/Analysis/transcript_usage/CCLE-OV/log2_gene_expression.pdf')

# gene expression vs. p-value
new_df = pd.DataFrame({'p': adjp_arr, 'gene_exp': np.log2(
    gene_exp.reindex(result.index).mean(axis=1).values)}, index=result.index)
sns.scatterplot(data=new_df, x='gene_exp', y='p')
