# Ha-Eun Hwangbo
# 2019/04/15

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import *
from statsmodels.stats.multitest import fdrcorrection


############ make annottaion file #############
def make_annotation():
    gencode = pd.read_csv('~/DATA/splicing/data/GENCODE/gencode.v29.protein_coding.transcripts.tsv',
                          sep='\t', header=0, index_col=1)
#                          names=['gene_id', 'transcript_id', 'gene_type', 'gene_name', 'transcript_type', 'transcript_name'])

    tmap = pd.read_csv('~/DATA/splicing/Analysis/Quant/merged/stringtie.BRCA_751.merged.gtf.tmap',
                       sep='\t', header=0, index_col=None)

    appris = pd.read_csv('~/DATA/splicing/data/APPRIS/appris_data.principal.txt',
                         sep='\t', index_col=2, header=None,
                         names=['gene_name', 'gene_id', 'transcrpt_id', 'ccds_id', 'code'])

    gencode['transcript_id2'] = pd.Series(
        gencode.index).str.split('.', expand=True)[0].values
    tmap['ref_id2'] = tmap['ref_id'].str.split('.', expand=True)[0]

    code_dict = appris['code'].to_dict()
    type_dict = gencode['transcript_type'].to_dict()

    tmap['ref_transcript_type'] = tmap['ref_id'].apply(
        lambda x: type_dict.get(x, np.nan))
    tmap['ref_transcript_code'] = tmap['ref_id2'].apply(
        lambda x: code_dict.get(x, np.nan))

    tmap = tmap.dropna(subset=['ref_transcript_type'], axis=0)
    tmap = tmap.drop(['TPM', 'FPKM', 'cov'], axis=1)
    tmap.to_csv(
        '~/DATA/splicing/Analysis/Quant/merged/stringtie.BRCA_751.merged.gtf.tmap_extended.tsv', sep='\t', index=False)


################################################################################

make_annotation()
group_info = pd.read_csv('~/DATA/splicing/Analysis/HR_leafviz/BRCA_free/groups/BRCA_free.HRD_sig3_median_groups.txt',
                         sep='\t', header=None, names=['sample', 'group'],
                         index_col=0, squeeze=True)
high_group = group_info[group_info == 'high'].index
low_group = group_info[group_info == 'low'].index


def manwhitney_test(df, alternative='two-sided'):
    stat_list = []
    p_list = []
    for gene in df.index:
        try:
            stat, p = mannwhitneyu(
                df.loc[gene, high_group], df.loc[gene, low_group],
                alternative=alternative)
            stat_list.append(stat)
            p_list.append(p)
        except ValueError:
            stat_list.append(np.nan)
            p_list.append(1)

    return np.array(stat_list), fdrcorrection(p_list)[1]


def ks_test(df):
    stat_list = []
    p_list = []
    for gene in df.index:
        try:
            stat, p = ks_2samp(
                df.loc[gene, high_group], df.loc[gene, low_group])
            p_list.append(p)
            stat_list.append(stat)
        except ValueError:
            stat_list.append(np.nan)
            p_list.append(1)

    return np.array(stat_list), fdrcorrection(p_list)[1]


tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                   header=0, index_col=0, sep='\t')
tmap = pd.read_csv('~/DATA/splicing/Analysis/Quant/merged/stringtie.BRCA_751.merged.gtf.tmap_extended.tsv',
                   sep='\t', header=0, index_col=4)

# change to tpm values
# exclude tier 3 genes
#tmap = tmap[tmap['ref_gene_id'].isin(tier[tier['Tier'] != 3].index)]


exp = pd.read_csv('~/DATA/splicing/Analysis/Quant/ballgown-all/allgene_394.transcript.fpkm.tsv',
                  sep='\t', header=0, index_col=0)
exp = exp.reindex(tmap.index)
merged = pd.concat([tmap['ref_gene_id'], exp], join='inner', axis=1)
gene_exp = merged.groupby('ref_gene_id').sum()

merged = merged.reset_index().set_index(['ref_gene_id', 'qry_id'])

prop = pd.DataFrame(np.nan, index=merged.index, columns=merged.columns)
for i in range(merged.shape[0]):
    prop.iloc[i, :] = merged.iloc[i, :] / gene_exp.loc[merged.index[i][0], :]

prop.shape

prop.to_csv('allgene_transcript_proportion.tsv', sep='\t')
prop.loc[(tier[tier['Tier'] != 3].index, slice(None)), :].to_csv('tier12_transcript_proportion.tsv', sep='\t')
################# diversity #####################
# number of transcripts with > 1% expression within each gene
grouped = prop[prop > 0.01].groupby('ref_gene_id')
#grouped.count().to_csv('~/DATA/splicing/Analysis/transcript_usage/diversity.tsv', sep='\t')
diversity = grouped.count()

stat_arr, adjp_arr = manwhitney_test(diversity, alternative='greater')
pd.DataFrame({'stat': stat_arr, 'adjp': adjp_arr}, index=diversity.index).sort_values('adjp').to_csv('/home/haeun/DATA/splicing/Analysis/transcript_usage/diversity_allgene_test_result.tsv', sep='\t')

for i in range(diversity.shape[0]):
    if adjp_arr[i] >= 0.05:
        continue
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_title('{0} (adjp = {1})'.format(diversity.index[i], adjp_arr[i]))
    sns.boxplot(x=group_info.sort_values(), y=diversity.iloc[i, :],
                ax=ax, palette=sns.husl_palette(2))
    sns.swarmplot(x=group_info.sort_values(), y=diversity.iloc[i, :],
                  ax=ax, color='black', alpha=0.7)
    fig.savefig(
        '/home/haeun/DATA/splicing/Analysis/transcript_usage/{0}_diversity2.png'.format(diversity.index[i]))
    plt.close()


################ nonfunctional transcripts #################

# not protein_coding
nonfunc = tmap[(tmap['ref_transcript_type'] != 'protein_coding') | (
    tmap['class_code'].isin(['i', 'y', 'p', 's']))].index

# not APPRIS (annotated transcirpt) or non-annotated transcript
nonfunc = tmap[(tmap['ref_transcript_code'].isna()) | (tmap['class_code'] != "=")].index
nonfunc.shape

tmap[(tmap['ref_transcript_code'].isna()) | (tmap['class_code'] != "=")]
grouped = prop.loc[(slice(None), nonfunc), :].groupby('ref_gene_id')
result = grouped.sum()

stat_arr, adjp_arr = manwhitney_test(result, alternative='greater')
pd.DataFrame({'stat': stat_arr, 'adjp': adjp_arr}, index=result.index).sort_values('adjp').to_csv('/home/haeun/DATA/splicing/Analysis/transcript_usage/notappris_allgene.mannwhitney_high_greater_result.tsv', sep='\t')

for i in range(result.shape[0]):
    if adjp_arr[i] >= 0.05:
        continue
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.set_title('{0} (adjp = {1:.5})'.format(result.index[i], adjp_arr[i]))
    sns.boxplot(x=group_info.sort_values(), y=result.iloc[i, :],
                ax=ax, palette=sns.husl_palette(2))
    sns.swarmplot(x=group_info.sort_values(), y=result.iloc[i, :],
                  ax=ax, color='black', alpha=0.7)
    fig.savefig(
        '/home/haeun/DATA/splicing/Analysis/transcript_usage/{0}_notcoding.png'.format(result.index[i]))
    plt.close()

pass_arr.shape
tmap[tmap['ref_gene_id'] == 'XRCC6']
prop.loc['XRCC6', :]
grouped.get_group('XRCC6')


############### Do some checking ##################

# log2 gene expression
sns.set()
fig, ax = plt.subplots(figsize=(5, 40))
np.log2(gene_exp + 0.00001).T.plot.box(vert=False, ax=ax)
plt.tight_layout()
fig.savefig(
    '/home/haeun/DATA/splicing/Analysis/transcript_usage/log2_gene_expression.pdf')

# gene expression vs. p-value
new_df = pd.DataFrame({'p': adjp_arr, 'gene_exp': np.log2(
    gene_exp.reindex(result.index).mean(axis=1).values)}, index=result.index)
sns.scatterplot(data=new_df, x='gene_exp', y='p')
