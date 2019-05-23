# Ha-Eun Hwangbo
# 2019/05/09

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import *
from statsmodels.stats.multitest import fdrcorrection

tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                   header=0, index_col=0, sep='\t')
tmap = pd.read_csv('/home/haeun/DATA/splicing/Analysis/Quant/CCLE-OV/merged/stringtie.CCLE_45.merged.gtf.tmap_extended.tsv',
                   sep='\t', header=0, index_col=4)
tmap = tmap[tmap['ref_gene_id'].isin(tier[tier['Tier'] != 3].index)]


PATH = "/home/haeun/DATA/splicing/Analysis/transcript_usage/CCLE-OV/"

prop = pd.read_csv(PATH + 'tier12_transcript_proportion.tsv',
                   sep='\t', header=0, index_col=[0, 1])

# ccle_name = pd.Series(prop.columns).str.split(".", expand=True)[1].str.replace(
#    '-', '').str.replace('_', '').str.upper() + '_BREAST'
#prop.columns = ccle_name

# use only BRCA-free samples
brca_status = pd.read_csv(
    PATH + 'CCLE-OV.group_info.tsv', sep='\t', header=0, index_col=0)

brca_free = brca_status[brca_status['BRCA_status'] == 'BRCA_free'].index


d2 = pd.read_csv('/home/haeun/DATA/splicing/data/Cell_lines/DepMap/D2_OV_combined_gene_dep_scores.tsv',
                 sep='\t', header=0, index_col=0)

shared_samples = np.intersect1d(
    np.intersect1d(d2.columns, prop.columns), brca_free)
d2 = d2[shared_samples]
prop = prop[shared_samples]


brca_status = brca_status.loc[shared_samples, :]
high_group = brca_status[brca_status['group'] == 'high'].index
low_group = brca_status[brca_status['group'] == 'low'].index
# not protein_coding
nonfunc = tmap[(tmap['ref_transcript_type'] != 'protein_coding') | (
    tmap['class_code'].isin(['i', 'y', 'p', 's']))].index
# not APPRIS (annotated transcirpt) or non-annotated transcript
nonfunc = tmap[(tmap['ref_transcript_code'].isna()) |
               (tmap['class_code'] != "=")].index

# principal = tmap[tmap['ref_transcript_code'].str.contains(
#    'PRINCIPAL') & (tmap['class_code'] == "=")].index

grouped = prop.loc[(slice(None), nonfunc), :].groupby('ref_gene_id')
result = grouped.sum()


def corr_test(a, b, test_gene):
    lists = []
    for gene in test_gene:
        try:
            r, p = pearsonr(a, b.loc[gene, :])
            lists.append([gene, r, p])
        except ValueError:
            lists.append([gene, np.nan, 1])
        except KeyError:
            pass
    df = pd.DataFrame(lists)
    df.columns = ['gene', 'r', 'p']
    df['adjp'] = fdrcorrection(df['p'])[1]
    return df


def mannwhitney_test(a, b, test_gene, alternative="two-sided"):
    lists = []
    for gene in test_gene:
        try:
            r, p = mannwhitneyu(
                a.loc[gene, :], b.loc[gene, :], alternative=alternative)
            lists.append([gene, r, p])
        except ValueError:
            lists.append([gene, np.nan, 1])
        except KeyError:
            pass
    df = pd.DataFrame(lists)
    df.columns = ['gene', 'r', 'p']
    df['adjp'] = fdrcorrection(df['p'])[1]
    return df


def fisher_test(a, b, test_gene, alternative='two-sided'):
    lists = []
    for gene in test_gene:
        # for gene in tier[tier['Tier'] != 3].index:
        try:
            high_dep = a.loc[gene, :].sum()
            low_dep = b.loc[gene, :].sum()
            high_no = high.shape[0] - high_dep
            low_no = low.shape[0] - low_dep
            oddsratio, p = fisher_exact([[high_dep, low_dep],
                                         [high_no, low_no]],
                                        alternative=alternative)
            lists.append([gene, oddsratio, p])
        except KeyError:
            pass

    df = pd.DataFrame(lists)
    df.columns = ['gene', 'or', 'p']
    df['adjp'] = fdrcorrection(df['p'])[1]
    return df


####### correlation test
for i in ['BRIP1', 'ATM', 'RAD51']:
#for i in ['ATM']:
    # for gene in d2.index:
    #df = corr_test(result.loc[i, :], d2, tier[tier['Tier'] != 3].index)
    df = corr_test(result.loc[i, :], d2, d2.index)
    df = df.dropna().sort_values('r', ascending=False)
    df.to_csv(PATH + 'dependency/{0}_notcoding_d2_corr.tsv'.format(i), sep='\t', index=False)
    print i
    #print df[df['adjp'] < 0.1]


df = corr_test(result.loc["BRIP1", :], d2, d2.index)
df = df.dropna().sort_values('r', ascending=False)
df.set_index('gene').reindex(tier[tier['Tier'] != 3].index).dropna()
(df.set_index('gene').reindex(tier[tier['Tier'] != 3].index).dropna()['r'] < 0).sum()
fig, ax = plt.subplots(figsize=(5,5))
sns.regplot(result.loc['ATM', :], d2.loc['POLD1', :], ax=ax)
ax.set_xlabel('ATM_notappris TU')
ax.set_ylabel('POLD D2 score')
fig.savefig(PATH + 'dependency/ATM_notappris_POLD1_relplot.png')

####### Mann-Whitney test
for i in ['BRIP1', 'ATM', 'UBE2V2', 'RAD51']:
    # for gene in d2.index:
    tmp = result.loc[i, :]
    low = tmp[tmp <= tmp.median()].index
    high = tmp[tmp > tmp.median()].index
    df = mannwhitney_test(d2.loc[:, low], d2.loc[:, high],
                          tier[tier['Tier'] != 3].index, alternative='less')
    print i
    print df[df['adjp'] < 0.1]


# dependency threshold -1 true/false
# drop genes where no sample is dependent on it
d2_tf = d2.loc[d2.index[(d2 < -1).sum(axis=1) > 0], :]
d2_tf = d2_tf < -1

####### fisher's exact test
for i in ['BRIP1', 'ATM', 'UBE2V2', 'RAD51']:
    tmp = result.loc[i, :]
    low = tmp[tmp <= tmp.median()].index
    high = tmp[tmp > tmp.median()].index
    df = fisher_test(d2_tf.loc[:, high],
                     d2_tf.loc[:, low],
                     tier[tier['Tier'] != 3].index)
    print i
    print df[df['adjp'] < 0.1]


############ sig3 vs. dependency
hrd_d2 = np.intersect1d(d2.columns, brca_free)

df = corr_test(brca_status.loc[hrd_d2, "sig3_rel"], d2.loc[:, hrd_d2], d2.index)
df = df.dropna().sort_values('r', ascending=False)
df.to_csv(PATH + 'dependency/sig3_rel_d2_corr.tsv', sep='\t', index=False)
############## proportion vs. sig3
shared = np.intersect1d(brca_status.index, result.columns)

target_list = ['BRIP1', 'ATM', 'UBE2V2', 'RAD51']
fig, axes = plt.subplots(2, 2, figsize=(10, 10))
axes = axes.flatten()
for i in range(len(target_list)):
    r, p = pearsonr(result.loc[target_list[i], shared],
                    brca_status.loc[shared, 'sig3_abs'])
    sns.regplot(result.loc[target_list[i], shared],
                brca_status.loc[shared, 'sig3_abs'],
                label="r={0:.3f}, p={1:.3f}".format(r, p), ax=axes[i])
    axes[i].legend()
fig.suptitle('CCLE notappris & sig3 correlation (n=40)')
fig.savefig(PATH + 'notappris_sig3_abs_corr.png')
plt.close()
