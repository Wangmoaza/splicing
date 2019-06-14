# Ha-Eun Hwangbo
# 2019/05/09

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import *
from statsmodels.stats.multitest import fdrcorrection

TUMOR = "BRCA"
PATH = "/home/haeun/DATA/splicing/Analysis/transcript_usage/CCLE-{0}/".format(
    TUMOR)

tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                   header=0, index_col=0, sep='\t')
tmap = pd.read_csv('/home/haeun/DATA/splicing/Analysis/Quant/CCLE-{0}/merged/stringtie.CCLE_56.merged.gtf.tmap_extended.tsv'.format(TUMOR),
                   sep='\t', header=0, index_col=4)
tmap = tmap[tmap['ref_gene_id'].isin(tier[tier['Tier'] != 3].index)]

prop = pd.read_csv(PATH + 'tier12_transcript_proportion.tsv',
                   sep='\t', header=0, index_col=[0, 1])

# use only BRCA-free samples
brca_status = pd.read_csv(PATH + 'CCLE-{0}.group_info.tsv'.format(TUMOR),
                          sep='\t', header=0, index_col=0)

brca_free = brca_status[brca_status['BRCA_status'] == 'BRCA_free'].index

d2 = pd.read_csv('/home/haeun/DATA/splicing/data/Cell_lines/DepMap/D2_{0}_combined_gene_dep_scores_minus.tsv'.format(TUMOR),
                 sep='\t', header=0, index_col=0)


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


def median_test(a, b, test_gene, group='median'):
    # a : transcript usage (series)
    # b : dependency (dataframe)
    if isinstance(group, tuple):
        lower_threshold = group[0]
        upper_threshold = group[1]
    elif isinstance(group, float):
        lower_threshold = group
        upper_threshold = group
    elif isinstance(group, str):
        if group == 'median':
            lower_threshold = a.median()
            upper_threshold = a.median()
        elif group == 'mean':
            lower_threshold = a.mean()
            upper_threshold = a.mean()
        elif group == 'quantile':
            lower_threshold = a.quantile(0.25)
            upper_threshold = a.quantile(0.75)
        else:
            print "Wrong group argument"
            return
    under_group = a[a < lower_threshold].index
    over_group = a[a >= upper_threshold].index
    lists = []
    for gene in test_gene:
        try:
            r = b.loc[gene, over_group].median(
            ) - b.loc[gene, under_group].median()
            lists.append([gene, r])
        except ValueError:
            lists.append([gene, np.nan])
        except KeyError:
            pass
    df = pd.DataFrame(lists)
    df.columns = ['gene', 'med_diff']
    return df


def transcript_type(t_type):
    """Returns transcript names of given type."""
    # not protein_coding
    if t_type == "notcoding":
        nonfunc = tmap[(tmap['ref_transcript_type'] != 'protein_coding') | (
            tmap['class_code'].isin(['i', 'y', 'p', 's']))].index
    # not APPRIS (annotated transcirpt) or non-annotated transcript
    elif t_type == "notappris":
        nonfunc = tmap[(tmap['ref_transcript_code'].isna())
                       | (tmap['class_code'] != "=")].index
    elif t_type == "principal":
        nonfunc = tmap[tmap['ref_transcript_code'].str.contains(
            'PRINCIPAL') & (tmap['class_code'] == "=")].index
    elif t_type == 'minor':
        nonfunc = tmap[~(tmap['ref_transcript_code'].str.contains(
            'PRINCIPAL') & (tmap['class_code'] == "="))].index
    else:
        print "try again"

    return nonfunc


shared_samples = np.intersect1d(
    np.intersect1d(d2.columns, prop.columns), brca_free)
d2 = d2[shared_samples]
prop = prop[shared_samples]
brca_status = brca_status.loc[shared_samples, :]

target_list = ['BRIP1', 'ATM', 'RAD51']
tmap.head()

t_type = 'minor'
nonfunc = transcript_type(t_type)
grouped = prop.loc[(slice(None), nonfunc), :].groupby('ref_gene_id')
result = grouped.sum()


# for minor sum (notappris / notcoding)
# correlation based
for i in target_list:
    # for i in ['ATM']:
    # for gene in d2.index:
    #df = corr_test(result.loc[i, :], d2, tier[tier['Tier'] != 3].index)
    df = corr_test(result.loc[i, :], d2, d2.index)
    df = df.dropna().sort_values('r', ascending=False)
    df.to_csv(
        PATH + 'dependency/{0}_{1}_d2_corr.tsv'.format(i, t_type), sep='\t', index=False)
    print i

# median difference based
for gene in target_list:
    for criteria in ['mean', 'median', 'quantile', (0.1, 0.3)]:
        df = median_test(result.loc[gene, :], d2, d2.index, group=criteria)
        df = df.dropna().sort_values('med_diff', ascending=False)
        if isinstance(criteria, tuple):
            df.to_csv(PATH + 'dependency/{0}_{1}_d2_{2}_{3}.tsv'.format(gene, t_type, criteria[0], criteria[1]),
                      sep='\t', index=False)
        else:
            df.to_csv(PATH + 'dependency/{0}_{1}_d2_{2}.tsv'.format(gene, t_type, criteria),
                      sep='\t', index=False)

# for individual transcript
# keep transcripts with > 5% mean proportion
result = prop.dropna(axis=0)[prop.dropna(axis=0).mean(axis=1) > 0.05]

for gene in target_list:
    for transcript in result.loc[gene, :].index:
        df = corr_test(result.loc[(gene, transcript), :], d2, d2.index)
        df = df.dropna().sort_values('r', ascending=False)
        df.to_csv(PATH + 'dependency/{0}_{1}_d2_corr.tsv'.format(gene, transcript),
                  sep='\t', index=False)
for gene in target_list:
    for transcript in result.loc[gene, :].index:
        for criteria in ['mean', 'median', 'quantile', (0.1, 0.3)]:
            df = median_test(
                result.loc[(gene, transcript), :], d2, d2.index, group=criteria)
            df = df.dropna().sort_values('med_diff', ascending=False)
            if isinstance(criteria, tuple):
                df.to_csv(PATH + 'dependency/{0}_{1}_d2_{2}_{3}.tsv'.format(gene, transcript, criteria[0], criteria[1]),
                          sep='\t', index=False)
            else:
                df.to_csv(PATH + 'dependency/{0}_{1}_d2_{2}.tsv'.format(gene, transcript, criteria),
                          sep='\t', index=False)


# target_repair = tier[tier['Pathway'].str.contains('Nucleotide excision')
target_repair = tier[tier['Tier'] != 3]
target_gene = ('BRIP1', "ENST00000577598.1")
df = corr_test(result.loc[target_gene, :], d2, target_repair.index)
df = df.dropna().sort_values('r', ascending=False)
df = df.set_index('gene').dropna().sort_values('r')
df[df['r'] > 0]
d2_gene = 'ERCC4'
for d2_gene in df[(df['p'] < 0.05) & (df['r'] > 0)].index:
    fig, ax = plt.subplots(figsize=(5, 5))
    sns.regplot(result.loc[target_gene, :], d2.loc[d2_gene, :], ax=ax,
                label='r={0:.3f} p={1:.3f}'.format(df.loc[d2_gene, 'r'], df.loc[d2_gene, 'p']))
    ax.set_xlabel('{0} {1} TU'.format(target_gene, t_type))
    ax.set_ylabel('{0} -D2 score'.format(d2_gene))
    ax.legend()
    fig.tight_layout()
    fig.savefig(
        PATH + 'dependency/{0}_{1}_{2}_relplot.png'.format(target_gene[0], target_gene[1], d2_gene))

# sig3 vs. dependency
hrd_d2 = np.intersect1d(d2.columns, brca_free)

df = corr_test(brca_status.loc[hrd_d2, "sig3_rel"],
               d2.loc[:, hrd_d2], d2.index)
df = df.dropna().sort_values('r', ascending=False)
df.to_csv(PATH + 'dependency/sig3_rel_d2_corr.tsv', sep='\t', index=False)


# proportion vs. sig3
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
fig.suptitle('CCLE {0} & sig3 correlation'.format(t_type))
fig.savefig(PATH + '{0}_sig3_abs_corr.png'.format(t_type))
plt.close()
