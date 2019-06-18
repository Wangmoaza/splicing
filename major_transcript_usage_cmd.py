# Haeun Hwangbo 2019-06-18

import argparse
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import *
from statsmodels.stats.multitest import fdrcorrection

def get_args():
    parser = argparse.ArgumentParser()s


def calculate_r(U_arr, n1, n2):
    return 1 - (2 * U_arr) / (n1 * n2)


def calculate_auc(U_arr, n1, n2):
    return U_arr / (n1 * n2)


def manwhitney_test(df, high_group, low_group, alternative='two-sided'):
    stat_list = []
    p_list = []
    for gene in df.index:
        try:
            stat, p = mannwhitneyu(
                df.loc[gene, low_group], df.loc[gene, high_group],
                alternative=alternative)
            stat_list.append(stat)
            p_list.append(p)
        except ValueError:
            stat_list.append(np.nan)
            p_list.append(1)
    U_arr = np.array(stat_list)
    n1, n2 = high_group.shape[0], low_group.shape[0]
    r_arr = calculate_r(U_arr, n1, n2)
    auc_arr = calculate_auc(U_arr, n1, n2)
    return U_arr, r_arr, auc_arr, p_list, fdrcorrection(p_list)[1]


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


def main():

    group_info = pd.read_csv('~/DATA/splicing/Analysis/HR_leafviz/BRCA_free/groups/BRCA_free.HRD_sig3_median_groups.txt',
                             sep='\t', header=None, names=['sample', 'group'],
                             index_col=0, squeeze=True)
    high_group = group_info[group_info == 'high'].index
    low_group = group_info[group_info == 'low'].index

    tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                       header=0, index_col=0, sep='\t')
    tmap = pd.read_csv('/home/haeun/DATA/splicing/Analysis/Quant/TCGA-BRCA/merged/stringtie.BRCA_751.merged.gtf.tmap_extended.tsv',
                       sep='\t', header=0, index_col=4)
    # exclude tier 3 genes
    tmap = tmap[tmap['ref_gene_id'].isin(tier[tier['Tier'] != 3].index)]
    prop = pd.read_csv('/home/haeun/DATA/splicing/Analysis/transcript_usage/TCGA-BRCA/tier12_transcript_proportion.tsv',
                       sep='\t', header=0, index_col=[0, 1])

    ################# diversity #####################
    # number of transcripts with > 1% expression within each gene
    grouped = prop[prop > 0.01].groupby('ref_gene_id')
    grouped.count().to_csv(
        '~/DATA/splicing/Analysis/transcript_usage/TCGA-BRCA/diversity.tsv', sep='\t')
    diversity = grouped.count()

    stat_arr, r_arr, auc_arr, p_arr, adjp_arr = manwhitney_test(
        diversity, high_group, low_group, alternative='less')
    pd.DataFrame({'stat': stat_arr, 'r': r_arr, 'auc': auc_arr, 'p': p_arr, 'adjp': adjp_arr}, index=diversity.index).sort_values('adjp').to_csv(
        '/home/haeun/DATA/splicing/Analysis/transcript_usage/TCGA-BRCA/diversity_tier12.mannwhitney_high_greater_result.tsv', sep='\t')

    col = [sns.husl_palette(10, l =.6)[1], sns.husl_palette(10, l =.6)[7]]
    sns.set(font_scale=3, style='white')
    for i in range(diversity.shape[0]):
        #if adjp_arr[i] >= 0.1 or diversity.index[i] not in tier[tier['Tier'] != 3].index:
        #    continue
        if diversity.index[i] not in ['BRIP1', 'UBE2V2', 'ATM', 'RAD51']:
            continue
        fig, ax = plt.subplots(figsize=(10, 9))
        #ax.set_title('{0} (adjp = {1:.5})'.format(diversity.index[i], adjp_arr[i]))
        ax.set_title(diversity.index[i])

        sns.boxplot(x=group_info.sort_values(), y=diversity.iloc[i, :],
                    ax=ax, palette=col, showfliers=False, width=.5, linewidth=2, notch=True)
        sns.stripplot(x=group_info.sort_values(), y=diversity.iloc[i, :],
                      ax=ax, color='0.2', alpha=0.6)
        ax.set_ylabel('diversity')
        fig.savefig(
            '/home/haeun/DATA/splicing/Analysis/transcript_usage/TCGA-BRCA/figures/{0}_diversity2.png'.format(diversity.index[i]))
        plt.close()


    ################ nonfunctional transcripts #################

    # not protein_coding
    nonfunc = tmap[(tmap['ref_transcript_type'] != 'protein_coding') | (
        tmap['class_code'].isin(['i', 'y', 'p', 's']))].index

    # not APPRIS (annotated transcirpt) or non-annotated transcript
    nonfunc = tmap[(tmap['ref_transcript_code'].isna())
                   | (tmap['class_code'] != "=")].index

    grouped = prop.loc[(slice(None), nonfunc), :].groupby('ref_gene_id')
    result = grouped.sum()

    result.to_csv(
        '/home/haeun/DATA/splicing/Analysis/transcript_usage/Cell-line/notappris_tier12.proportion.tsv', sep='\t')
    stat_arr, r_arr, auc_arr, adjp_arr = manwhitney_test(
        result, high_group, low_group, alternative='less')
    pd.DataFrame({'stat': stat_arr, 'r': r_arr, 'auc': auc_arr, 'adjp': adjp_arr}, index=result.index).sort_values('adjp').to_csv(
        '/home/haeun/DATA/splicing/Analysis/transcript_usage/Cell-line/notappris_tier12.mannwhitney_high_greater_result.tsv', sep='\t')

    for i in range(result.shape[0]):
        if adjp_arr[i] >= 0.05 or diversity.index[i] not in tier[tier['Tier'] != 3].index:
            continue
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.set_title('{0} (adjp = {1:.5})'.format(result.index[i], adjp_arr[i]))
        sns.boxplot(x=group_info.sort_values(), y=result.iloc[i, :],
                    ax=ax, palette=sns.husl_palette(2))
        sns.strippplot(x=group_info.sort_values(), y=result.iloc[i, :],
                       ax=ax, color='black', alpha=0.2)
        fig.savefig(
            '/home/haeun/DATA/splicing/Analysis/transcript_usage/Cell-line/figures/{0}_notappris2.png'.format(result.index[i]))
        plt.close()
