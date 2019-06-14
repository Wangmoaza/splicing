# Haeun Hwangbo 2019-06-05

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import *
from statsmodels.stats.multitest import fdrcorrection


TUMOR = "BRCA"
PATH = "/home/haeun/DATA/splicing/Analysis/transcript_usage/TCGA-{0}/".format(TUMOR)

tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                   header=0, index_col=0, sep='\t')
tmap = pd.read_csv('/home/haeun/DATA/splicing/Analysis/Quant/TCGA-{0}/merged/stringtie.BRCA_751.merged.gtf.tmap_extended.tsv'.format(TUMOR),
                   sep='\t', header=0, index_col=4)
tmap = tmap[tmap['ref_gene_id'].isin(tier[tier['Tier'] != 3].index)]

prop = pd.read_csv(PATH + 'tier12_transcript_proportion.tsv',
                   sep='\t', header=0, index_col=[0, 1])

# use only BRCA-free samples
group_info = pd.read_csv('~/DATA/splicing/Analysis/DEG/sample_list/all_438_samples.txt',
                         sep='\t', header=None, index_col=0,
                         names=['sample', 'group', 'BRCA_status'])

sig9 = pd.read_csv('~/DATA/splicing/Analysis/Signature/mSigDB.sig9.txt', sep='\t', header=None, index_col=1,
                   names=['sig', 'sample', 'sig9_rel', 'stage'])

mut = pd.read_csv('~/DATA/splicing/Analysis/Signature/MutationalPatterns/TCGA-BRCA.96subs_matrix.txt',
                  sep='\t', header=0, index_col=0).sum()
sig9['sig9_abs'] = sig9['sig9_rel'] * mut.reindex(sig9.index)

group_info = pd.concat([group_info, sig9[['sig9_rel', 'sig9_abs']]], axis=1, join='inner')

#prop.loc['BRIP1']
#fig, ax = plt.subplots(figsize=(10, 6))
#                prop.loc[('BRIP1', 'ENST00000577598.5'), :], ax=ax)
#fig.savefig(PATH + 'BRIP1_major_vs_minor_prop.png')
#fig, ax = plt.subplots(figsize=(10, 6))
#sns.distplot(prop.loc[('BRIP1', 'ENST00000259008.6'), :])
#fig.savefig(PATH + 'BRIP1-201_dist.png')
#fig, ax = plt.subplots(figsize=(10, 6))
#sns.distplot(prop.loc[('BRIP1', 'ENST00000577598.5'), :])
#fig.savefig(PATH + 'BRIP1-202_dist.png')

gene = "RAD51"
transcript = 'ENST00000525066.5'
shared_samples = np.intersect1d(group_info[group_info['sig9_rel'] > 0].index, prop.columns)
fig, ax = plt.subplots(figsize=(5, 5))
r, p = pearsonr(prop.loc[(gene, transcript), shared_samples],
                group_info.loc[shared_samples, 'sig9_abs'])
sns.regplot(prop.loc[(gene, transcript), shared_samples],
            group_info.loc[shared_samples, 'sig9_abs'], ax=ax,
            label='r={0:.3f} p={1:.3f}'.format(r, p))
ax.legend()
fig.savefig(PATH + 'sig9/{0}-{1}_sig9_nonzero_regplot.png'.format(gene, transcript))
