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
tmap = pd.read_csv('~/DATA/splicing/Analysis/Quant/merged/stringtie.BRCA_751.merged.tier_123.gtf.tmap_extended.tsv',
                   sep='\t', header=0, index_col=4)

fpkm = pd.read_csv('~/DATA/splicing/Analysis/Quant/ballgown-all/allgene_394.transcript.fpkm.tsv',
                   sep='\t', header=0, index_col=0)
fpkm = fpkm.reindex(tmap[tmap['ref_gene_id'].isin(
    tier[tier['Tier'] == 1].index)].index)

tpm = pd.read_csv('~/DATA/splicing/Analysis/transcript_usage/novel_merged.txt',
                  sep='\t', header=0, index_col=None)

tpm = tpm.set_index(tpm['gene_ENST'].str.split(
    '-', 1, expand=True)[0]).drop('gene_ENST', axis=1)
tpm = tpm.reindex(tmap[tmap['ref_gene_id'].isin(
    tier[tier['Tier'] == 1].index)].index)
tpm = tpm[fpkm.columns]


fpkm_prop = pd.read_csv('~/DATA/splicing/Analysis/transcript_usage/tier12_transcript_proportion.tsv',
                        sep='\t', header=0, index_col=1)
fpkm_prop = fpkm_prop.drop('ref_gene_id', axis=1)
fpkm_prop = fpkm_prop.reindex(tmap[tmap['ref_gene_id'].isin(
    tier[tier['Tier'] == 1].index)].index)

tpm_prop = pd.read_csv('~/DATA/splicing/Analysis/transcript_usage/novel_transcript_tier_1-3.txt',
                       sep='\t', header=0, index_col=None)
tpm_prop = tpm_prop.set_index(tpm_prop['gene_ENST'].str.split(
    '-', 1, expand=True)[0]).drop('gene_ENST', axis=1)
tpm_prop = tpm_prop.reindex(tmap[tmap['ref_gene_id'].isin(
    tier[tier['Tier'] == 1].index)].index)
tpm_prop = tpm_prop[fpkm_prop.columns]
tpm_prop = tpm_prop.dropna()
fpkm_prop.loc['ENST00000392927.7', :]
fpkm_prop['TCGA-LL-A50Y']

tpm_prop.loc['ENST00000392927.7', :]
for transcript in tpm_prop.index:
    r, p = pearsonr(fpkm_prop.loc[transcript, :], tpm_prop.loc[transcript, :])
    if np.isnan(r):
        print transcript
    else:
        print r
