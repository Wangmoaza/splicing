# Ha-Eun Hwangbo
# 2019/04/24

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import *

tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                   header=0, index_col=0, sep='\t')
tmap = pd.read_csv('~/DATA/splicing/Analysis/Quant/merged/stringtie.BRCA_751.merged.gtf.tmap_extended.tsv',
                   sep='\t', header=0, index_col=4)

det = pd.read_csv('~/DATA/splicing/Analysis/DET/allgene_DET_results.tsv',
                  sep='\t', header=0, index_col=2)

det = pd.concat([tmap['ref_gene_id'], det], axis=1,
                join='inner').drop('geneNames', axis=1)

det.to_csv('~/DATA/splicing/Analysis/DET/allgene_DET_results.tsv', sep='\t')

det[det['ref_gene_id'] == 'RAD51']
det[det['ref_gene_id'].isin(tier[tier['Tier'] != 3].index)].to_csv(
    '~/DATA/splicing/Analysis/DET/allgene_DET_results_tier12.tsv', sep='\t')

exp = pd.read_csv('~/DATA/splicing/Analysis/transcript_usage/novel_merged.txt',
                  sep='\t', header=0, index_col=None)
exp = exp.set_index(exp['gene_ENST'].str.split(
    '-', 1, expand=True)[0]).drop('gene_ENST', axis=1)
exp = exp.reindex(tmap.index)
merged = pd.concat([tmap['ref_gene_id'], exp], join='inner', axis=1)
