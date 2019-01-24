import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pval_qqplot import *

gwas = pd.read_csv('~/2TB_disk/splicing/data/GWAS/OV/28346442-GCST004480-EFO_1001516.h.tsv',
                   sep='\t', header=0, index_col=0, usecols=['variant_id', 'p_value'], na_values=-99)

repair = pd.read_csv('~/2TB_disk/splicing/Analysis/GWAS-enrich/OV_GWAS_SNPs.tier12_overlap.bed',
                     sep='\t', header=None)[13].drop_duplicates().values

gwas = gwas[~gwas.index.duplicated(keep='last')]

fig, ax = plt.subplots(figsize=(5, 5))
pval_qqplot(gwas['p_value'].dropna(), ax=ax, label='Genome-Wide')
pval_qqplot(gwas.reindex(repair).dropna()['p_value'], ax=ax, label="Mapped to Repair Genes")
ax.legend()
fig.savefig('/home/haeun/2TB_disk/splicing/Analysis/GWAS-enrich/OV_qqplot.png')
#fig.savefig('/home/haeun/2TB_disk/splicing/Analysis/GWAS-enrich/BRCA_qqplot.pdf')
plt.show()

gwas[gwas.index.duplicated()]
