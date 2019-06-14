# Haeun Hwangbo 2019-06-12

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from scipy.stats import *
from statsmodels.stats.multitest import fdrcorrection


def corr_test(a, b, test_gene):
    lists = []
    for gene in test_gene:
        try:
            r, p = spearmanr(a, b.loc[gene, :])
            lists.append([gene, r, p])
        except ValueError:
            lists.append([gene, np.nan, 1])
        except KeyError:
            pass
    df = pd.DataFrame(lists)
    df.columns = ['gene', 'r', 'p']
    df['adjp'] = fdrcorrection(df['p'])[1]
    return df


TUMOR = "PDC"
PATH = "/home/haeun/DATA/splicing/Analysis/transcript_usage/{0}/".format(TUMOR)

tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                   header=0, index_col=0, sep='\t')
tmap = pd.read_csv('/home/haeun/DATA/splicing/Analysis/Quant/{0}/merged/stringtie.PDX_24.merged.tmap_extended.tsv'.format(TUMOR),
                   sep='\t', header=0, index_col=4)
tmap = tmap[tmap['ref_gene_id'].isin(tier[tier['Tier'] != 3].index)]

prop = pd.read_csv(PATH + 'tier12_transcript_proportion.tsv',
                   sep='\t', header=0, index_col=[0, 1])

exp = pd.read_csv('/home/haeun/DATA/splicing/Analysis/DEG/PDC/PDX.log_cpm.txt',
                  sep='\t', header=0, index_col=0)

exp = exp[prop.columns]
exp.head()
target_repair = tier[tier['Pathway'].str.contains('Nucleotide excision')]
#target_repair = tier[tier['Tier'] != 3]
#target_gene = "BRIP1"
target_gene = ('BRIP1', "ENST00000577598.5")

df = corr_test(prop.loc[target_gene, :], exp, target_repair.index)
df = df.dropna().set_index('gene').sort_values('r', ascending=False)

df
df.to_csv(PATH + 'expression/BRIP_ENST00000577598.5_NER_exp_corr.tsv', sep='\t')

for d2_gene in target_repair.index:
    try:
        fig, ax = plt.subplots(figsize=(5, 5))
        sns.regplot(prop.loc[target_gene, :], exp.loc[d2_gene, :], ax=ax,
                    label='r={0:.3f} p={1:.3f}'.format(df.loc[d2_gene, 'r'], df.loc[d2_gene, 'p']))
        ax.set_xlabel('{0} TU'.format(target_gene))
        ax.set_ylabel('{0} logCPM'.format(d2_gene))
        ax.legend()
        fig.tight_layout()
        # fig.savefig(
        #     PATH + 'dependency/{0}_{1}_{2}_z_relplot.png'.format(target_gene, t_type, d2_gene))
        fig.savefig(
           PATH + 'expression/figures/{0}_{1}_{2}_exp_relplot.png'.format(target_gene[0], target_gene[1], d2_gene))
        plt.close()
    except KeyError:
        print d2_gene
