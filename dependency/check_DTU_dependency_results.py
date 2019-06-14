# Haeun Hwangbo 2019/05/15
import pandas as pd
import os
from statsmodels.stats.multitest import fdrcorrection


def isfloat(value):
    try:
        float(value)
        return True
    except ValueError:
        return False


PATH = '/home/haeun/DATA/splicing/Analysis/transcript_usage/PDC/dependency/gsea/'
dir_list = os.listdir(PATH)
new_df = pd.DataFrame(
    columns=['gene', 'minor_type', 'method', 'db', 'Term', 'nes', 'pval', 'fdr'])

for gsea in dir_list:
    if 'repair' not in gsea:
        continue
    if 'tsv' in gsea:
        continue
    if 'ENST' in gsea or 'PDX' in gsea:
        continue
    if 'd2-z' not in gsea:
        continue
    # if 'BRIP1' not in gsea:
    #     continue
    # if "CCLE" in gsea:
    #     continue
    # if "BRIP1_ENST00000259008.2" in gsea:
    #     continue
    try:
        df = pd.read_csv(PATH + gsea + '/gseapy.prerank.gene_sets.report.csv',
                          header=0, index_col=None)
    except IOError:
        continue
    slide = 0
    tokens = gsea.split('_')
    df['gene'] = tokens[0]
    if tokens[1] == "PDX":
        slide += 1
        df['minor_type'] = tokens[1] + '_' + tokens[1 + slide]
    else:
        df['minor_type'] = tokens[1]
    if isfloat(tokens[4 + slide]):
        df['method'] = tokens[3 + slide] + '_' + tokens[4 + slide]
        slide += 1
    else:
        df['method'] = tokens[3 + slide]
    df['db'] = tokens[4 + slide]
    df = df[(df['pval'] < 0.05)]
    new_df = new_df.append(
        df[['gene', 'minor_type', 'method', 'db', 'Term', 'nes', 'pval', 'fdr']])

new_df[new_df['method'] == "corr"].sort_values(['gene', 'minor_type', 'method', 'db', 'nes', 'pval'])
new_df[new_df['method'] == "corr"].sort_values(['gene', 'minor_type', 'method', 'db', 'nes', 'pval']).to_csv(PATH + 'sum_corr_repair_gsea.tsv', sep='\t', index=False)
new_df.sort_values(['gene', 'minor_type', 'method', 'db', 'pval']).to_csv(PATH + 'BRIP1_repair_gsea.tsv', sep='\t', index=False)
