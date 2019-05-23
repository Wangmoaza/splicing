# Haeun Hwangbo 2019/05/15
import pandas as pd
import os

PATH = '/home/haeun/DATA/splicing/Analysis/transcript_usage/CCLE-OV/dependency/gsea/'
dir_list = os.listdir(PATH)
new_df = pd.DataFrame(
    columns=['gene', 'minor_type', 'db', 'Term', 'nes', 'pval', 'fdr'])
for gsea in dir_list:
    if 'repair' not in gsea:
        continue
    if 'tsv' in gsea:
        continue
    if 'sig3' not in gsea:
        continue
    df = pd.read_csv(PATH + gsea + '/gseapy.prerank.gene_sets.report.csv',
                     header=0, index_col=None)
    #print '****************************'
    #print gsea
    #a = df.index.str.contains('DNA')
    #b = df.index.str.contains('REPAIR')
    #c = df.index.str.contains('CHECKPOINT')
    #d = df.index.str.contains('DAMAGE')
    #e = df.index.str.contains('HOMOLOGOUS')
    #condition = a | b | c | d | e
    tokens = gsea.split('_')
    df['gene'] = tokens[0]
    df['minor_type'] = tokens[1]
    df['db'] = tokens[4]
    df = df[(df['pval'] < 0.05)]
    new_df = new_df.append(
        df[['gene', 'db', 'minor_type', 'Term', 'nes', 'pval', 'fdr']])

new_df = new_df[['gene', 'minor_type', 'db', 'Term', 'nes', 'pval', 'fdr']]
new_df
new_df.sort_values(['gene', 'minor_type', 'db']).to_csv(PATH + 'repair_sig_gsea.tsv',
                                                        sep='\t', index=False)

new_df.sort_values('gene')
