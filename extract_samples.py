import pandas as pd
import numpy as np


htseq = pd.read_csv('~/DATA/splicing/data/TCGA/TCGA-BRCA/gene-htseq/TCGA-BRCA.htseq_counts.tsv',
                    sep='\t', header=0, index_col=0)
sample_info = pd.read_csv('~/DATA/splicing/Analysis/DEG/sample_list/all_438_samples.txt',
                          sep='\t', header=None, index_col=0)
protein_coding = pd.read_csv('~/DATA/splicing/data/GENCODE/gencode.v22.protein_coding.id_name_mapping.txt',
                             sep='\t', header=None, index_col=0, squeeze=True)

# unwrap log2
htseq = np.power(2, htseq) - 1

# extract tumor samples
our_samples = []
with open('/home/haeun/DATA/splicing/Analysis/Quant-all/htseq_tier123_sample_list.txt', 'r') as f:
    for line in f.readlines():
        our_samples.append(line.rsplit('-', 3)[0])

htseq = htseq.loc[:, our_samples]
htseq.columns = pd.Series(htseq.columns).str.rsplit(
    '-', 1, expand=True)[0].values
# extract samples in sample_info
htseq = htseq.loc[:, sample_info.index]

# filter genes
filter = htseq[htseq.agg(sum, axis=1) > 100].index
htseq = htseq.loc[np.intersect1d(filter, protein_coding.index), :]
htseq.to_csv(
    '~/DATA/splicing/Analysis/DEG/all_438.HRD_sig3_median_groups.protein_coding.htseq.txt', sep='\t')
