import pandas as pd
import numpy as np

WORK_DIR = "/home/haeun/DATA/splicing/"

tcga = pd.read_csv(WORK_DIR + 'data/TCGA/TCGA-BRCA/TCGA-BRCA_white_samples.tsv', sep='\t', header=0, index_col=0,
                   usecols=['sample_submitter_id', 'case_submitter_id', 'sample_type_id'])
tcga = tcga[tcga['sample_type_id'] == 1]
tcga = tcga['case_submitter_id'].unique()
rnaseq = []
with open(WORK_DIR + 'Analysis/Quant-all/stringtie.BRCA_751.files.txt', 'r') as f:
    for line in f.readlines():
        rnaseq.append(line.strip().split(
            '/')[-1].split('.')[0].rsplit('-', 4)[0])
somatic = []
with open(WORK_DIR + 'data/TCGA/TCGA-BRCA/mutect2_sample_list.txt', 'r') as f:
    for line in f.readlines():
        somatic.append(line.strip())

germline = []
with open(WORK_DIR + 'data/TCGA/TCGA_all_germline_sample_list.txt', 'r') as f:
    for line in f.readlines():
        germline.append(line.strip().split('\t')[1])

hrd = pd.read_csv(WORK_DIR + 'Analysis/Signature/HRDscore/TCGA.HRD_withSampleID.txt',
                  sep='\t', header=0, index_col=0).index
sig = pd.read_csv(WORK_DIR + 'Analysis/Signature/mSignatureDB/mSigDB_whole.sample.txt',
                  sep='\t', header=0, index_col=None)
sig = sig['sample'].unique()

sqtl = []
with open(WORK_DIR + 'data/leafcutter_result/leafcutter_samples.txt', 'r') as f:
    for line in f.readlines():
        sqtl.append(line.strip())

df = pd.DataFrame(0, index=tcga, columns=[
                  'RNA-seq', 'somatic_mut', 'germline_mut', 'sQTL', 'mut_sig', 'HRDscore'])

for col, arr in zip(df.columns, [rnaseq, somatic, germline, sqtl, sig, hrd]):
    df.loc[np.intersect1d(df.index, np.unique(np.array(arr))), col] = 1

df.index.name = 'participant_ID'
df.to_csv(
    WORK_DIR + 'data/TCGA/TCGA-BRCA/TCGA-BRCA_white.data_availability.txt', sep='\t')

df[df[['RNA-seq', 'somatic_mut', 'mut_sig',
       'HRDscore']].agg(np.all, axis=1)].shape
