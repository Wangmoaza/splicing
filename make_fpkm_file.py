import os
import numpy as np
import pandas as pd

PATH = '/home/haeun/DATA/splicing/Analysis/Quant/CCLE-BRCA/fpkms/'
columns = []
mat = []

for name in os.listdir(PATH):
    if "fpkm" in name:
        columns.append(name.split('.')[1].replace('-', '').replace('_', '').upper() + "_BREAST")
        with open(PATH + name, 'r') as f:
            f.readline()
            tmp = [float(i.strip()) for i in f.readlines()]
        mat.append(tmp)

mat = np.array(mat).T

# get transcript names
with open(PATH + 'transcript_names.txt', 'r') as f:
    f.readline()
    idx = [i.strip() for i in f.readlines()]

df = pd.DataFrame(mat, index=idx, columns=columns)
df.to_csv(PATH + 'CCLE-BRCA.allgene_fpkm.tsv', sep='\t')

# make tier12_fpkm
tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                   header=0, index_col=0, sep='\t')
tmap = pd.read_csv('/home/haeun/DATA/splicing/Analysis/Quant/CCLE-BRCA/merged/stringtie.CCLE_56.merged.gtf.tmap_extended.tsv',
                   sep='\t', header=0, index_col=4)
tmap = tmap[tmap['ref_gene_id'].isin(tier[tier['Tier'] != 3].index)]

df.reindex(tmap.index).to_csv(PATH + 'CCLE-BRCA.tier12_fpkm.tsv', sep='\t')
