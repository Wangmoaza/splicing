import pandas as pd
import numpy as np

# import data
snv = pd.read_csv('~/DATA/splicing/data/TCGA/TCGA-BRCA/SNV/TCGA-BRCA.mutect2_snv.tsv',
                  sep='\t', header=0, index_col=None)
tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                   sep='\t', header=0, index_col=0)
samples = pd.read_csv('~/DATA/splicing/Analysis/DEG/sample_list/all_438_samples.txt',
                      sep='\t', header=None, index_col=0).index.values
so_term = pd.read_csv('~/DATA/splicing/data/TCGA/TCGA-BRCA/SNV/VEP_description.txt',
                      sep='\t', header=0, index_col=0)

# extract tier 123 genes
snv = snv[snv['gene'].isin(tier.index)]

# filter out panel of normals variants (keep only PASS)
snv = snv[snv['filter'] == "PASS"]

# extract tumor samples
snv['Participant_ID'] = None
snv['tissue_type'] = None
snv.loc[:, ['Participant_ID', 'tissue_type']] = snv['Sample_ID'].str.rsplit('-', 1, expand=True).values
snv = snv[snv['tissue_type'].str.startswith('01')]

# extract all_432 samples (BRCA-event/HRD-high, BRCA-free/HRD-high, BRCA-free/HRD-low)
snv = snv[snv['Participant_ID'].isin(samples)]

# scale variant effect
impact_dic = {'HIGH': 4, "MODERATE": 3, "LOW": 2, "MODIFIER": 1}
so_term['impact_scale'] = so_term['IMPACT'].apply(lambda x: impact_dic[x])
scale_dic = so_term['impact_scale'].to_dict()

# choose representative variant
snv.loc[:, 'impact_scale'] = 0
for i in snv.index:
    try:
        snv.loc[i, 'impact_scale'] = scale_dic[snv.loc[i, 'effect']]
    # if mutliiple variant effect, choose the highest impact
    except KeyError:
        snv.loc[i, 'impact_scale'] = np.max([scale_dic[k] for k in snv.loc[i, 'effect'].split(';')])

# if multiple variants in a gene, choose highest impact variant
remove_list = []
for i, group in snv.groupby(by=['Participant_ID', 'gene']):
    if group.shape[0] > 1:
        remove_list += list(group.index)
        remove_list.remove(group['impact_scale'].idxmax())
snv = snv.drop(index=remove_list)

# reshape dataframe to column: samples, row: genes
snv_pivot = snv.pivot(index="gene", columns="Participant_ID", values="impact_scale")
snv_pivot = snv_pivot.reindex(index=tier.index, columns=samples).fillna(0).astype("int32")

# save dataframes
snv_pivot.to_csv('~/DATA/splicing/Analysis/SNV/all_438.tier123.snv_effect_v2.txt', sep='\t')
snv.to_csv('~/DATA/splicing/Analysis/SNV/all_438.tier123.snv_effect_v2.melted.txt', sep='\t', index=False)
