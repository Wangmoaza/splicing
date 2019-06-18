# Haeun Hwangbo 2019-06-18
# Refined from major_transcript_usage.py

import pandas as pd
import numpy as np


STUDY = "CCLE-BRCA"
PATH = "/home/haeun/DATA/splicing/Analysis/transcript_usage/{0}/".format(STUDY)


def make_annotation():
    gencode_version = 29
    appris_version = 38
    if "CCLE" in STUDY:
        gencode_version = 19
        appris_version = 19

    gencode = pd.read_csv('~/DATA/splicing/data/GENCODE/gencode.v{0}.protein_coding.transcripts.tsv'.format(gencode_version),
                          sep='\t', header=0, index_col=1)
#                          names=['gene_id', 'transcript_id', 'gene_type', 'gene_name', 'transcript_type', 'transcript_name'])

    tmap = pd.read_csv('~/DATA/splicing/Analysis/Quant/TCGA-BRCA/merged/stringtie.TCGA-BRCA.merged_strict.gtf.tmap',
                       sep='\t', header=0, index_col=None)

    appris = pd.read_csv('~/DATA/splicing/data/APPRIS/appris_data.principal.hg{0}.txt'.format(appris_version),
                         sep='\t', index_col=2, header=None,
                         names=['gene_name', 'gene_id', 'transcrpt_id', 'ccds_id', 'code'])

    gencode['transcript_id2'] = pd.Series(
        gencode.index).str.split('.', expand=True)[0].values
    tmap['ref_id2'] = tmap['ref_id'].str.split('.', expand=True)[0]

    code_dict = appris['code'].to_dict()
    type_dict = gencode['transcript_type'].to_dict()

    tmap['ref_transcript_type'] = tmap['ref_id'].apply(
        lambda x: type_dict.get(x, np.nan))
    tmap['ref_transcript_code'] = tmap['ref_id2'].apply(
        lambda x: code_dict.get(x, np.nan))
    # drop transcrpits from non-protein coding genes
    tmap = tmap.dropna(subset=['ref_transcript_type'], axis=0)
    tmap = tmap.drop(['TPM', 'FPKM', 'cov'], axis=1)
    tmap.to_csv(
        '~/DATA/splicing/Analysis/Quant/TCGA-BRCA/merged/stringtie.TCGA-BRCA.merged_strict.gtf.tmap_extended.tsv', sep='\t', index=False)
    return tmap


def calc_prop(tmap=None):
    tier = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                       header=0, index_col=0, sep='\t')
    if tmap is None:
        tmap = pd.read_csv('/home/haeun/DATA/splicing/Analysis/Quant/CCLE-BRCA/merged/stringtie.CCLE_56.merged.gtf.tmap_extended.tsv',
                           sep='\t', header=0, index_col=4)
    else:
        tmap = tmap.set_index('qry_id')

    exp = pd.read_csv(PATH + 'strict_merged.txt',
                      sep='\t', header=0, index_col=None)

    # remove gene names from transcripts
    exp = exp.set_index(exp['gene_ENST'].str.split(
        '-', 1, expand=True)[0]).drop('gene_ENST', axis=1)

    # change name to CCLE convention
    ccle_name = pd.Series(exp.columns).str.split(".", expand=True)[1].str.replace(
        '-', '').str.replace('_', '').str.upper() + '_BREAST'
    exp.columns = ccle_name

    # keep transcripts in tmap
    exp = exp.reindex(tmap.index)

    merged = pd.concat([tmap['ref_gene_id'], exp], join='inner', axis=1)
    gene_exp = merged.groupby('ref_gene_id').sum()
    merged = merged.reset_index().set_index(['ref_gene_id', 'qry_id'])
    prop = pd.DataFrame(np.nan, index=merged.index, columns=merged.columns)
    for i in range(merged.shape[0]):
        prop.iloc[i, :] = merged.iloc[i, :] / gene_exp.loc[merged.index[i][0], :]

    prop = prop.sort_index(axis=1)
    prop.to_csv(PATH + 'strict_allgene_transcript_proportion.tsv', sep='\t')
    prop.loc[(tier[tier['Tier'] != 3].index, slice(None)), :].to_csv(
        PATH + 'strict_tier12_transcript_proportion.tsv', sep='\t')

################################################################################


#tmap = make_annotation()
calc_prop()
