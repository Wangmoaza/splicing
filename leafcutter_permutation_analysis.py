import pandas as pd
import numpy as np
import scipy.stats as stats
import matplotlib.pyplot as plt
import seaborn as sns


def import_data():
    df = pd.read_csv('~/DATA/splicing/Analysis/HR_leafviz/BRCA_free/BRCA_free_cluster_significance.permute_adjusted_388.txt',
                     sep='\t', header=0)

    tier_genes = pd.read_csv('~/DATA/splicing/data/KEGG/KEGG_plus_curated_genes_tier_ensembl.txt',
                             sep='\t', header=0, index_col=0)

    df = df.dropna(subset=['genes'])
    #df = df.set_index(['genes', 'cluster'])
    df.loc[:, 'tier'] = 0
    df = df[~df['genes'].str.contains(',')]
    df = df.set_index('cluster')
    for i in tier_genes.index:
        df.loc[df[df['genes'] == i].index, 'tier'] = tier_genes.loc[i, 'Tier']

    effect_size = pd.read_csv('~/DATA/splicing/Analysis/HR_leafviz/BRCA_free/results/BRCA_free_effect_sizes.txt',
                              sep='\t', header=0)

    effect_size[['chr', 'start', 'end', 'cluster']] = effect_size['intron'].str.split(':', expand=True)
    effect_size.loc[:, 'cluster'] = effect_size['chr'] + ":" + effect_size['cluster']
    grouped = effect_size.groupby('cluster')
    max_deltapsi = grouped['deltapsi'].apply(lambda x: np.max(np.abs(x)))
    cluster_start = grouped['start'].min()
    cluster_end = grouped['end'].max()
    merged = pd.concat([df, max_deltapsi, cluster_start, cluster_end], axis=1, join='inner')
    return merged


df = import_data()
df.columns = ['p', 'p.cluster_adjust', 'p.group_adjust', 'gene', 'tier', 'max_deltapsi', 'start', 'end']
df.to_csv('~/DATA/splicing/Analysis/HR_leafviz/BRCA_free/BRCA_free_cluster_significance.permute_adjusted_388_extended.txt', sep='\t')
df[(df['tier'] != 0) & (df['p.group_adjust'] < 0.05)].sort_values('max_deltapsi', ascending=False).to_csv('~/DATA/splicing/Analysis/HR_leafviz/BRCA_free/BRCA_free_cluster_significance.permute_adjusted_388_extended.FDR_0.05_tier123.txt', sep='\t')


for name, group in grouped:
    print name
    print group['deltapsi'].max()
fig, ax = plt.subplot()
sns.distplot(df[df['tier'] != 0]['p.group_adjust'], kde=False, norm_hist=True, bins=10)
sns.distplot(df[df['tier'] == 0]['p.group_adjust'], kde=False, norm_hist=True, bins=10)

grouped = df.groupby('genes')
new_df = pd.DataFrame(index=df['genes'].unique(), columns=['FDR<0.05', 'FDR<0.01', 'tier'])

for name, group in grouped:
    new_df.loc[name, 'FDR<0.05'] = int(np.any(group['p.group_adjust'] < 0.05))
    new_df.loc[name, 'FDR<0.01'] = int(np.any(group['p.group_adjust'] < 0.01))
    try:
        new_df.loc[name, 'tier'] = tier_genes.loc[name, 'Tier']
    except KeyError:
        new_df.loc[name, 'tier'] = 0

a11 = new_df[new_df['tier'] != 0]['FDR<0.01'].sum()
a12 = new_df[new_df['tier'] != 0]['FDR<0.01'].shape[0] - a11
a21 = new_df[new_df['tier'] == 0]['FDR<0.01'].sum()
a22 = new_df[new_df['tier'] == 0]['FDR<0.01'].shape[0] - a21
print a11, a12, a21, a22
stats.chi2_contingency([[a11, a12], [a21, a22]])

for i in df[df['p.group_adjust'] < 0.01]['genes'].unique():
    print i

df.head(10)
