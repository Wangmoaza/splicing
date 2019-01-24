import os
import pandas as pd
import numpy as np
from statsmodels.stats.multitest import fdrcorrection

PREFIX = "All"
PATH = "/home/haeun/2TB_disk/splicing/Analysis/HR_leafviz/{0}/results/".format(PREFIX)


def return_qval(p_vals):
    rejected, q = fdrcorrection(p_vals)
    return q


def get_permute_q():
    file_list = sorted([i for i in os.listdir(PATH) if i.endswith('_significance.txt')])

    # merge permuted p-values
    df = pd.DataFrame()
    for fname in file_list:
        ser = pd.read_csv(PATH + fname,
                          sep='\t', header=0, index_col=0, na_values='NA',
                          usecols=['cluster', 'p'], squeeze=True)
        if "permute" in fname:
            ser.name = fname[11:14]
        else:
            ser.name = 'original'
        df = pd.concat([df, ser], axis=1)

    # apply fdr correction
    df = df.dropna(axis=0)
    qvals = np.zeros(df.shape)
    for i in range(df.shape[0]):
        qvals[i] = return_qval(df.iloc[i, :].values)
    adj_df = pd.DataFrame(qvals, index=df.index, columns=df.columns)
    return adj_df['original']


def main():
    q = get_permute_q()
    ori_df = pd.read_csv(PATH + 'All_cluster_significance.txt', sep='\t', header=0, index_col=0)
    new_df = pd.concat([ori_df[['p', 'p.adjust', 'genes']], q], axis=1, join='inner')
    new_df = new_df[['p', 'p.adjust', 'original', 'genes']]
    new_df.columns = ['p', 'p.cluster_adjust', 'p.group_adjust', 'genes']
    new_df.index.name = "cluster"
    new_df.sort_values('p.group_adjust').to_csv(PATH + 'All_cluster_significance.permute_adjusted.txt', sep='\t')


if __name__ == "__main__":
    main()
