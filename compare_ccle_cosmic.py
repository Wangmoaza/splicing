import pandas as pd
import numpy as np
import os

PATH = "/home/haeun/2TB_disk/splicing/data"


def import_cell_list(study, tissue):
    """ Import cell line names. """
    ls = os.listdir('{0}/{1}/tmp_{2}'.format(PATH, study.upper(), tissue.upper()))
    return [i.split('.')[0] for i in ls]


def compare(tissue):
    tissue = tissue.upper()
    output_file = PATH + '/CCLE_COSMIC_{0}_genotype_concordance.txt'.format(tissue)
    with open(output_file, 'w') as f:
        f.write('{0}\t{1}\t{2}\t{3}\n'.format('Cell Line', 'Total', 'Same Call', 'Percentage'))

    ccle_list = import_cell_list('ccle', tissue)
    cosmic_list = import_cell_list('cosmic', tissue)

    overlap_cells = set(ccle_list) & set(cosmic_list)
    print tissue
    print "CCLE: {0} COSMIC: {1} overlap: {2}".format(len(ccle_list), len(cosmic_list), len(overlap_cells))
    for i in overlap_cells:
        ccle_ser = pd.read_csv('{0}/CCLE/tmp_{1}/{2}.input'.format(PATH, tissue.upper(), i),
                               sep='\t', header=0, index_col=0, squeeze=True)
        cosmic_ser = pd.read_csv('{0}/COSMIC/tmp_{1}/{2}.input'.format(PATH, tissue.upper(), i),
                                 sep='\t', header=0, index_col=0, squeeze=True)

        ccle_ser.name = 'CCLE'
        cosmic_ser.name = 'COSMIC'
        # drop duplicated index
        cosmic_ser = cosmic_ser[~cosmic_ser.index.duplicated(keep='first')]
        ccle_ser = ccle_ser[~ccle_ser.index.duplicated(keep='first')]

        # merge two sets
        merged = pd.concat([ccle_ser, cosmic_ser], axis=1, sort=True).dropna()

        concord = merged.agg(lambda x: x[0] == x[1], axis=1).values
        total_len = len(concord)
        same_len = np.count_nonzero(concord)
        with open(output_file, 'a') as f:
            f.write('{0}\t{1}\t{2}\t{3}\n'.format(i, total_len, same_len, np.round(same_len / float(total_len), 5)))


def main():
    for i in ['breast', 'pancreas', 'ovary']:
        compare(i)


main()
