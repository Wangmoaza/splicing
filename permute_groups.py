import pandas as pd
import numpy as np

ser = pd.read_csv('/home/haeun/DATA/splicing/Analysis/HR_leafviz/BRCA_free/groups/BRCA_free.HRD_sig3_median_groups.txt',
                  sep='\t', header=None, index_col=0, squeeze=True)

for i in range(1000):
    np.random.shuffle(ser.values)
    ser.to_csv('/home/haeun/DATA/splicing/Analysis/HR_leafviz/BRCA_free/groups/BRCA_free.HRD_sig3_median_groups.permute{0:03}.txt'.format(i),
               sep='\t', header=False)
