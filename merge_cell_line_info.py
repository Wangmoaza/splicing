import pandas as pd
import numpy as np

def remove_special(string):
    string = str(string)
    return ''.join(e for e in string if e.isalnum())

ccle = pd.read_csv('~/2TB_disk/splicing/data/CCLE/CCLE_sample_info_file_2012-10-18.txt', sep='\t', header=0, index_col=0)
ccle = ccle.astype(str)
cosmic = pd.read_excel('~/2TB_disk/splicing/data/GDSC/Cell_Lines_Details.xlsx', sheet_name='COSMIC tissue classification', header=0)
cosmic = cosmic.astype(str)
cosmic['index-name'] = cosmic['Line'].str.upper().apply(remove_special)
cosmic = cosmic.set_index('index-name')
cosmic = cosmic.loc[~cosmic.index.duplicated(), :]


ccle_name1 = pd.DataFrame({'CCLE_Name': ccle.index, 'Gender': ccle['Gender'].values}, index=ccle['Cell line primary name'].str.upper().apply(remove_special).values)
ccle_name1.index.name = 'CCLE Name'
ccle_name1 = ccle_name1[~ccle_name1.index.duplicated()]
ccle_name1.index.is_unique
merged = pd.concat([cosmic, ccle_name1], axis=1, sort=True)
merged.columns = ['COSMIC_Name', 'COSMIC_ID', 'Site', 'Histology', 'CCLE_Name', 'Gender']
merged.index.name = 'Line'
for idx in merged[merged['Site'].isna()].index:
    merged.loc[idx, 'Site'] = merged.loc[idx, 'CCLE_Name'].split('_', 1)[1].lower()

merged.sort_values('Site').to_csv('~/2TB_disk/splicing/data/Cell_lines_CCLE_COSMIC_merged_info.txt', sep='\t')
