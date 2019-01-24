import os
import pandas as pd
import gzip

DIR_PATH = "/home/haeun/2TB_disk/splicing/data"
FILE_LIST = [i for i in os.listdir(DIR_PATH + '/COSMIC/genotypes') if i.endswith('.gz')]

def get_info(infofile):
    info = pd.read_csv(infofile, sep='\t', header=0)
    info = info.dropna(subset=['COSMIC_ID'])
    info = info.set_index('COSMIC_ID')
    info = info.astype(str)
    info.index = info.index.astype(int).astype(str)
    return info

def changeFileName(infofile):
    site_dic = dict()
    info = get_info(infofile)

    for i in FILE_LIST:
        try:
            tokens = info.loc[i.split('_')[0], ['Line', 'Site', 'Gender']].values
            new_name = '{0}_{1}'.format(tokens[0], tokens[1].upper())
            site_name = tokens[1].upper()
            gender = tokens[2]
            try:
                _ = site_dic[site_name]
            except KeyError:
                site_dic[site_name] = 1
                os.system('mkdir ' + DIR_PATH + '/COSMIC/genotypes/' + site_name)


            os.system('mv {0}/COSMIC/genotypes/{1} {0}/COSMIC/genotypes/{2}/{3}'.format(DIR_PATH, i, site_name, new_name))
            # make gender list
            with open('{0}/COSMIC/genotypes/{1}/list.gender.txt'.format(DIR_PATH, site_name),'a') as gen_file:
                if gender=='M':
                    gen_file.write(new_name+'\t1\n')
                elif gender=='F':
                    gen_file.write(new_name+'\t2\n')
                else:
                    gen_file.write(new_name+'\tother\n')

            # make filename list
            with open('{0}/COSMIC/genotypes/{1}/filename.txt'.format(DIR_PATH, site_name),'a') as name_file:
                name_file.write(new_name+'\n')

        except KeyError:
            print i, "not in cell line info."

changeFileName(DIR_PATH + '/Cell_lines_CCLE_COSMIC_merged_info.txt')
