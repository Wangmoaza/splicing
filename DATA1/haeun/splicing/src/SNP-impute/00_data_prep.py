import os

DIR_PATH = "/home/omics/DATA1/haeun/splicing/data/CCLE/variant"
FILE_LIST = os.listdir(DIR_PATH + '/birdseed')

def changeFileName(inputfile):
    site_dic = dict()
    with open(inputfile, 'r') as f:
        f.readline() # skip header
        for line in f.readlines():
            token = line.strip().split('\t')
            id_name = token[1]  # cell line id
            try:
                site_name = id_name.split('_', 1)[1]
            except IndexError:
                site_name = "UNKNOWN"
            # make new directory on first account
            try:
                _ = site_dic[site_name]
            except KeyError:
                site_dic[site_name] = 1
                os.system('mkdir ' + DIR_PATH + '/birdseed/' + site_name)

            # change file name
            if token[0] + '.birdseed.data.txt' in FILE_LIST:
                os.system('mv {0}/birdseed/{1} {0}/birdseed/{2}/{3}'.format(DIR_PATH, token[0] + '.birdseed.data.txt', site_name, id_name))
                # make gender list
                with open('{0}/birdseed/{1}/list.gender.txt'.format(DIR_PATH, site_name),'a') as gen_file:
                    if token[2]=='M':
                        gen_file.write(id_name+'\t1\n')
                    elif token[2]=='F':
                        gen_file.write(id_name+'\t2\n')
                    else:
                        gen_file.write(id_name+'\tother\n')

                # make filename list
                with open('{0}/birdseed/{1}/filename.txt'.format(DIR_PATH, site_name),'a') as name_file:
                    name_file.write(id_name+'\n')


changeFileName(DIR_PATH + '/CCLE_SNP.Arrays.sif_2013-12-03.txt')
