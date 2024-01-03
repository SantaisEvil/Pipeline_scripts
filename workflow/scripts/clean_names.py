import os
import re

def clean_names(in_dir):
    files = os.listdir(in_dir)
    for file in files:
        split_file = file.split('_')
        #print(re.search('S\d+',split_file[1]))
        split_file = [x for x in split_file if not (re.search('S\d+',x))]
        split_file = ('_').join(split_file)
        split_file = split_file.replace('_001','')
        #os.system('mv ' + os.path.join(in_dir,file) + ' ' + os.path.join(in_dir,split_file))
        print('mv ' + os.path.join(in_dir,file) + ' ' + os.path.join(in_dir,split_file))

if __name__=='__main__':
    clean_names('/media/asangani1/NSD2_Project/2023-11-28_rerun/fastqs')