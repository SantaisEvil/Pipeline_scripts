import os

samples = snakemake.input
out_dir = snakemake.output[0]
root_dir = snakemake.params['dir']
norm_dir = snakemake.params['norm_dir']

try:
    os.system('mkdir pre-analysis/normalized_files')
except:
    pass

for sample in samples:
    sample_name = sample.split('.')[0].split('/')[-1]
    print('ln -s ' +root_dir+ sample + ' '+ root_dir+norm_dir+sample_name + '.bw')
    #os.system('ln -s ' +root_dir+ sample + ' '+ root_dir+norm_dir+sample_name + '.bw')
    os.system('cp '+root_dir+ sample + ' '+ root_dir+norm_dir+sample_name + '.bw')