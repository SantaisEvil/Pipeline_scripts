import os
import pandas as pd

norm_file = snakemake.input['norm']
bam_file = snakemake.input['bam_file']
inp_files = snakemake.input['inp_files']
out_dir = snakemake.output[0]
sample = bam_file.split('/')[1]

ctrl_file = [x for x in inp_files if 'gst' in x.lower() and x.split('_')[0]==bam_file.split('_')[0]][0]

norm_df = pd.read_csv(norm_file)
print()
bam_sf = norm_df[(norm_df['sample'].str.contains(bam_file.split('/')[1])) & ~(norm_df['sample'].str.contains('gst'))]['scaling_bam'].item()
ctrl_sf = norm_df[norm_df['sample'].str.contains(ctrl_file.split('/')[1])]['scaling_bam'].item()
print(bam_sf)
print(ctrl_sf)
#print('bam - ' + str(bam_sf) + 'ctrl - ' + str(ctrl_sf))

print('spiker.py --spikeIn -p 4 --csf ' + str(ctrl_sf) + ' --tsf ' + str(bam_sf) + ' -t ' + bam_file + ' -c ' + ctrl_file + ' -o pre-analysis/' +sample + '/normalized_files/' + sample)
os.system('spiker.py --spikeIn -p 4 --csf ' + str(ctrl_sf) + ' --tsf ' + str(bam_sf) + ' -t ' + bam_file + ' -c ' + ctrl_file + ' -o pre-analysis/' +sample + '/normalized_files/' + sample)