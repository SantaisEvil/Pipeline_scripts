import os
#import pandas as pd
import matplotlib.pyplot as plt
from matplotlib_venn import venn3

#ar_files = snakemake.input
#out_file = snakemake.output[0]
#tmp_dir = snakemake.params[0]
in_dir = '/media/asangani1/NSD2_Project/2023-11-28/pre-analysis2/'
file = 'macs/narrow/NA_peaks.narrowPeak'
ar_files = ['LNCaP_Input_TM','LNCaP_AR_WT','LNCaP_AR_Y1092A']
print(os.path.join(in_dir,ar_files[0],file))
ar_files = [os.path.join(in_dir,x,file) for x in ar_files]
print(ar_files)

tmp_dir = '/home/sven/Data/temp_files/'
out_file = os.path.join(tmp_dir,'out.png')

#samples = [x.split('_')[2].split('/')[0] for x in ar_files]
samples = ['LNCaP_Input_TM','LNCaP_AR_WT','LNCaP_AR_Y1092A']
print(samples)
a_samp = samples[0]
b_samp = samples[1]
c_samp = samples[2]
a = ar_files[0]
b = ar_files[1]
c = ar_files[2]

os.system('bedtools intersect -wa -a '+ a+' -b '+b + ' > ' + tmp_dir + a_samp+'_'+b_samp)
os.system('bedtools intersect -wa -a '+ a+' -b '+c + ' > ' + tmp_dir + a_samp+'_'+c_samp)
os.system('bedtools intersect -wa -a '+ b+' -b '+c + ' > ' + tmp_dir + b_samp+'_'+c_samp)
os.system('bedtools intersect -a '+ b+' -b '+c + ' ' + a + ' > ' + tmp_dir + 'all')

os.system('bedtools intersect -v -a '+ a+' -b '+ tmp_dir + 'all' + ' > ' + tmp_dir + a_samp)
os.system('bedtools intersect -v -a '+ b+' -b '+ tmp_dir + 'all' ' > ' + tmp_dir + b_samp)
os.system('bedtools intersect -v -a '+ c+' -b '+ tmp_dir + 'all' + ' > ' + tmp_dir + c_samp)

only_a = len(open(tmp_dir + a_samp,'r').readlines())
only_b = len(open(tmp_dir + b_samp,'r').readlines())
only_c = len(open(tmp_dir + c_samp,'r').readlines())

a_b = len(open(tmp_dir + a_samp+'_'+b_samp,'r').readlines())
a_c = len(open(tmp_dir + a_samp+'_'+c_samp,'r').readlines())
b_c = len(open(tmp_dir + b_samp+'_'+c_samp,'r').readlines())

all = len(open(tmp_dir + 'all','r').readlines())

venn3(subsets=(only_a,only_b,a_b,only_c,a_c,b_c,all),set_labels=samples)
plt.savefig(out_file)