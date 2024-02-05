import os
import pandas as pd

norm_file = pd.read_csv(snakemake.input['norm_file'])
bam_file = snakemake.input['bam_file']
out_dir = snakemake.output[0]
bl = snakemake.params['blacklist']
log_file = snakemake.log[0]
threads = snakemake.threads
sample = bam_file.split('/')[1]
#broad_marks = snakemake.params['broad']
#narrow_marks = snakemake.params['narrow']
#ctrl_bam = snakemake.params['ctrl_bam']

## Check delimiters here
exp = bam_file.split('/')[1].split('-')[1]
# if exp.upper() in broad_marks:
#     broad =1
# else:
#     broad = 0

#ctrl = [x for x in os.listdir('pre-analysis') if 'Cas9_'+exp in x][0]
ctrl = [x for x in norm_file['sample'] if 'Cas9-'+exp in x][0]
ctrl_file = 'pre-analysis/'+ctrl+'/bowtie2_dros/aligned.primary.rmdup.bam'

#exp_out_file = out_dir + sample + '.bw'
print(out_dir)
ctrl_out_file = '/'.join([x for x in out_dir.split('/') if '.' not in x]) + '/' + ctrl + '.bw'
print(ctrl_out_file)
#os.system('mkdir ' + out_dir)
if 'gst' in sample.lower():
    os.system('touch ' + out_file)
else:
    ## This is for the single control per experiment
    dros_exp = norm_file[norm_file['sample'].str.contains(sample)]['dros_perc'].item()
    dros_ctrl = norm_file[norm_file['sample'].str.contains(ctrl)]['dros_perc'].item()
    if dros_exp>dros_ctrl:
        scaling_factor = 1
        ctrl_scale = dros_exp/dros_ctrl
    else:
        scaling_factor = dros_ctrl/dros_exp
        ctrl_scale = 1
    ## This is for the mapr analysis and one control for multiple samples
    # scaling_factor = norm_file[norm_file['sample'].str.contains(sample)]['scaling_new'].item()
    # if 'ct' in sample.lower():
    #     ctr_scale = norm_file[norm_file['sample']=='CT_RNAseH1_GST']['scaling_new'].item()
    #     ctrl_file = 'pre-analysis/CT_RNAseH1_GST/bowtie2/aligned.primary.rmdup.bam'
    # elif 'dmso' in sample.lower():
    #     ctr_scale = norm_file[norm_file['sample']=='DMSO_RNAseH1_GST']['scaling_new'].item()
    #     ctrl_file = 'pre-analysis/DMSO_RNAseH1_GST/bowtie2/aligned.primary.rmdup.bam'
    # else:
    #     raise Exception('Check you are running the right program!')
    print(scaling_factor)
    print(ctrl_scale)

    ## This is for spiker based normalization, broken now

    # if broad:
    #     print('Running broad for ' + exp)
    #     print('spiker.py --spikeIn --broad  --csf {} --tsf {} -t {} -c {} -o {} -p {}'.format(ctrl_scale,scaling_factor,bam_file,ctrl_file,out_dir + '/NA',threads))
    #     os.system('spiker.py --spikeIn --bw  --broad --q-peak 0.15 --q-link 0.15 --csf {} --tsf {} -t {} -c {} -o {} -p {}'.format(ctrl_scale,scaling_factor,bam_file,ctrl_file,out_dir+'/NA',threads))
    # else:
    #     print('Running narrow for ' + exp)
    #     print('spiker.py --spikeIn --csf {} --tsf {} -t {} -c {} -o {} -p {}'.format(ctrl_scale,scaling_factor,bam_file,ctrl_file,out_dir + '/NA',threads))
    #     os.system('spiker.py --spikeIn --bw  --q-peak 0.15 --q-link 0.15 --csf {} --tsf {} -t {} -c {} -o {} -p {}'.format(ctrl_scale,scaling_factor,bam_file,ctrl_file,out_dir+'/NA',threads))
    

    ## This is for deep tools/ samtools based normalization
    #os.system('samtools view -b -o ' + out_file + ' --subsample ' + str(scaling_factor) + ' ' + bam_file)
    print('bamCoverage --binSize 1 --scaleFactor '+ str(scaling_factor) +' -b ' + bam_file  + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_dir + ' -p ' + str(threads) + ' > ' + log_file)
    os.system('bamCoverage --binSize 1 --scaleFactor '+ str(scaling_factor) +' -b ' + bam_file + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_dir + ' -p ' + str(threads) + ' > ' + log_file)
    print('bamCoverage --binSize 1 --scaleFactor '+ str(ctrl_scale) +' -b ' + ctrl_file + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + ctrl_out_file + ' -p ' + str(threads) + ' > ' + log_file)
    os.system('bamCoverage --binSize 1 --scaleFactor '+ str(ctrl_scale) +' -b ' + ctrl_file + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + ctrl_out_file + ' -p ' + str(threads) + ' > ' + log_file)
    #print('bamCompare --binSize 1 -b1 ' + bam_file + ' -b2 '+ ctrl_bam + ' --scaleFactorsMethod None --scaleFactors '+ str(scaling_factor) + ':' + str(ctr_scale) + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_file + ' -p ' + str(threads) + ' > ' + log_file)
    #os.system('bamCompare --binSize 1 -b1 ' + bam_file + ' -b2 '+ ctrl_bam + ' --scaleFactorsMethod None --scaleFactors '+ str(scaling_factor) + ':' + str(ctr_scale) + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_file + ' -p ' + str(threads) + ' > ' + log_file)