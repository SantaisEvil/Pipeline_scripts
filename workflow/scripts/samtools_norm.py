import os
import pandas as pd

norm_file = pd.read_csv(snakemake.input['norm_file'])
bam_file = snakemake.input['bam_file']
out_file = snakemake.output[0]
bl = snakemake.params['blacklist']
log_file = snakemake.log[0]
threads = snakemake.threads
sample = bam_file.split('/')[1]
index_file = bam_file + '.bai'
ctrl_bam = snakemake.params['ctrl_bam']

if 'gst' in sample.lower():
    os.system('touch ' + out_file)
else:
    scaling_factor = norm_file[norm_file['sample'].str.contains(sample)]['scaling_sam_view'].item()
    if 'ct' in sample.lower():
        ctr_scale = norm_file[norm_file['sample']=='CT_RNAseH1_GST']['scaling_bam'].item()
    elif 'dmso' in sample.lower():
        ctr_scale = norm_file[norm_file['sample']=='DMSO_RNAseH1_GST']['scaling_bam'].item()
    else:
        raise Exception('Check you are running the right program!')
    print(scaling_factor)
    print(ctr_scale)
    #os.system('samtools view -b -o ' + out_file + ' --subsample ' + str(scaling_factor) + ' ' + bam_file)
    #print('bamCoverage --binSize 10 --scaleFactor '+ str(scaling_factor) +' -b ' + bam_file + ' ' + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_file + ' -p ' + str(threads) + ' > ' + log_file)
    #os.system('bamCoverage --binSize 10 --scaleFactor '+ str(scaling_factor) +' -b ' + bam_file + ' ' + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_file + ' -p ' + str(threads) + ' > ' + log_file)
    #print('bamCompare --binSize 1 -b1 ' + bam_file + ' -b2 '+ ctrl_bam + ' --operation subtract--scaleFactorsMethod None --scaleFactors '+ str(scaling_factor) + ':' + str(ctr_scale) + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_file + ' -p ' + str(threads) + ' > ' + log_file)
    #os.system('bamCompare --binSize 1 -b1 ' + bam_file + ' -b2 '+ ctrl_bam + ' --operation subtract --scaleFactorsMethod None --scaleFactors '+ str(scaling_factor) + ':' + str(ctr_scale) + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_file + ' -p ' + str(threads) + ' > ' + log_file)
    #print('samtools view -b -u -o ' + out_file + ' --subsample ' + str(scaling_factor) + ' ' + bam_file)
    #os.system('samtools view -b -u -o ' + out_file + ' --subsample ' + str(scaling_factor) + ' ' + bam_file)
    #print('bamCompare --binSize 1 -b1 ' + bam_file + ' -b2 '+ ctrl_bam + ' --scaleFactorsMethod None --scaleFactors --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_file + ' -p ' + str(threads) + ' > ' + log_file)
    #os.system('samtools index -b '+ bam_file + ' -o ' + index_file)
    #os.system('bamCompare --binSize 1 -b1 ' + bam_file + ' -b2 '+ ctrl_bam + ' --scaleFactorsMethod None --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_file + ' -p ' + str(threads) + ' > ' + log_file)
    print('bamCoverage --binSize 10 -b ' + bam_file + ' ' + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_file + ' -p ' + str(threads) + ' > ' + log_file)
    os.system('bamCoverage --binSize 10 -b ' + bam_file + ' ' + ' --centerReads --normalizeUsing CPM -bl ' + bl + ' -o ' + out_file + ' -p ' + str(threads) + ' > ' + log_file)