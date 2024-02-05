import os
import pandas as pd

stats_df = pd.DataFrame()
inputs = list(set(snakemake.input))
#cells = snakemake.params['cells']
#exps = snakemake.params['exps']

for sample in inputs:
    if 'gst' in sample.lower():
        continue
    print(sample)
    df = pd.read_csv(sample,sep='\t')
    df['sample'] = sample.split('/')[1]
    df['sample_path'] = sample
    if stats_df.empty: 
        stats_df = df.copy(deep=True)
    else:
        stats_df = pd.concat([stats_df,df])
    #print(stats_df)

samples = [x.split('/')[1] for x in inputs]
#histones = [x for x in samples if 'H3K' in x or 'H4K' in x]
#tfs = [x for x in samples if x not in histones]

## For every control and experiment, calculate the scaling factor.
## If it's cas9 then don't look at cell lines, otherwise cell lines



max_dros = stats_df['n_exogenous'].max()
stats_df['dros_perc'] = 100*stats_df['n_exogenous']/(stats_df['n_exogenous']+stats_df['n_sample'])
stats_df['scaling_factor'] = stats_df['n_exogenous'].apply(lambda x: max_dros/x)
scale_max = (stats_df['n_exogenous']/stats_df['n_sample']).max()
#seq_min = align_stats['Sequencing Depth'].min()
#stats_df['seq_norm'] = align_stats
stats_df['scaling_sam'] = stats_df['n_exogenous']/(stats_df['n_sample']*scale_max)

stats_df.to_csv(snakemake.output[0])
