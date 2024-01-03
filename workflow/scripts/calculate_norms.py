import os
import pandas as pd

stats_df = pd.DataFrame()
inputs = list(set(snakemake.input))
for sample in inputs:
    if 'gst' in sample.lower():
        continue
    print(sample)
    df = pd.read_csv(sample,sep='\t')
    df['sample'] = sample
    if stats_df.empty: 
        stats_df = df.copy(deep=True)
    else:
        stats_df = pd.concat([stats_df,df])
    print(stats_df)

max_dros = stats_df['n_exogenous'].max()

stats_df['scaling_factor'] = stats_df['n_exogenous'].apply(lambda x: max_dros/x)
scale_max = (stats_df['n_exogenous']/stats_df['n_sample']).max()
#seq_min = align_stats['Sequencing Depth'].min()
#stats_df['seq_norm'] = align_stats
stats_df['scaling_sam'] = stats_df['n_exogenous']/(stats_df['n_sample']*scale_max)

stats_df.to_csv(snakemake.output[0])