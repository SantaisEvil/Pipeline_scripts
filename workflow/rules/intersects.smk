samples,reads = glob_wildcards("data/{sample}_{read}.fastq.gz")
samples = list(set(samples))
ar = [x for x in samples if 'ar' in x.lower()]
print(ar)

rule venn3:
    input:
        expand('pre-analysis/{sample}/macs/narrow/NA_peaks.narrowPeakFiltered',sample=ar)
    output:
        'analysis/ar_venn.png'
    params:
        'analysis/intersect_temps/'
    threads:
        1
    conda:
        'ngs'
    script:
        '../scripts/intersect.py'