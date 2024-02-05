rule findmotifsgenome:
    input:
        'pre-analysis/{sample}/macs/narrow/NA_peaks.narrowPeakFiltered'
    output:
        directory('analysis/homer/{sample}')
    threads:
        4
    shell:
        'findMotifsGenome.pl {input} hg38 {output} -size 200 -p {threads} -S 0 -nomotif'