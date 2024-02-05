rule bed2fasta:
    input:
        'pre-analysis/{sample}/macs/narrow/NA_peaks.narrowPeakFiltered'
    output:
        'analysis/meme/{sample}.fa'
    threads:
        4
    params:
        genome_fa = config['ref']['fasta']
    shell:
        'bed2fasta -o {output} {input} {params.genome_fa}'

rule streme:
    input:
        'analysis/meme/{sample}.fa'
    output:
        'analysis/meme/{sample}/streme.html'
    threads:
        4
    shell:
        'streme -o {output} --p {input}'