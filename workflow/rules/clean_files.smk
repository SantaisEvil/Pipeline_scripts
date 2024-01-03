samples,reads = glob_wildcards("data/{sample}_{read}_001.fastq.gz")
samples = list(set(samples))

rule clean_alignment:
    input:
        expand('pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam',sample=samples)
    output:
        'pre-analysis/{sample}/bowtie2/cleaned.txt'
    threads:
        16
    shell:
        '''
        rm pre-analysis/*/bowtie2/*.sam
        rm pre-analysis/*/bowtie2/*markDup*
        rm pre-analysis/*/bowtie2/*Aligned.sort*
        touch pre-analysis/*/bowtie2/cleaned.txt
        '''