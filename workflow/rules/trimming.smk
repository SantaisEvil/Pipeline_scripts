rule atria_trim_PE:
    input:
        r1 = 'data/{sample}_R1.fastq.gz',
        r2 = 'data/{sample}_R2.fastq.gz'
    output:
        r1 = 'data_trimmed/{sample}_R1.atria.fastq.gz',
        r2 = 'data_trimmed/{sample}_R2.atria.fastq.gz'
    threads:
        4
    log:
        'data_trimmed/logs/{sample}_atria.log'
    shell:
        'atria --read1 {input.r1} --read2 {input.r2} -o data_trimmed 2> {log}'


# rule markDup:
#     input:
#         bamFile= "pre-analysis/{sample}/bowtie2/Aligned.sortedByCoord.out.bam"     
#     output:
#         bamMarked= "pre-analysis/{sample}/bowtie2/aligned_markDup.bam",
#         dupStats="pre-analysis/{sample}/bowtie2/dupMetrics.tsv" 
#     params:
#         samtools="-ASO=coordinate --TAGGING_POLICY All"
#     threads: 4

#     log:
#         "pre-analysis/{sample}/logs/MarkDuplicates.log"

#     shell:
#         "gatk MarkDuplicates -I {input.bamFile} -O  {output.bamMarked} -M {output.dupStats} 2> {log}"
