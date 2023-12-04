ruleorder: bowtie2_PE > bowtie2_SE

rule bowtie2_PE:
    priority: 2

    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"
        
        
    output:
        aligned_sam="pre-analysis/{sample}/bowtie2/unfiltered_aligned.sam"
    log:
        alignment_stats = "pre-analysis/{sample}/bowtie2/unfiltered_aligned.txt",
        met = "pre-analysis/{sample}/logs/bowtie2.txt"
    threads: 12
    params:
        bowtie2 = "bowtie2 --local --no-mixed --no-discordant ",
        index = config["ref"]["bowtie2"]
    conda:
        "ngsmo"

    shell:
        """
        {params.bowtie2} -p {threads} --met-file {log.met} -x {params.index} -1 {input.r1} -2 {input.r2} -S {output.aligned_sam}  &> {log.alignment_stats}
        
        """
rule bowtie2_SE:
    priority: 2

    input:
        "data/{sample}_R1.fastq.gz"
        
        
    output:
        aligned_sam="pre-analysis/{sample}/bowtie2/unfiltered_aligned.sam"
    log:
        "pre-analysis/{sample}/bowtie2/unfiltered_aligned.txt"
    threads: 6
    params:
        bowtie2 = "bowtie2 --no-mixed --no-discordant -p",
        index = config["ref"]["bowtie2"]
    conda:
        "ngsmo"

    shell:
        """
        {params.bowtie2} {threads} -x {params.index} -U {input} -S {output.aligned_sam}  &> {log}
        
        """

rule bwa:
    priority: 2
    
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"
        
    output:
        aligned_sam="pre-analysis/{sample}/bwa/unfiltered_aligned.sam"
    threads: 5
    params:
        bwa = "-M",
        index = config["ref"]["bwa"]
    conda:
        "ngsmo"

    shell:
        """
        bwa mem {params.index} {input.r1} {input.r2} -t {threads} -M -S  > {output.aligned_sam} 
        """
rule STAR:
    priority: 2
    
    input:
        r1 = "data/{sample}_R1.fastq.gz",
        r2 = "data/{sample}_R2.fastq.gz"

    output:
        aligned_bam="pre-analysis/{sample}/STAR/Aligned.sortedByCoord.out.bam"

    log:
        detailed_log="pre-analysis/{sample}/STAR/Log.out",
        mapping_summary="pre-analysis/{sample}/STAR/Log.final.out"

    
    threads: 8

    params:
        star = "--outSAMtype BAM SortedByCoordinate \
        --alignEndsType EndToEnd --alignIntronMax 1 \
        --outSAMattributes NH NM MD --outSAMunmapped Within --outSAMmapqUnique 50 --limitBAMsortRAM 24000000000 \
        --readFilesCommand zcat",
        
        index = config["ref"]["star"],
        outdir = "pre-analysis/{sample}/STAR/"
    conda:
        "ngsmo"

    shell:
        """
        STAR --genomeDir {params.index} --readFilesIn {input.r1} {input.r2} --runThreadN {threads} {params.star} --outFileNamePrefix {params.outdir}
    
        """