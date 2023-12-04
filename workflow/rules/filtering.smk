rule markDup:
    input:
        bamFile= "pre-analysis/{sample}/bowtie2/Aligned.sortedByCoord.out.bam"     
    output:
        bamMarked= "pre-analysis/{sample}/bowtie2/aligned_markDup.bam",
        dupStats="pre-analysis/{sample}/bowtie2/dupMetrics.tsv" 
    params:
        samtools="-ASO=coordinate --TAGGING_POLICY All"
    threads: 4

    conda:
        "ngsmo"
    log:
        "pre-analysis/{sample}/logs/MarkDuplicates.log"

    shell:
        "gatk MarkDuplicates -I {input.bamFile} -O  {output.bamMarked} -M {output.dupStats} 2> {log}"


rule markDupFiltered:
    input:
        bamFile= "pre-analysis/{sample}/bowtie2/aligned_markDup.bam"     
    output:
        bamMarked= "pre-analysis/{sample}/bowtie2/aligned_markDupFiltered.bam",
        dupStats="pre-analysis/{sample}/bowtie2/dupMetricsFiltered.tsv" 
    params:
        samtools="-ASO=coordinate --TAGGING_POLICY All"
    threads: 4

    conda:
        "ngsmo"
    log:
        "pre-analysis/{sample}/logs/MarkDuplicates.log"

    shell:
        "gatk MarkDuplicates -I {input.bamFile} -O  {output.bamMarked} -M {output.dupStats} 2> {log}"

rule samtools:
    input:
        "pre-analysis/{sample}/bowtie2/unfiltered_aligned.sam"
        
    output:
        "pre-analysis/{sample}/bowtie2/unfiltered_aligned_stats.tsv"
    threads: 4
    params:
        bwa = "-M",
        index = config["ref"]["bowtie2"]
    conda:
        "ngsmo"

    shell:
        """
        samtools flagstat -O tsv -@ {threads} {input} > {output}
        
        """

rule samToBamSorted:
    input:
        samFile="pre-analysis/{sample}/bowtie2/unfiltered_aligned.sam"         
    output:
        bamFile= "pre-analysis/{sample}/bowtie2/Aligned.sortedByCoord.out.bam"
    params:
        samtools="sort --write-index -O BAM"
    threads: 4

    conda:
        "ngsmo"

    shell:
        "samtools {params.samtools} -@ {threads}  -o {output.bamFile}  {input.samFile} "

rule filterBam:
    input:
        "pre-analysis/{sample}/bowtie2/aligned_markDup.bam"      
    output:
        filteredBam="pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam",
        index="pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam.bai"
    params:
        filter="view -h -f bam --num-filter /3340",
        index="index --show-progress "
    threads: 4

    conda:
        "sambamba"

    shell:
        """
        sambamba {params.filter} -t {threads} {input} > {output.filteredBam}
        sambamba {params.index} -t {threads} {output.filteredBam} {output.index}
        
        """
rule removeOverlap:
    input:
        bamFile = "pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam",
        bedsToRemove = config["ref"]["bedsToRemove"]

    output:
        filteredBam="pre-analysis/{sample}/bowtie2/aligned.primary.rmdupAndOverlap.bam"
        # index="pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam.bai"
    params:
        "-v"
    threads: 4

    conda:
        "ngsmo"

    shell:
        """
        bedtools intersect -a {input.bamFile} -b {input.bedsToRemove} {params} > {output} 
        
        """


