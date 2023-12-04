ruleorder: fastqc  > assembleStats

rule fastqc:
    priority: 1
    
    input:
        rawread="data/{sample}_{read}.fastq.gz"
    output:
        zip ="pre-analysis/{sample}/fastqc/{sample}_{read}_fastqc.zip",
        html="pre-analysis/{sample}/fastqc/{sample}_{read}_fastqc.html"
    threads: 1
    
    conda:
        "ngsmo"

    params:
        path="pre-analysis/{sample}/fastqc/"
    shell:
        "fastqc {input.rawread} --threads {threads} -o {params.path}"

rule assembleStats:
    priority: 10
    
    output:
        "pre-analysis/alignment_stats.csv"
    params:
        files = "pre-analysis/*",
        dup_files = "pre-analysis/*/bowtie2/dupMetrics.tsv",
        rmdup_files = "pre-analysis/*/bowtie2/dupMetricsFiltered.tsv"

    script:
        "../scripts/bowtie_alignment_stats.R"
    