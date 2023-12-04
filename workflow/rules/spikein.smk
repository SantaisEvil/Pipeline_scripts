samples,reads = glob_wildcards("data/{sample}_{read}.fastq.gz")
spiker_files = ['split_human.bam','split_exogenous.bam','split_both.bam','split_neither.bam','split.report.txt']
samples = list(set(samples))
non_gst = [x for x in samples if 'gst' not in x.lower()]


rule bowtie2_PE_dros:
    input:
        r1 = "data_trimmed/{sample}_R1.atria.fastq.gz",
        r2 = "data_trimmed/{sample}_R2.atria.fastq.gz"
        
        
    output:
        aligned_sam="pre-analysis/{sample}/bowtie2_dros/unfiltered_aligned.sam"
    log:
        alignment_stats = "pre-analysis/{sample}/bowtie2_dros/unfiltered_aligned.txt",
        met = "pre-analysis/{sample}/logs/bowtie2.txt"
    threads: 8
    params:
        bowtie2 = "bowtie2 --local --no-mixed --no-discordant ",
        index = config["ref"]["bowtie2"]
    conda:
        "ngs"

    shell:
        """
        {params.bowtie2} -p {threads} --met-file {log.met} -x {params.index} -1 {input.r1} -2 {input.r2} -S {output.aligned_sam}  &> {log.alignment_stats}
        
        """


rule markDup_dros:
    input:
        bamFile= "pre-analysis/{sample}/bowtie2_dros/Aligned.sortedByCoord.out.bam"     
    output:
        bamMarked= "pre-analysis/{sample}/bowtie2_dros/aligned_markDup.bam",
        dupStats="pre-analysis/{sample}/bowtie2_dros/dupMetrics.tsv" 
    params:
        samtools="-ASO=coordinate --TAGGING_POLICY All"
    threads: 4

    log:
        "pre-analysis/{sample}/logs/MarkDuplicates.log"

    shell:
        "gatk MarkDuplicates -I {input.bamFile} -O  {output.bamMarked} -M {output.dupStats} 2> {log}"


rule markDupFiltered_dros:
    input:
        bamFile= "pre-analysis/{sample}/bowtie2_dros/aligned_markDup.bam"     
    output:
        bamMarked= "pre-analysis/{sample}/bowtie2_dros/aligned_markDupFiltered.bam",
        dupStats="pre-analysis/{sample}/bowtie2_dros/dupMetricsFiltered.tsv" 
    params:
        samtools="-ASO=coordinate --TAGGING_POLICY All"
    threads: 4

    log:
        "pre-analysis/{sample}/logs/MarkDuplicates.log"

    shell:
        "gatk MarkDuplicates -I {input.bamFile} -O  {output.bamMarked} -M {output.dupStats} 2> {log}"

rule samtools_dros:
    input:
        "pre-analysis/{sample}/bowtie2_dros/unfiltered_aligned.sam"
        
    output:
        "pre-analysis/{sample}/bowtie2_dros/unfiltered_aligned_stats.tsv"
    threads: 4
    params:
        bwa = "-M",
        index = config["ref"]["bowtie2"]
    conda:
        "ngs"

    shell:
        """
        samtools flagstat -O tsv -@ {threads} {input} > {output}
        
        """

rule samToBamSorted_dros:
    input:
        samFile="pre-analysis/{sample}/bowtie2_dros/unfiltered_aligned.sam"         
    output:
        bamFile= "pre-analysis/{sample}/bowtie2_dros/Aligned.sortedByCoord.out.bam"
    params:
        samtools="sort --write-index -O BAM"
    threads: 4

    conda:
        "ngs"

    shell:
        "samtools {params.samtools} -@ {threads}  -o {output.bamFile}  {input.samFile} "

rule filterBam_dros:
    input:
        "pre-analysis/{sample}/bowtie2_dros/aligned_markDup.bam"      
    output:
        filteredBam="pre-analysis/{sample}/bowtie2_dros/aligned.primary.rmdup.bam",
        index="pre-analysis/{sample}/bowtie2_dros/aligned.primary.rmdup.bam.bai"
    params:
        filter="view -h -f bam --num-filter /3340",
        index="index --show-progress "
    threads: 4

    conda:
        "ngs"

    shell:
        """
        sambamba {params.filter} -t {threads} {input} > {output.filteredBam}
        sambamba {params.index} -t {threads} {output.filteredBam} {output.index}
        
        """

rule bamToBw_dros:
    input:
        "pre-analysis/{sample}/bowtie2_dros/aligned.primary.rmdup.bam"      
    output:
        "pre-analysis/ucsc/" + config["ref"]["build"] + "/" + "{sample}.bw"
    params:
        main="--binSize 1 --normalizeUsing CPM --centerReads",
        blacklist= config["ref"]["blacklist"]
    threads: 4

    log:
        "pre-analysis/{sample}/logs/BigWig.log"

    conda:
        "ngs"

    shell:
        "bamCoverage {params.main} -bl {params.blacklist} -b {input} -o {output}  -p {threads}  > {log} "

rule split_bam:
    input:
        'pre-analysis/{sample}/bowtie2_dros/aligned.primary.rmdup.bam'
    output:
        'pre-analysis/{sample}/spiker/split.report.txt'
    threads:
        4
    params:
        prefix = 'pre-analysis/{sample}/spiker/split'
    conda:
        'spiker'
    shell:
        'split_bam.py --threads {threads} -i {input} -o {params.prefix}'

rule calc_stats:
    input:
        expand('pre-analysis/{sample}/spiker/split.report.txt',sample=samples)
    output:
        'pre-analysis/normalization.txt'
    threads:
        1
    script:
        '../scripts/calculate_norms.py'

rule normalize_bams:
    input:
        norm = 'pre-analysis/normalization.txt',
        bam_file = 'pre-analysis/{sample}/spiker/split_human.sorted.bam',
        inp_files = expand('pre-analysis/{sample}/spiker/split_human.sorted.bam',sample=samples)
    output:
        'pre-analysis/{sample}/normalized_files/{sample}.treat.pileup.SpikeIn_scaled.bigWig'
    threads:
        4
    conda:
        'spiker'
    script:
        '../scripts/normalize_bams.py'

rule clean_bigwigs:
    input:
        expand('pre-analysis/{sample}/normalized_files/{sample}.treat.pileup.SpikeIn_scaled.bigWig',sample=non_gst)
    output:
        directory('pre-analysis/normalized_bigwigs/')
    params:
        dir = '/media/asangani1/CDK7_project/2023-11-28/',
        norm_dir ='pre-analysis/normalized_files/'
    threads:
        1
    script:
        '../scripts/clean_bigwigs_spikein.py'

rule assembleStats_dros:
    priority: 10
    
    output:
        "pre-analysis/alignment_stats_dros.csv"
    params:
        files = "pre-analysis/*",
        dup_files = "pre-analysis/*/bowtie2_dros/dupMetrics.tsv",
        rmdup_files = "pre-analysis/*/bowtie2_dros/dupMetricsFiltered.tsv"

    script:
        "../scripts/bowtie_alignment_stats.R"

        
