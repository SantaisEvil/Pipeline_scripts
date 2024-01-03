samples,reads = glob_wildcards("data/{sample}_{read}.fastq.gz")
spiker_files = ['split_human.bam','split_exogenous.bam','split_both.bam','split_neither.bam','split.report.txt']
samples = list(set(samples))
non_gst = [x for x in samples if 'gst' not in x.lower()]


# rule bowtie2_PE_dros:
#     input:
#         r1 = "data_trimmed/{sample}_R1.atria.fastq.gz",
#         r2 = "data_trimmed/{sample}_R2.atria.fastq.gz"
        
        
#     output:
#         aligned_sam="pre-analysis/{sample}/bowtie2_dros/unfiltered_aligned.sam"
#     log:
#         alignment_stats = "pre-analysis/{sample}/bowtie2_dros/unfiltered_aligned.txt",
#         met = "pre-analysis/{sample}/logs/bowtie2.txt"
#     threads: 8
#     params:
#         bowtie2 = "bowtie2 --local --no-mixed --no-discordant ",
#         index = config["ref"]["bowtie2_spikein"]
#     conda:
#         "ngs"

#     shell:
#         """
#         {params.bowtie2} -p {threads} --met-file {log.met} -x {params.index} -1 {input.r1} -2 {input.r2} -S {output.aligned_sam}  &> {log.alignment_stats}
        
#         """


# rule markDup_dros:
#     input:
#         bamFile= "pre-analysis/{sample}/bowtie2_dros/Aligned.sortedByCoord.out.bam"     
#     output:
#         bamMarked= "pre-analysis/{sample}/bowtie2_dros/aligned_markDup.bam",
#         dupStats="pre-analysis/{sample}/bowtie2_dros/dupMetrics.tsv" 
#     params:
#         samtools="-ASO=coordinate --TAGGING_POLICY All"
#     threads: 4

#     log:
#         "pre-analysis/{sample}/logs/MarkDuplicates.log"

#     shell:
#         "gatk MarkDuplicates -I {input.bamFile} -O  {output.bamMarked} -M {output.dupStats} 2> {log}"


# rule markDupFiltered_dros:
#     input:
#         bamFile= "pre-analysis/{sample}/bowtie2_dros/aligned_markDup.bam"     
#     output:
#         bamMarked= "pre-analysis/{sample}/bowtie2_dros/aligned_markDupFiltered.bam",
#         dupStats="pre-analysis/{sample}/bowtie2_dros/dupMetricsFiltered.tsv" 
#     params:
#         samtools="-ASO=coordinate --TAGGING_POLICY All"
#     threads: 4

#     log:
#         "pre-analysis/{sample}/logs/MarkDuplicates.log"

#     shell:
#         "gatk MarkDuplicates -I {input.bamFile} -O  {output.bamMarked} -M {output.dupStats} 2> {log}"

# rule samtools_dros:
#     input:
#         "pre-analysis/{sample}/bowtie2_dros/unfiltered_aligned.sam"
        
#     output:
#         "pre-analysis/{sample}/bowtie2_dros/unfiltered_aligned_stats.tsv"
#     threads: 4
#     params:
#         bwa = "-M",
#         index = config["ref"]["bowtie2"]
#     conda:
#         "ngs"

#     shell:
#         """
#         samtools flagstat -O tsv -@ {threads} {input} > {output}
        
#         """

# rule samToBamSorted_dros:
#     input:
#         samFile="pre-analysis/{sample}/bowtie2_dros/unfiltered_aligned.sam"         
#     output:
#         bamFile= "pre-analysis/{sample}/bowtie2_dros/Aligned.sortedByCoord.out.bam"
#     params:
#         samtools="sort --write-index -O BAM"
#     threads: 4

#     conda:
#         "ngs"

#     shell:
#         "samtools {params.samtools} -@ {threads}  -o {output.bamFile}  {input.samFile} "

# rule filterBam_dros:
#     input:
#         "pre-analysis/{sample}/bowtie2_dros/aligned_markDup.bam"      
#     output:
#         filteredBam="pre-analysis/{sample}/bowtie2_dros/aligned.primary.rmdup.bam",
#         index="pre-analysis/{sample}/bowtie2_dros/aligned.primary.rmdup.bam.bai"
#     params:
#         filter="view -h -f bam --num-filter /3340",
#         index="index --show-progress "
#     threads: 4

#     conda:
#         "ngs"

#     shell:
#         """
#         sambamba {params.filter} -t {threads} {input} > {output.filteredBam}
#         sambamba {params.index} -t {threads} {output.filteredBam} {output.index}
        
#         """

# rule bamToBw_dros:
#     input:
#         "pre-analysis/{sample}/bowtie2_dros/aligned.primary.rmdup.bam"      
#     output:
#         "pre-analysis/ucsc/" + config["ref"]["build"] + "/" + "{sample}.bw"
#     params:
#         main="--binSize 1 --normalizeUsing CPM --centerReads",
#         blacklist= config["ref"]["blacklist"]
#     threads: 4

#     log:
#         "pre-analysis/{sample}/logs/BigWig.log"

#     conda:
#         "ngs"

#     shell:
#         "bamCoverage {params.main} -bl {params.blacklist} -b {input} -o {output}  -p {threads}  > {log} "

# rule split_bam:
#     input:
#         'pre-analysis/{sample}/bowtie2_dros/aligned.primary.rmdup.bam'
#     output:
#         'pre-analysis/{sample}/spiker/split.report.txt'
#     threads:
#         4
#     params:
#         prefix = 'pre-analysis/{sample}/spiker/split'
#     conda:
#         'spiker'
#     shell:
#         'split_bam.py --threads {threads} -i {input} -o {params.prefix}'

rule calc_stats:
    input:
        expand('pre-analysis2/{sample}/spiker/split.report.txt',sample=samples)
    output:
        'pre-analysis2/normalization.txt'
    threads:
        1
    script:
        '../scripts/calculate_norms.py'

# rule normalize_bams:
#     input:
#         norm = 'pre-analysis/normalization.txt',
#         bam_file = 'pre-analysis/{sample}/spiker/split_human.sorted.bam',
#         inp_files = expand('pre-analysis/{sample}/spiker/split_human.sorted.bam',sample=samples)
#     output:
#         'pre-analysis/{sample}/normalized_files/{sample}.treat.pileup.SpikeIn_scaled.bigWig'
#     threads:
#         4
#     conda:
#         'spiker'
#     script:
#         '../scripts/normalize_bams.py'

# rule clean_bigwigs:
#     input:
#         expand('pre-analysis/{sample}/normalized_files/{sample}.treat.pileup.SpikeIn_scaled.bigWig',sample=non_gst)
#     output:
#         directory('pre-analysis/normalized_bigwigs/')
#     params:
#         dir = '/media/asangani1/CDK7_project/2023-11-28/',
#         norm_dir ='pre-analysis/normalized_files/'
#     threads:
#         1
#     script:
#         '../scripts/clean_bigwigs_spikein.py'

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

        
# rule samtools_normalize:
#     input:
#         norm_file = 'pre-analysis/normalization.txt',
#         bam_file = 'pre-analysis/{sample}/spiker/split_human.sorted.bam',
#         align_file = 'pre-analysis/alignment_stats_dros.csv'
#     output:
#         'pre-analysis/{sample}/normalized_bams/normalized.bam'
#     threads:
#         4
#     conda:
#         'ngs'
#     script:
#         '../scripts/samtools_norm.py'

rule norm_bigwigs:
    input:
        norm_file = 'pre-analysis/normalization.txt',
        bam_file = "pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam"      
    output:
        "pre-analysis/ucsc_norm/" + config["ref"]["build"] + "/" + "{sample}.bw"
    params:
        blacklist= config["ref"]["blacklist"]
    threads: 4
    log:
        "pre-analysis/{sample}/logs/BigWig_norm.log"

    conda:
        "ngs"

    script:
        '../scripts/samtools_norm.py'

# rule bamcompare_ct:
#     input:
#         bam_file = 'pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam',
#         norm_file = 'pre-analysis/normalization.txt'
#     output:
#         'pre-analysis/ucsc_compare_ct/' + config["ref"]["build"] + "/" + "{sample}.bw"
#     params:
#         ctrl_bam = 'pre-analysis/CT_RNAseH1_GST/bowtie2/aligned.primary.rmdup.bam',
#         blacklist= config["ref"]["blacklist"]
#     log:
#         "pre-analysis/{sample}/logs/BigWig_norm.log"
#     threads:
#         8
#     conda:
#         'ngs'
#     script:
#         '../scripts/samtools_norm.py'

# rule bamcompare_dmso:
#     input:
#         bam_file = 'pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam',
#         norm_file = 'pre-analysis/normalization.txt'
#     output:
#         'pre-analysis/ucsc_compare_dmso/' + config["ref"]["build"] + "/" + "{sample}.bw"
#     params:
#         ctrl_bam = 'pre-analysis/CT_RNAseH1_GST/bowtie2/aligned.primary.rmdup.bam',
#         blacklist= config["ref"]["blacklist"]
#     log:
#         "pre-analysis/{sample}/logs/BigWig_norm.log"
#     threads:
#         8
#     conda:
#         'ngs'
#     script:
#         '../scripts/samtools_norm.py'

# rule make_hub_norm:   
#     output:
#         trackDb = "pre-analysis/ucsc_norm/" + config["ref"]["build"] + "/" + "trackDb.txt",
#         hub = "pre-analysis/ucsc_norm/hub.txt",
#         genomes = "pre-analysis/ucsc_norm/genomes.txt"
#     params:
#         bwDir= "pre-analysis/ucsc_norm/" + config["ref"]["build"],
#         build = config["ref"]["build"]
#     conda:
#         'analysis'
#     script:
#         "../scripts/make_trackDB.py"

def macs_output_spikein(Sample):
    modes = []
    for index,i in enumerate(set(Sample)):
        if i == "igG":
            pass
        else:
            for mark in config["macs"]["broad"]:
                if mark in i.upper():
                    if 'CT' in i.upper():
                        output = multiext(f"pre-analysis/{i}/macs/broad/CT_","peaks.broadPeak","peaks.gappedPeak","peaks.xls")
                    elif 'DMSO' in i.upper():
                        output = multiext(f"pre-analysis/{i}/macs/broad/DMSO_","peaks.broadPeak","peaks.gappedPeak","peaks.xls")
                    mode = "broad"
                    break
                else:
                    if 'CT' in i.upper():
                        output =  multiext(f"pre-analysis/{i}/macs/narrow/CT_","peaks.narrowPeak","summits.bed","peaks.xls")
                    elif 'DMSO' in i.upper():
                        output =  multiext(f"pre-analysis/{i}/macs/narrow/DMSO_","peaks.narrowPeak","summits.bed","peaks.xls")
                    mode = "narrow"
            modes.append(output)
            print(f"{i} is {mode}")
    print(modes)
    return(modes)

rule macs_broad_spikein_CT:
    input:
        "pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam"
    output:
        multiext("pre-analysis/{sample}/macs/broad/CT_","peaks.broadPeak","peaks.gappedPeak","peaks.xls")

    params:
        macs="callpeak --broad -g hs --broad-cutoff 0.1 -f AUTO -n CT -c /media/asangani1/CDK7_project/2023-11-28/pre-analysis/CT_RNAseH1_GST/bowtie2/aligned.primary.rmdup.bam",
        outdir="pre-analysis/{sample}/macs/broad",
        blacklist= config["ref"]["blacklist"]

    threads: 8

    log:
        "pre-analysis/{sample}/logs/macs3.log"

    conda:
        "ngs"

    shell:
        """
        macs3 {params.macs}  -t {input} --outdir {params.outdir}  2> {log}
        sed -i -r '/(chrX|chrM|chrUn)/d' {params.outdir}/CT_peaks.broadPeak
        bedtools intersect -a {params.outdir}/CT_peaks.broadPeak -b {params.blacklist} -v > {params.outdir}/CT_peaks.broadPeakFiltered

        """

rule macs_broad_spikein_DMSO:
    input:
        "pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam"
    output:
        multiext("pre-analysis/{sample}/macs/broad/DMSO_","peaks.broadPeak","peaks.gappedPeak","peaks.xls")

    params:
        macs="callpeak --broad -g hs --broad-cutoff 0.1 -f AUTO -n DMSO -c /media/asangani1/CDK7_project/2023-11-28/pre-analysis/DMSO_RNAseH1_GST/bowtie2/aligned.primary.rmdup.bam",
        outdir="pre-analysis/{sample}/macs/broad",
        blacklist= config["ref"]["blacklist"]

    threads: 8

    log:
        "pre-analysis/{sample}/logs/macs3.log"

    conda:
        "ngs"

    shell:
        """
        macs3 {params.macs}  -t {input} --outdir {params.outdir}  2> {log}
        sed -i -r '/(chrX|chrM|chrUn)/d' {params.outdir}/DMSO_peaks.broadPeak
        bedtools intersect -a {params.outdir}/DMSO_peaks.broadPeak -b {params.blacklist} -v > {params.outdir}/DMSO_peaks.broadPeakFiltered
        """

def macs_spikein:
    input:
        "pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam"
    output:
        