samples,reads = glob_wildcards("data/{sample}_{read}.fastq.gz")
samples = list(set(samples))
non_gst = [x for x in samples if 'gst' not in x.lower()]
ctrls = [x for x in samples if 'dmso' in x.lower()]
exps = [x for x in samples if 'ct' in x.lower()]
ar = [x for x in samples if ('ar' in x.lower() or 'p300_wt' in x.lower()) and 'ar_tm' not in x.lower()]
phf19 = [x for x in samples if '_EZH2_' in x]
c9_samples = [x for x in samples if 'BrA' in x or 'Cas9' in x]
h3k4me1 = [x for x in c9_samples if 'h3k4me1' in x.lower()]
h3k27me3 = [x for x in c9_samples if 'h3k27me3' in x.lower()]
h3k36me2 = [x for x in c9_samples if 'h3k36me2' in x.lower()]

### Reference point analysis

# rule genome_coding_profile:
#     input:
#         'pre-analysis/ucsc_sam_view_cov/' + config["ref"]["build"] + '/{sample}.bw'
#     output:
#         'analysis/matrices/{sample}.gz'
#     threads:
#         8
#     params:
#         coding_genes = '/media/asangani2/scripts/workflows/resources/HG38/coding_genesUCSC.bed',
#         #coding_genes = 'pre-analysis/LNCaP_AR_Y1092A/macs/narrow/NA_peaks.narrowPeakFiltered',
#         ref = 'TSS',
#         a = 3000,
#         b = 3000,
#         bs = 50
#     conda:
#         'ngs'
#     shell:
#         'computeMatrix reference-point --referencePoint {params.ref} --binSize {params.bs} -S {input} -R {params.coding_genes} -o {output} -a {params.a} -b {params.b} -p {threads} --skipZeros' 

# rule relabel_matrices:
#     input:
#         'analysis/matrices/{sample}.gz'
#     output:
#         'analysis/matrices/{sample}.relabeled.gz'
#     threads:
#         1
#     conda:
#         'ngs'
#     shell:
#         """
#         sample=$(cut -d '/' -f3 <<< {input})
#         sample=$(cut -d '_' -f1 <<<$sample)
#         echo $sample
#         computeMatrixOperations relabel -m {input} -o {output} --groupLabels $sample
#         """

# rule rbind_matrices:
#     input:
#         expand('analysis/matrices/{sample}.relabeled.gz',sample=non_gst)
#     output:
#         'analysis/matrices/combined.gz'
#     threads:
#         1
#     conda:
#         'ngs'
#     shell:
#         """
#         computeMatrixOperations rbind -m {input} -o {output} 
#         """

# rule plot_profiles:
#     input:
#         'analysis/matrices/combined.gz'
#         #expand('analysis/matrices/{sample}.gz',sample=samples[0])
#     output:
#         'analysis/plot_profile.png'
#     threads:
#         1
#     params:
#         opt = '--perGroup'
#     conda:
#         'ngs'
#     shell:
#         'plotProfile -m {input} -o {output}'

### Scaled regions analysis

rule genome_coding_profile_scale:
    input:
        'pre-analysis/ucsc_norm_2/' + config["ref"]["build"] + '/{sample}.bw'
    output:
        'analysis/matrices/{sample}.scaled.gz'
    threads:
        8
    params:
        #coding_genes = '/media/asangani2/resources/hg38/coding_genes_ucsc/coding_all.bed',
        coding_genes = '/media/asangani2/resources/custom_bed_files/1K_NSD2_INdependent_ARsites.bed',
        #coding_genes = '~/landing/sknmc_phf19ko_up.bed',
        ref = 2000,
        a = 0,
        b = 0,
        bs = 10
    conda:
        'ngs'
    shell:
        'computeMatrix scale-regions -S {input} --binSize {params.bs} -R {params.coding_genes} -o {output} -a {params.a} -b {params.b} -p {threads} --regionBodyLength {params.ref} --skipZeros' 

rule relabel_matrices_scale:
    input:
        'analysis/matrices/{sample}.scaled.gz'
    output:
        'analysis/matrices/{sample}.scaled.relabeled.gz'
    threads:
        1
    conda:
        'ngs'
    shell:
        """
        sample=$(cut -d '/' -f3 <<< {input})
        sample=$(cut -d '-' -f1 <<<$sample)
        echo $sample
        computeMatrixOperations relabel -m {input} -o {output} --groupLabels $sample
        """

rule h3k27me3:
    input:
        expand('analysis/matrices/{sample}.scaled.relabeled.gz',sample=h3k27me3)
    output:
        'analysis/matrices/h3k27me3.scaled.gz'
    threads:
        1
    conda:
        'ngs'
    shell:
        """
        computeMatrixOperations rbind -m {input} -o {output} 
        """

rule rbind_matrices_scale_h3k4me1:
    input:
        expand('analysis/matrices/{sample}.scaled.relabeled.gz',sample=h3k4me1)
    output:
        'analysis/matrices/h3k4me1.scaled.gz'
    threads:
        1
    conda:
        'ngs'
    shell:
        """
        computeMatrixOperations rbind -m {input} -o {output} 
        """

rule rbind_matrices_scale_h3k36me2:
    input:
        expand('analysis/matrices/{sample}.scaled.relabeled.gz',sample=h3k36me2)
    output:
        'analysis/matrices/h3k36me2.scaled.gz'
    threads:
        1
    conda:
        'ngs'
    shell:
        """
        computeMatrixOperations rbind -m {input} -o {output} 
        """

rule plot_profiles_scale_1:
    input:
        'analysis/matrices/h3k4me1.scaled.gz'
    output:
        'analysis/h3k4me1.png'
    threads:
        1
    params:
        opt = '--perGroup'
    conda:
        'ngs'
    shell:
        'plotProfile -m {input} -o {output} --outFileNameData analysis/h3k4me1.tsv'

rule plot_profiles_scale:
    input:
        'analysis/matrices/h3k27me3.scaled.gz'
    output:
        'analysis/h3k27me3.png'
    threads:
        1
    params:
        opt = '--perGroup'
    conda:
        'ngs'
    shell:
        'plotProfile -m {input} -o {output} --outFileNameData analysis/h3k27me3.tsv'

rule plot_profiles_scale_2:
    input:
        'analysis/matrices/h3k36me2.scaled.gz'
    output:
        'analysis/h3k36me2.png'
    threads:
        1
    params:
        opt = '--perGroup'
    conda:
        'ngs'
    shell:
        'plotProfile -m {input} -o {output} --outFileNameData analysis/h3k36me2.tsv'