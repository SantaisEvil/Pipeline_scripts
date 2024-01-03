samples,reads = glob_wildcards("data/{sample}_{read}.fastq.gz")
samples = list(set(samples))
non_gst = [x for x in samples if 'gst' not in x.lower()]
ctrls = [x for x in samples if 'dmso' in x.lower()]
exps = [x for x in samples if 'ct' in x.lower()]

### Reference point analysis

# rule genome_coding_profile:
#     input:
#         'pre-analysis/ucsc_norm/' + config["ref"]["build"] + '/{sample}.bw'
#     output:
#         'analysis/matrices/{sample}.gz'
#     threads:
#         8
#     params:
#         coding_genes = '/media/asangani2/scripts/workflows/resources/HG38/coding_genesUCSC.bed',
#         ref = 'TSS',
#         a = 2000,
#         b = 2000,
#         bs = 25
#     conda:
#         'ngs'
#     shell:
#         'computeMatrix reference-point --referencePoint {params.ref} --binSize {params.bs} -S {input} -R {params.coding_genes} -o {output} -a {params.a} -b {params.b} -p {threads}' 

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
#     shell:25
#         'plotProfile -m {input} -o {output}'

### Scaled regions analysis

rule genome_coding_profile_scale:
    input:
        'pre-analysis/ucsc_compare_ct/' + config["ref"]["build"] + '/{sample}.bw'
    output:
        'analysis/matrices/{sample}.scaled.gz'
    threads:
        8
    params:
        coding_genes = '/media/asangani2/resources/hg38/coding_genes_ucsc/coding_all.bed',
        ref = 1000,
        a = 1000,
        b = 1000,
        bs = 50
    conda:
        'ngs'
    shell:
        'computeMatrix scale-regions -S {input} --binSize {params.bs} -R {params.coding_genes} -o {output} -a {params.a} -b {params.b} -p {threads} --regionBodyLength {params.ref}' 

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
        sample=$(cut -d '_' -f1 <<<$sample)
        echo $sample
        computeMatrixOperations relabel -m {input} -o {output} --groupLabels $sample
        """

rule rbind_matrices_scale:
    input:
        expand('analysis/matrices/{sample}.scaled.relabeled.gz',sample=non_gst)
    output:
        'analysis/matrices/combined.scaled.gz'
    threads:
        1
    conda:
        'ngs'
    shell:
        """
        computeMatrixOperations rbind -m {input} -o {output} 
        """

rule plot_profiles_scale:
    input:
        'analysis/matrices/combined.scaled.gz'
    output:
        'analysis/plot_profile.scaled.png'
    threads:
        1
    params:
        opt = '--perGroup'
    conda:
        'ngs'
    shell:
        'plotProfile -m {input} -o {output}'