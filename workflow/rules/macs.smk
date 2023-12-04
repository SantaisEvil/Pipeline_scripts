# ruleorder: macs_narrow > macs_broad
rule macs_narrow:
    input:
        "pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam"
    
    output:
        multiext("pre-analysis/{sample}/macs/narrow/NA_","peaks.narrowPeak","summits.bed","peaks.xls")
    
    params:
        macs="callpeak -q 0.05 -f BAMPE -g hs -c pre-analysis/IgG/bowtie2/aligned.primary.rmdup.bam",
        outdir="pre-analysis/{sample}/macs/narrow",
        blacklist= config["ref"]["blacklist"]        

    threads: 6

    log:
        "pre-analysis/{sample}/logs/macs3.log"


    conda:
        "macs3"

    shell:
        """
        macs3 {params.macs}  -t {input} --outdir {params.outdir}  2> {log} 
        sed -i -r '/(chrX|chrM|chrUn)/d' {params.outdir}/NA_peaks.narrowPeak
        bedtools intersect -a {params.outdir}/NA_peaks.narrowPeak -b {params.blacklist} -v > {params.outdir}/NA_peaks.narrowPeakFiltered
        """



rule macs_broad:
    input:
        "pre-analysis/{sample}/bowtie2/aligned.primary.rmdup.bam"

    output:
        multiext("pre-analysis/{sample}/macs/broad/NA_","peaks.broadPeak","peaks.gappedPeak","peaks.xls")

    params:
        macs="callpeak --broad -g hs --broad-cutoff 0.1 -f BAMPE -c pre-analysis/IgG/bowtie2/aligned.primary.rmdup.bam",
        outdir="pre-analysis/{sample}/macs/broad",
        blacklist= config["ref"]["blacklist"]

    threads: 6

    log:
        "pre-analysis/{sample}/logs/macs3.log"

    conda:
        "macs3"

    shell:
        """
        macs3 {params.macs}  -t {input} --outdir {params.outdir}  2> {log}
        sed -i -r '/(chrX|chrM|chrUn)/d' {params.outdir}/NA_peaks.broadPeak
        bedtools intersect -a {params.outdir}/NA_peaks.broadPeak -b {params.blacklist} -v > {params.outdir}/NA_peaks.broadPeakFiltered

        """

#Automate mode (Narrow/Broad) selection based on Encode and our previous findings.
def macs_output(Sample):
    modes = []
    for index,i in enumerate(set(Sample)):
        if i == "igG":
            pass
        else:
            for mark in config["macs"]["broad"]:
                
                if mark in i.upper():
                    output = multiext(f"pre-analysis/{i}/macs/broad/NA_","peaks.broadPeak","peaks.gappedPeak","peaks.xls")
                    mode = "broad"
                    
                    break
                else:
                    output =  multiext(f"pre-analysis/{i}/macs/narrow/NA_","peaks.narrowPeak","summits.bed","peaks.xls")
                    mode = "narrow"
            modes.append(output)
            print(f"{i} is {mode}")
    print(modes)
    return(modes)





