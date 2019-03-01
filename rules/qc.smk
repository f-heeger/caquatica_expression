rule fastqc:
    input: "reads/{sample}.fastq.gz" 
    output: "QC/{sample}_fastqc.html"
    threads: 6
    conda:
        "../envs/fastqc.yaml"
    shell:
        "fastqc -o QC -t {threads} {input}" % config

rule multiqc:
    input: expand("QC/{sample}_fastqc.html", sample=samples.keys())
    output: "QC/multiqc_report.html", "QC/multiqc_data/multiqc_fastqc.txt"
    conda:
        "../envs/multiqc.yaml"
    shell:
        "multiqc -f --interactive -o QC QC/*_fastqc.zip" % config
