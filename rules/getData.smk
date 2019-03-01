rule getReads:
    output: "reads/{sampleId}.fastq.gz"
    log: "logs/{sampleId}_download.log"
    params: sraId=lambda wildcards: samples[wildcards.sampleId]["sraId"]
    conda:
        "../envs/sra.yaml"
    shell:
        "fastq-dump -X 10000 --gzip -Z {params.sraId} > {output} 2> {log}"
        #TODO remove maximum
