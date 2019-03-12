

rule salmon_index:
    input:
        "Clavariopsis_aquatica_WDA-00-1.transcripts.fa"
    output:
        directory("salmon/transcriptome_index")
    log:
        "logs/salmon/transcriptome_index.log"
    threads: 2
    params:
        # optional parameters
        extra=""
    wrapper:
        "0.31.1/bio/salmon/index"

rule salmon_quant_reads:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        r="reads/{sample}.fastq.gz",
        index="salmon/transcriptome_index"
    output:
        quant = 'salmon/{sample}/quant.sf',
        lib = 'salmon/{sample}/lib_format_counts.json'
    log:
        'logs/salmon/{sample}.log'
    params:
        # optional parameters
        libtype ="A",
        zip_ext = "gz", # req'd for bz2 files ('bz2'); optional for gz files('gz')
        extra="--gcBias"
    threads: 2
    wrapper:
        "0.31.1/bio/salmon/quant"
