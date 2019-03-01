

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

def salmon_quant_reads_input(wildcards):
    cSamples = []
    for sId, sData in samples.items():
        if sData["condNum"] == wildcards.condition:
            cSamples.append(sId)
    return {"r": ["reads/%s.fastq.gz" % s for s in  cSamples], "index": "salmon/transcriptome_index"}

rule salmon_quant_reads:
    input:
        # If you have multiple fastq files for a single sample (e.g. technical replicates)
        # use a list for r1 and r2.
        unpack(salmon_quant_reads_input)
    output:
        quant = 'salmon/{condition}/quant.sf',
        lib = 'salmon/{condition}/lib_format_counts.json'
    log:
        'logs/salmon/{condition}.log'
    params:
        # optional parameters
        libtype ="A",
        zip_ext = "gz", # req'd for bz2 files ('bz2'); optional for gz files('gz')
        extra="--gcBias"
    threads: 2
    wrapper:
        "0.31.1/bio/salmon/quant"
