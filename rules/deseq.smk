
rule trans2gene:
    input: "Clavariopsis_aquatica_WDA-00-1.transcripts.fa"
    output: "tanscript2gene.csv"
    run:
        with open(output[0], "w") as out:
            for line in open(input[0]):
                if line[0] == ">":
                    out.write(",".join(line.strip().split(" ")))

rule deseq_load:
    input:
        quant = 'salmon/{condition}/quant.sf'
    output:
        
    log:
        'logs/salmon/{condition}.log'
    params:
    script:
        "scripts/deseq_load.R"
