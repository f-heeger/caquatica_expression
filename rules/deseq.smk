
#adapted from Johannes Koesters https://github.com/snakemake-workflows/rna-seq-star-deseq2

rule trans2gene:
    input: "Clavariopsis_aquatica_WDA-00-1.transcripts.fa"
    output: "transcript2gene.csv"
    run:
        with open(output[0], "w") as out:
            out.write("TAXNAME,GENEID\n")
            for line in open(input[0]):
                if line[0] == ">":
                    out.write(",".join(line[1:].strip().split(" "))+"\n")

def deseq_load_input(wildcards):
    files = ["samples.tsv", "transcript2gene.csv"]
    if wildcards.model == "medium":
        for row in samples.values():
            if row["culture"] == "lq" and row["phase"] == "exp":
                files.append("salmon/%(ID)s/quant.sf" % row)
    elif wildcards.model == "phase":
        for row in samples.values():
            if row["culture"] == "lq" and row["medium"] == "straw":
                files.append("salmon/%(ID)s/quant.sf" % row)
    elif wildcards.model == "culture":
        for row in samples.values():
            if row["medium"] == "straw":
                files.append("salmon/%(ID)s/quant.sf" % row)
    elif wildcards.model == "all":
        for row in samples.values():
            files.append("salmon/%(ID)s/quant.sf" % row)
    else:
        print("ERROR: Unknown model '%s'" % wildcards.model)
        return None
    return files
    
rule deseq_load:
    input:
        deseq_load_input
    output:
        "deseq2/{model}.rds", "deseq2/{model}_disp.png"
    log:
        'logs/deseq_load/{model}.log'
    conda:
        "../envs/deseq.yaml"
    script:
        "../scripts/deseq_load.R"

rule deseq_pca:
    input:
        "deseq2/{model}.rds"
    output:
        "results/{model}_pca.pdf"
    conda:
        "../envs/deseq.yaml"
    log:
        'logs/deseq_pca/{model}.log'
    script:
        "../scripts/deseq_pca.R"

#rule deseq2:
#    input:
#        "deseq2/{model}.rds"
#    output:
#        disp_medium="disp_medium.png", diff_medium="diff_medium.pdf", sm_tab="expr/expr_straw-malt.tsv"
#        table=report("results/diffexp/{contrast}.diffexp.tsv", "../report/diffexp.rst"),
#        ma_plot=report("results/diffexp/{contrast}.ma-plot.svg", "../report/ma.rst"),
#    params:
#        contrast=get_contrast
#    conda:
#        "../envs/deseq.yaml"
#    log:
#        "logs/deseq2/{contrast}.diffexp.log"
#    threads: get_deseq2_threads
#    script:
#"../scripts/deseq2.R"
