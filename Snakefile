

rule all:
    input: "QC/multiqc_report.html", expand("results/{model}_pca.pdf", model=["culture","all"])


include: "rules/common.smk"
include: "rules/getData.smk"
include: "rules/qc.smk"
include: "rules/salmon.smk"
include: "rules/deseq.smk"
