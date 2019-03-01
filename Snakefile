

rule all:
    input: "QC/multiqc_report.html", expand("salmon/{condition}/quant.sf", condition=["1","2"])


include: "rules/common.smk"
include: "rules/getData.smk"
include: "rules/qc.smk"
include: "rules/salmon.smk"
