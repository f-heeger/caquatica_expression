log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("tximport")
library("DESeq2")
library("readr")


samples = read.table(snakemake@input[[1]], header=F, sep="\t")
colnames(samples) = c("sId", "sName", "sraId", "condNum", "culture", "medium", "phase")

samples$medium = relevel(samples$medium, ref="malt")
samples$culture = relevel(samples$culture, ref="lq")
samples$phase = relevel(samples$phase, ref="exp")

if (snakemake@wildcards[["model"]] == "culture") {
        tSamples = subset(samples, medium=="straw")
#    } else if (snakemake@wildcards[["model"]] == "medium") {
#        tSamples = subset(samples, culture=="lq" & phase=="exp")
#    } else if (snakemake@wildcards[["model"]] == "phase") {
#        tSamples = subset(samples, culture=="lq" & medium=="straw")
    } else if (snakemake@wildcards[["model"]] == "all") {
        tSamples = subset(samples, culture=="lq")
    } else {
        stop(paste("Unknown model: ", snakemake@wildcards[["model"]]))
}


names = tSamples$sId
files = paste("salmon/", names, "/quant.sf", sep="")
names(files) = names

#snakemake@input[[seq(3:length(snakemake@input)]]

tx2gene = read_csv(snakemake@input[[2]])

txi = tximport(files, type="salmon", tx2gene=tx2gene)

if (snakemake@wildcards[["model"]] == "culture") {
    formStr ="~ culture"
} else {
    formStr ="~ medium + phase"
}

dds <- DESeqDataSetFromTximport(txi, colData = tSamples,
                        design = as.formula(formStr))
dds=DESeq(dds)

saveRDS(dds, file=snakemake@output[[1]])

png(snakemake@output[[2]])
plotDispEsts(dds)
dev.off()
