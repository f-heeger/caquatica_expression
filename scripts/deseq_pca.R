log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library("DESeq2")
library("ggplot2")

dds <- readRDS(snakemake@input[[1]])
vsd <- varianceStabilizingTransformation(dds, blind=FALSE)

if (snakemake@wildcards[["model"]] == "all") {
    d=plotPCA(vsd, intgroup=c("medium", "phase"), returnData=T)
    ggplot(d, aes(x=PC1, y=PC2, shape=phase, color=medium)) + geom_point()
} else if (snakemake@wildcards[["model"]] == "culture") {
    d=plotPCA(vsd, intgroup=c("culture", "phase"), returnData=T)
    ggplot(d, aes(x=PC1, y=PC2, shape=phase, color=culture)) + geom_point()
} else {
        stop(paste("Unknown model: ", snakemake@wildcards[["model"]]))
}
ggsave(snakemake@output[[1]])


