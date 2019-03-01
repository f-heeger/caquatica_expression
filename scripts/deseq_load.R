library("tximport")
library("DESeq2")
#library("readr")


samples = read.table("samples.tsv")
colnames(samples) = c("sId", "sName", "sraId", "condNum", "culture", "medium", "phase")

samples$medium = relevel(samples$medium, ref="malt")
samples$culture = relevel(samples$culture, ref="lq")
samples$phase = relevel(samples$phase, ref="exp")


files = snakemake@input
names(files) = samples$run #!!!

tx2gene = read_csv() #!!!

txi = tximport(files, type="salmon", tx2gene=tx2gene)

ddsTxi <- DESeqDataSetFromTximport(txi,
                                   colData = samples,
                                   design = ~ condition)
