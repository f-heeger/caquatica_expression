export PATH=/home/heeger/anaconda3/bin:$PATH

#conda create -c bioconda -n snakemake5.4.2 snakemake=5.4.2

source activate snakemake5.4.2

snakemake --use-conda $@

conda deactivate
