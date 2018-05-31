import codecs
from urllib.request import urlopen
import math
import pickle
import json
from html.parser import HTMLParser
from itertools import zip_longest

from snakemake.utils import min_version, R
from Bio import SeqIO
from Bio.KEGG.REST import *
from Bio.KEGG.KGML import KGML_parser
from Bio.Graphics.KGML_vis import KGMLCanvas

from GoTools import parseGoOboDict
from KeggTools import CachedKeggKoToPathwayMap
from KeggTools import CachedKeggKoToEnzymeMap
from KeggTools import CachedKeggPathwayIdToNameMap
from KeggTools import KeggReactionIdToEcMap
from MultiCachedDict import PersistantDict

configfile: "config.json"

rule all:
    input: expand("activation/ca_active_{geneset}_{comp}.tsv", geneset=["go", "cazy"], comp=["straw-malt", "alder-straw", "alder-malt", "stat-exp", "solid-liquidExp", "solid-liquidSta"]), expand("ca_{gSet}_activation.tsv", gSet=["go", "cazy", "kegg"]), expand("expr/{dir}_{comp}.tsv", dir=["up", "down", "diff"], comp=["straw-malt", "alder-straw", "alder-malt", "stat-exp", "solid-liquidExp", "solid-liquidSta"]), expand("expr/diffAnno_{comp}.tsv", comp=["straw-malt", "alder-straw", "alder-malt", "stat-exp", "solid-liquidExp", "solid-liquidSta"]), "expr/expr_strawSS-malt.tsv", "expr/expr_alder-strawSS.tsv", expand("keggPlots/{comp}.ca_active_kegg.expressionplots", comp=["straw-malt", "alder-straw", "alder-malt", "stat-exp", "solid-liquidExp", "solid-liquidSta"]), "overlap/genes.tsv", "overlap/geneSets.tsv"

rule fastqc:
    input: "../ClaaquEProfiling/Raw Data/{sample}.fastq.gz" 
    output: "QC/{sample}_fastqc.html"
    threads: 6
    shell: "%(fastqc)s -o QC -t {threads} \"{input}\"" % config

rule multiqc:
    input: expand("QC/{sample}_fastqc.html", sample=config["fileNames"].values())
    output: "QC/multiqc_report.html", "QC/multiqc_data/multiqc_fastqc.txt"
    shell:
        "%(multiqc)s -f --interactive -o QC QC/*_fastqc.zip" % config

rule combineIntoMatrix:
    input: expand("../rsem2/{sample}.genes.results", sample=config["sampleNames"].values())
    output: fList="rsem.resultFiles.txt", mtx="genes.counts.matrix"
    shell:
        "ls {input} > {output.fList}; %(rsem2matrix)s --rsem_files {output.fList} --mode counts > {output.mtx}" % config

rule deSeq2:
    input: counts="genes.counts.matrix", samples="samples.tsv"
    output: pcas="pcas.pdf", disp_medium="disp_medium.png", disp_phase="disp_phase.png", disp_cultureE="disp_cultureExp.png", disp_cultureS="disp_cultureSta.png", diff_medium="diff_medium.pdf", diff_phase="diff_phase.pdf", diff_cultureE="diff_cultureExp.pdf", diff_cultureS="diff_cultureSta.pdf", am_tab="expr/expr_alder-malt.tsv", sm_tab="expr/expr_straw-malt.tsv", as_tab="expr/expr_alder-straw.tsv", sm_tab_ss="expr/expr_strawSS-malt.tsv", as_tab_ss="expr/expr_alder-strawSS.tsv", se_tab="expr/expr_stat-exp.tsv", sle_tab="expr/expr_solid-liquidExp.tsv", sls_tab="expr/expr_solid-liquidSta.tsv"
    run:
        R("""
        library(DESeq2, quietly=T)
        library(ggplot2, quietly=T)

        colData <- read.table("{input.samples}", row.names=1, header=T)
        colData$medium <- relevel(colData$medium, ref="malt")
        colData$culture <- relevel(colData$culture, ref="lq")
        colData$phase <- relevel(colData$phase, ref="exp")

        rawCountData <- read.table("{input.counts}", row.names=1, header=T)
        intCountData = sapply(rawCountData, as.integer)
        countData = intCountData[,rownames(colData)]
        
        row.names(countData) = row.names(rawCountData)
        

        dds <- DESeqDataSetFromMatrix(countData = countData, colData = colData, design = ~ culture+medium)

        dds_phase <- DESeqDataSetFromMatrix(countData = countData[,c(1,2,3,13,14)], 
                                            colData = data.frame(phase=colData[c(1,2,3,13,14), 3]), 
                                            design = ~ phase)

        dds_culture_exp <- DESeqDataSetFromMatrix(countData = countData[,c(10,11,12,13,14)], 
                                                  colData = data.frame(culture=colData[c(10,11,12,13,14), 1]), 
                                              design = ~ culture)
        dds_culture_sta <- DESeqDataSetFromMatrix(countData = countData[,c(1,2,3,10,11,12)], 
                                                  colData = data.frame(culture=colData[c(1,2,3,10,11,12), 1]), 
                                              design = ~ culture)

        dds_lq_exp <- DESeqDataSetFromMatrix(countData = countData[,c(4,5,6,7,8,9,13,14)],
                                             colData = data.frame(medium=colData[c(4,5,6,7,8,9,13,14), 2]),
                                             design = ~ medium)
        dds_lq_exp_ss <- DESeqDataSetFromMatrix(countData = countData[,c(4,5,6,7,8,9,15,16)],
                                             colData = data.frame(medium=colData[c(4,5,6,7,8,9,15,16), 2]),
                                             design = ~ medium)

        vsd <- varianceStabilizingTransformation(dds, blind=FALSE)
        d=plotPCA(vsd, intgroup=c("medium", "culture"), returnData=T)
        p1=ggplot(d, aes(x=PC1, y=PC2, shape=culture, color=medium)) + geom_point()
        vsd <- varianceStabilizingTransformation(dds_phase, blind=FALSE)
        p2=plotPCA(vsd, intgroup=c("phase"))
        vsd <- varianceStabilizingTransformation(dds_culture_exp, blind=FALSE)
        p3=plotPCA(vsd, intgroup=c("culture"))
        vsd <- varianceStabilizingTransformation(dds_culture_sta, blind=FALSE)
        p3=plotPCA(vsd, intgroup=c("culture"))
        vsd <- varianceStabilizingTransformation(dds_lq_exp, blind=FALSE)
        p4=plotPCA(vsd, intgroup=c("medium"))
        pdf("{output.pcas}")
        print(p1)
        print(p2)
        print(p3)
        print(p4)
        dev.off()

        #MEDIUM RESULT
        dds_lq_exp <- DESeq(dds_lq_exp, betaPrior=T)

        res_alder_malt = results(dds_lq_exp, alpha=0.05, contrast=c("medium","alder","malt"))
        res_straw_malt = results(dds_lq_exp, alpha=0.05, contrast=c("medium","straw","malt"))
        res_alder_straw = results(dds_lq_exp, alpha=0.05, contrast=c("medium","alder","straw"))


        png("{output.disp_medium}")
        plotDispEsts(dds_lq_exp)
        dev.off()
        pdf("{output.diff_medium}")
        plotMA(dds_lq_exp)
        dev.off()
        write.table(res_alder_malt, "{output.am_tab}", quote=F, sep="\t")
        write.table(res_straw_malt, "{output.sm_tab}", quote=F, sep="\t")
        write.table(res_alder_straw, "{output.as_tab}", quote=F, sep="\t")
        
        #SUB-SAMPLE MEDIUM RESULT
        dds_lq_exp_ss <- DESeq(dds_lq_exp_ss, betaPrior=T)

        res_straw_malt = results(dds_lq_exp_ss, alpha=0.05, contrast=c("medium","straw","malt"))
        res_alder_straw = results(dds_lq_exp_ss, alpha=0.05, contrast=c("medium","alder","straw"))


        write.table(res_straw_malt, "{output.sm_tab_ss}", quote=F, sep="\t")
        write.table(res_alder_straw, "{output.as_tab_ss}", quote=F, sep="\t")
        
        #PHASE RESULT
        dds_phase = DESeq(dds_phase, betaPrior=T)
        res_stat_exp = results(dds_phase, alpha=0.05)
        
        png("{output.disp_phase}")
        plotDispEsts(dds_phase)
        dev.off()
        pdf("{output.diff_phase}")
        plotMA(dds_phase)
        dev.off()
        write.table(res_stat_exp, "{output.se_tab}", quote=F, sep="\t")
        
        #CULTURE RESULT
        dds_culture_exp = DESeq(dds_culture_exp, betaPrior=T)
        res_solid_liquidExp = results(dds_culture_exp, alpha=0.05)
        
        png("{output.disp_cultureE}")
        plotDispEsts(dds_culture_exp)
        dev.off()
        pdf("{output.diff_cultureE}")
        plotMA(dds_culture_exp)
        dev.off()
        write.table(res_solid_liquidExp, "{output.sle_tab}", quote=F, sep="\t")
        
        dds_culture_sta = DESeq(dds_culture_sta, betaPrior=T)
        res_solid_liquidSta = results(dds_culture_sta, alpha=0.05)
        
        png("{output.disp_cultureS}")
        plotDispEsts(dds_culture_sta)
        dev.off()
        pdf("{output.diff_cultureS}")
        plotMA(dds_culture_sta)
        dev.off()
        write.table(res_solid_liquidSta, "{output.sls_tab}", quote=F, sep="\t")
        """)


rule diffExpr:
    input: "expr/expr_{comp}.tsv"
    output: diff="expr/diff_{comp}.tsv", up="expr/up_{comp}.tsv", down="expr/down_{comp}.tsv"
    run:
        with open(output.diff, "w") as diff, open(output.up, "w") as up, open(output.down, "w") as down:
            with open(input[0]) as inStream:
                next(inStream) #header
                for line in inStream:
                    gId, mean, change, lfc, start, p, padj = line.strip("\n").split("\t")
                    if padj != "NA" and float(padj)<0.05 and abs(float(change))>1:
                        diff.write(line)
                        if float(change)>0:
                            up.write(line)
                        else:
                            down.write(line)

rule subsampleRes:
    input: all="expr/expr_straw-malt.tsv", subset="expr/expr_strawSS-malt.tsv"
    output: "subsampleChange.txt"
    run:
        allUp = set()
        allDown = set()
        allNC = set()
        with open(input.all) as inStream:
            header = next(inStream)
            for line in inStream:
                pId, mean, change, se, start, pvalue, padj = line.strip("\n").split("\t")
                if padj != "NA" and float(padj)<0.05:
                    if float(change)>1:
                        allUp.add(pId)
                    elif float(change)<-1:
                        allDown.add(pId)
                    else:
                        allNC.add(pId)
                else:
                    allNC.add(pId)
        ssUp = set()
        ssDown = set()
        ssNC = set()
        with open(input.subset) as inStream:
            header = next(inStream)
            for line in inStream:
                pId, mean, change, se, start, pvalue, padj = line.strip("\n").split("\t")
                if padj != "NA" and float(padj)<0.05:
                    if float(change)>1:
                        ssUp.add(pId)
                    elif float(change)<-1:
                        ssDown.add(pId)
                    else:
                        ssNC.add(pId)
                else:
                    ssNC.add(pId)
        assert len(allUp.union(allDown).union(allNC)) == len(ssUp.union(ssDown).union(ssNC))
        
        with open(output[0], "w") as out:
            out.write("both up:\t%i\n" % len(allUp & ssUp))
            out.write("both down:\t%i\n" % len(allDown & ssDown))
            out.write("both no change:\t%i\n" % len(allNC & ssNC))
            out.write("total:\t%i\n" % len(allUp.union(allDown).union(allNC)))
            out.write("percentage same:\t%f%%\n" % ((len(allUp & ssUp)+len(allDown & ssDown)+len(allNC & ssNC))/len(allUp.union(allDown).union(allNC))*100))

rule signalP:
    input: "../funannotate/caqua/predict_results/Clavariopsis_aquatica_WDA-00-1.proteins.fa"
    output: "signalp/Clavariopsis_aquatica_WDA-00-1.signalp.out"
    shell:
        "%(signalp)s -t euk -f short {input} > {output}" % config

rule collectAnno:
    input: tsv="../Clavariopsis_aquatica_WDA-00-1.annotations.txt"
    output: go="ca_proteins_go.tsv", cog="ca_proteins_cog.tsv", cazy="ca_proteins_cazy.tsv"
    run:
        go = {}
        cazy = {}
        cog = {}
        with open(input.tsv) as inStream:
            header = next(inStream)
            arr = header.split("\t")
            assert arr[12] == "COG"
            assert arr[13] == "GO Terms"
            assert arr[17] == "CAZyme"
            for line in inStream:
                arr = line.split("\t")
                gId = arr[0]
                if arr[1] != "CDS":
                    continue
                cogIds = arr[12].split(";")
                goIds = arr[13].split(";")
                cazyIds = arr[17].split(";")
                for gId in goIds:
                    if len(gId) == 0:
                        continue
                    try:
                        go[gId].append(arr[0])
                    except KeyError:
                        go[gId] = [arr[0]]
                for cId in cazyIds:
                    if len(cId) == 0:
                        continue
                    try:
                        cazy[cId].append(arr[0])
                    except KeyError:
                        cazy[cId] = [arr[0]]
                for cId in cogIds:
                    if len(cId) == 0:
                        continue
                    try:
                        cog[cId].append(arr[0])
                    except KeyError:
                        cog[cId] = [arr[0]]
        
        with open(output.go, "w") as out:
            for gId, genes in go.items():
                out.write("%s\t%s\n" % (gId, ",".join(genes)))
        with open(output.cog, "w") as out:
            for cId, genes in cog.items():
                out.write("%s\t%s\n" % (cId, ",".join(genes)))
        with open(output.cazy, "w") as out:
            for cId, genes in cazy.items():
                out.write("%s\t%s\n" % (cId, ",".join(genes)))

rule getKeggKO:
    input: kegg="../kegg/user_ko.tsv"
    output: "ca_proteins_kegg.tsv"
    params: db="/home/heeger/data/dbs/keggKo2Pathway.db"
    run:
        ko2path = CachedKeggKoToPathwayMap(params.db)
        path = {}
        for line in open(input.kegg):
            arr = line.strip("\n").split("\t")
            if len(arr) == 1:
                continue
            gene, ko = arr
            try:
                pathSet = ko2path[ko]
            except KeyError:
                continue
            for p in pathSet:
                try:
                    path[p].append(gene.split("-")[0])
                except KeyError:
                    path[p] = [gene.split("-")[0]]
        with open(output[0], "w") as out:
            for p, gset in path.items():
                out.write("%s\t%s\n" % (p, ",".join(gset)))

rule anno:
    input: fun="../Clavariopsis_aquatica_WDA-00-1.annotations.txt", kegg="ca_proteins_kegg.tsv", signalp="signalp/Clavariopsis_aquatica_WDA-00-1.signalp.out"
    output: json="anno.json", pickle="anno.pic"
    run:
        go = {}
        cazy = {}
        interpro = {}
        name = {}
        kegg = {}
        sigP = {}

        for line in open(input.fun):
            gene, feature, contig, start, stop, strand, tName, product, busco, pfam, interPro, eggNog, cog, goTerms, secreted, membrane, protease, cazyme, notes, translation = line.split("\t")
            if gene not in go:
                go[gene] = set()
            if len(goTerms) >0:
                go[gene] |= set(goTerms.split(";"))
            if gene not in cazy:
                cazy[gene] = set()
                
            if len(cazyme) > 0:
                cazy[gene] |= set(cazyme.split(";"))
            if gene not in interpro:
                interpro[gene] = set()
            if len(interPro) > 0:
                interpro[gene] |= set(interPro.split(";"))
            name[gene] = tName
            kegg[gene] = []

        for line in open(input.kegg):
            path, genes = line.strip("\n").split("\t")
            for gene in genes.split(","):
                kegg[gene].append(path)

        for line in open(input.signalp):
            if line[0] == "#":
                continue
            arr = line.strip("\n").split()
            sigP[arr[0].split("-")[0]] = arr[9]

        anno = {"go": go, "cazy": cazy, "interpro": interpro, "name": name, "kegg": kegg, "signalp": sigP}

        pickle.dump(anno, open(output.pickle, "wb"))
        json.dump(anno, open(output.json, "wt"), indent=2)

rule mgsa:
    input: anno="ca_proteins_{geneset}.tsv", exp="expr/expr_{comp}.tsv"
    output: tab="activation/ca_active_{geneset}_{comp}.tsv"
    run:
        R("""
        library(mgsa, quietly=T)
        set.seed(42)
        tab=read.table("{input.anno}", sep="\t", stringsAsFactors=F, as.is=T)
        setNames = tab[,1]
        sets = strsplit(tab[,2], ',')
        names(sets) = setNames


        expTab=read.table("{input.exp}", sep="\t", row.names=NULL)
        observations = expTab$row.names[!is.na(expTab$padj) & expTab$padj<0.05 & abs(expTab$log2FoldChange)>1]
        print(length(observations))
        result=mgsa(observations, sets)
        outTab=cbind(rownames(setsResults(result)), setsResults(result))
        colnames(outTab) = c("setName", colnames(outTab)[2:5])
        write.table(outTab, "{output.tab}", row.names=F, sep="\t", quote=F)
        """)

rule getGoSummary:
    input: expand("activation/ca_active_{{goSet}}_{comp}.tsv", comp=["straw-malt", "alder-straw", "alder-malt", "stat-exp", "solid-liquidExp", "solid-liquidSta"])
    output: "ca_{goSet,go.*}_activation.tsv"
    run:
        reader = codecs.getreader("utf-8")
        go = parseGoOboDict(reader(urlopen("http://purl.obolibrary.org/obo/go.obo")))
        with open(output[0], "w") as out:
            out.write("condition\tGO ID\tGO name\tin Population\tin Study Set\tactivation\n")
            for inFile in input:
                comp = inFile.rsplit("_", 1)[-1].split(".")[0]
                with open(inFile) as inStream:
                    headerLine = next(inStream)
                    for line in inStream:
                        goId, inPop, inSet, activation, _ = line.strip().split("\t")
                        if float(activation) > 0.6:
                            out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (comp, goId, go[goId].name, inPop, inSet, activation))

rule getCazySummary:
    input: expand("activation/ca_active_cazy_{comp}.tsv", comp=["straw-malt", "alder-straw", "alder-malt", "stat-exp", "solid-liquidExp", "solid-liquidSta"])
    output: "ca_cazy_activation.tsv"
    run:
#        famName={}
#        reader = codecs.getreader("utf-8")
#        for line in reader(urlopen("http://csbl.bmb.uga.edu/dbCAN/download/FamInfo.txt")):
#            fam, ccd, cls, note, activ = line.strip("\n").split("\t")
#            famName[fam] = activ
        with open(output[0], "w") as out:
            out.write("condition\tCAZy Fam\tactivity\tin Population\tin Study Set\tactivation\n")
            for inFile in input:
                comp = inFile.rsplit("_", 1)[-1].split(".")[0]
                with open(inFile) as inStream:
                    headerLine = next(inStream)
                    for line in inStream:
                        fam, inPop, inSet, activation, _ = line.strip().split("\t")
                        if float(activation) > 0.6:
                            out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (comp, fam, "--", inPop, inSet, activation))

rule getEcNum:
    input: "../kegg/user_ko.tsv"
    output: "ca_proteins_ec.tsv"
    params: db="/home/heeger/data/dbs/KeggKoToEcMap.db"
    run:
        ko2ec = CachedKeggKoToEnzymeMap(params.db)
        with open(output[0], "w") as out:
            for line in open(input[0]):
                arr = line.strip("\n").split("\t")
                if len(arr) == 1:
                    continue
                gene, ko = arr
                try:
                    ecSet = ko2ec[ko]
                except KeyError:
                    continue
                out.write("%s\t%s\n" % (gene, ",".join(ecSet)))

rule getKeggSummary:
    input: expand("activation/ca_active_kegg_{comp}.tsv", comp=["straw-malt", "alder-straw", "alder-malt", "stat-exp", "solid-liquidExp", "solid-liquidSta"])
    output: "ca_kegg_activation.tsv"
    params: db="/home/heeger/data/dbs/KeggPathwayIdToNameMap.db"
    run:
        kegg = CachedKeggPathwayIdToNameMap(params.db)
        with open(output[0], "w") as out:
            out.write("condition\tKEGG ID\tKEGG pathway name\tin Population\tin Study Set\tactivation\n")
            for inFile in input:
                comp = inFile.rsplit("_", 1)[-1].split(".")[0]
                with open(inFile) as inStream:
                    headerLine = next(inStream)
                    for line in inStream:
                        keggId, inPop, inSet, activation, _ = line.strip().split("\t")
                        if float(activation) > 0.6:
                            out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (comp, keggId, kegg[keggId], inPop, inSet, activation))

rule KeggUserData:
    input: ko="../kegg/user_ko.tsv", exp="expr/expr_{comp}.tsv", geneset="ca_proteins_kegg.tsv", act="activation/ca_active_kegg_{comp}.tsv"
    output: dynamic("keggUserData/{comp}_{path}.userData.tsv")
    run:
        gId2ko={}
        for line in open(input.ko):
            arr = line.strip("\n").split("\t")
            if len(arr) > 1:
                tId, ko = arr
                gId2ko[tId.split("-")[0]] = ko
        expr={}
        diff={}
        with open(input.exp) as diffIn:
            header = next(diffIn)
            for line in diffIn:
                gId, mean, change, lfc, start, p, adjp = line.strip("\n").split("\t")
                if adjp != "NA":
                    expr[gId] = float(change)
                    if float(adjp)<0.05 and abs(float(change))>1:
                        diff[gId] = True
                    else:
                        diff[gId] = False
        active = set()
        with open(input.act) as inStream:
            header = next(inStream)
            for line in inStream:
                sId, _ = line.split("\t", 1)
                active.add(sId)
            
        gSet={}
        for line in open(input.geneset):
            path, geneList = line.strip().split("\t")
            #if path in active:
            gSet[path] = geneList.split(",")

        for path in gSet:
            if path in ["map01110", "map00958"]:
                print("skipping path: %s. Its too big." % path)
                continue
            koExpr = {}
            koDiff = {}
            for gId in gSet[path]:
                try:
                    tKo = gId2ko[gId]
                    tExpr = expr[gId]
                except KeyError:
                    continue
                try:
                    koExpr[tKo].append(tExpr)
                except KeyError:
                    koExpr[tKo] = [tExpr]
                try:
                    koDiff[tKo].append(diff[gId])
                except KeyError:
                    koDiff[tKo] = [diff[gId]]
            if len(koExpr) > 0:
                lim = math.ceil(max([abs(v) for vList in koExpr.values() for v in vList]))
                scale = 510.0/(lim*2)
                colors = {}
                border = {}
                for ko in koExpr:
                    if len(koExpr[ko]) > 1:
                        if all([e<0 for e in koExpr[ko]]) or all([e>0 for e in koExpr[ko]]):
                            colors[ko] = dezToColor(((sum(koExpr[ko])/len(koExpr[ko]))+lim)*scale)
                            if any(koDiff[ko]):
                                border[ko] = "black"
                            else:
                                border[ko] = "grey"
                        else:
                            colors[ko] = "blue"
                            border[ko] = "white"
                    else:
                        colors[ko] = dezToColor((koExpr[ko][0]+lim)*scale)
                        border[ko] = "grey"
                    if koDiff[ko][0]:
                        border[ko] = "black"
                with open("keggUserData/%s_%s.userData.tsv" % (wildcards.comp, path), "w") as out:
                    for ko, col in colors.items():
                        out.write("%s\t%s,%s\n" % (ko, col, border[ko]))

rule keggPlot:
    input: act="activation/ca_active_kegg_{comp}.tsv", exp="expr/expr_{comp}.tsv", geneset="ca_proteins_kegg.tsv", ec="ca_proteins_ec.tsv", ko="../kegg/user_ko.tsv"
    output: "keggPlots/{comp}.ca_active_kegg.expressionplots"
    log: "logs/keggPlot_{comp}.log"
    params: minProp = 0.6,
    run:
#        reac2ec = KeggReactionIdToEcMap(cachePath="keggReacToEc.tsv")
#        
#        gId2ec={}
#        for line in open(input.ec):
#            arr = line.strip().split("\t")
#            if len(arr) == 1:
#                continue
#            gId, ecNum = arr
#            gId2ec[gId.split("-")[0]] = ecNum
#            
        gId2ko={}
        for line in open(input.ko):
            arr = line.strip("\n").split("\t")
            if len(arr) > 1:
                tId, ko = arr
                gId2ko[tId.split("-")[0]] = ko
        prop={}
        with open(input.act) as actIn:
            header = next(actIn)
            for line in actIn:
                path, pop, stud, act, err = line.strip().split("\t")
                if float(act) > params.minProp:
                    prop[path] = float(act)
        gSet={}
        allGenes={}
        for line in open(input.geneset):
            path, geneList = line.strip().split("\t")
            gSet[path] = geneList.split(",")
            for gId in gSet[path]:
                allGenes[gId] = None
        expr={}
        diff={}
        with open(input.exp) as diffIn:
            header = next(diffIn)
            for line in diffIn:
                gId, mean, change, lfc, start, p, adjp = line.strip("\n").split("\t")
                if adjp != "NA":
                    if gId in allGenes:
                        expr[gId] = float(change)
                        if float(adjp)<0.05 and abs(float(change))>1:
                            diff[gId] = True
                        else:
                            diff[gId] = False
        with open(log[0], "w") as logStream:
            for path in gSet:
                if path in ["ko01110", "ko00958"]:
                    print("skipping path: %s. Its too big." % path)
                    continue
                if path not in prop and path not in ["ko00500", "ko04146"]:
                    #skip non-significant pathes
                    continue
                koExpr = {}
                koDiff = {}
                for gId in gSet[path]:
                    try:
                        tKo = gId2ko[gId]
                    except KeyError:
                        continue
                    try:
                        tExpr = expr[gId]
                    except KeyError:
                        logStream.write("No expressions for %s\n" % gId)
                        continue
                    if tKo not in koExpr:
                        koExpr[tKo] = []
                        koDiff[tKo] = []
                    koExpr[tKo].append(tExpr)
                    koDiff[tKo].append(diff[gId])
                with open("keggPlots/%s.ca_active_kegg.%s.tsv" % (wildcards.comp, path), "w") as out:
                    for gId in gSet[path]:
                        out.write("%s\t%s\t%s\n" % (gId, gId2ko[gId], expr.get(gId, "-")))
                
                if len(koExpr) == 0:
                    continue
                    #nothing we can do for this pathway
            
                pathway = KGML_parser.read(kegg_get(path.replace("map", "ko"), "kgml"))
                reacExpr={}
                reacDiff={}
                for reac in pathway.orthologs:
                    for entry in reac.name.split(" "):
                        rko = entry.split(":")[1]
                        if rko in koExpr:
                            try:
                                reacExpr[reac.id].extend(koExpr[rko])
                            except KeyError:
                                reacExpr[reac.id] = koExpr[rko]
                                reacDiff[reac.id] = []
                            reacDiff[reac.id].extend(koDiff[rko])
                colors = {}
                border = {}
                lim = math.ceil(max([abs(v) for vList in reacExpr.values() for v in vList]))
                scale = 510.0/(lim*2)
                for rId in reacExpr:
                    if len(reacExpr[rId]) > 1:
                        if all([e<0 for e in reacExpr[rId]]) or all([e>0 for e in reacExpr[rId]]):
                            colors[rId] = dezToColor(((sum(reacExpr[rId])/len(reacExpr[rId]))+lim)*scale)
                            if any(reacDiff[rId]):
                                border[rId] = "#000000"
                            else:
                                border[rId] = "#878787"
                        else:
                            colors[rId] = "#0000ff"
                            border[rId] = "#ffffff"
                    else:
                        colors[rId] = dezToColor((reacExpr[rId][0]+lim)*scale)
                        border[rId] = "#878787"
                    if reacDiff[rId][0]:
                        border[rId] = "#000000"
                for reac in pathway.orthologs:
                    if reac.id in colors:
                        reac.graphics[0].bgcolor = colors[reac.id]
                        reac.graphics[0].fgcolor = border[reac.id]
                            
                #remove compound labels
                for compound in pathway.compounds:
                    for graphic in compound.graphics:
            #            print("remove label: "+graphic.name)
                        graphic.name = ""
                canvas = KGMLCanvas(pathway, import_imagemap=True)
                canvas.fontsize=10
                outPath="keggPlots/%s.ca_active_kegg.%s.pdf" % (wildcards.comp, path)
                canvas.draw(outPath)
        open(output[0], "w")

rule collectExp:
    input: "singleGenePlots/{group}_geneList.tsv", expand("expr/expr_{comp}.tsv", comp=["straw-malt", "alder-malt"])
    output: "singleGenePlots/{group}_expression.tsv"
    run:
        gene_list = set()
        g2n = {}
        
        for line in open(input[0]):
            name, gene = line.strip().split("\t")
            gene_list.add(gene)
            g2n[gene] = name
        
        with open(output[0], "w") as out:
            out.write("gene\tname\tcomp\tchange\tpvalue\n")
            for inFile in input[1:]:
                comp = inFile.rsplit(".", 1)[0].rsplit("_", 1)[-1]
                with open(inFile) as inStream:
                    next(inStream) #header
                    for line in inStream:
                        gId, mean, change, lfc, start, p, padj = line.strip("\n").split("\t")
                        if gId in gene_list:
                            out.write("%s\t%s\t%s\t%s\t%s\n" % (gId, g2n[gId], comp, change, padj))

rule plotExp:
    input: "singleGenePlots/{group}_expression.tsv"
    output: "singleGenePlots/{group}_expression.pdf"
    run:
        R("""

        library(ggplot2)

        d=read.table("{input}", sep="\t", header=T)
        sig=rep("", length(d$pvalue))
        sig[d$pvalue<0.05] = "*"
        sig[d$pvalue<0.01] = "**"
        sig[d$pvalue<0.001] = "***" 
        d$sig = sig
        ggplot(d, aes(gene, change, fill=comp, label=sig)) + geom_bar(position="dodge", stat="identity", width=.6) + coord_flip() + geom_text(aes(label=sig), position=position_dodge(.6)) + ggtitle("{wildcards.group}") + scale_x_discrete(breaks=d$gene, labels=paste(d$gene, d$name, sep="\n"))


        #number of comparisons:
        c=length(levels(d$comp))
        #number of rows
        n=dim(d)[1]
        ggsave("{output}", width=7, height=max(7,n/c), limitsize = F)

        """)


rule plotExpOverview_data:
    input: expand("singleGenePlots/{group}_expression.tsv", group=["cazy-AA1", "peroxi", "cyp450", "cazy-AA9", "cazy-CE1", "cazy-GH5", "cazy-GH7", "cazy-GH10", "cazy-GH11", "kegg-ko00052", "kegg-ko00040", "kegg-ko00500", "kegg-ko00640", "kegg-ko04146"])
    output: "singleGenePlots/expOverview.tsv"
    run:
        num = {}
        order = []
        for inFile in input:
            group = inFile.split("/", 1)[1].rsplit("_", 1)[0]
            num[group] = {}
            order.append(group)
            for comp in ["alder-malt", "straw-malt"]:
                num[group][comp] = {"up": 0, "down": 0}
            with open(inFile) as inStream:
                header = next(inStream)
                for line in inStream:
                    gene, name, comp, change, pvalue = line.strip("\n").split("\t")
                    if change!="NA" and pvalue != "NA" and abs(float(change)) > 1 and float(pvalue)<0.05:
                        if float(change) > 0:
                            num[group][comp]["up"] += 1
                        else:
                            num[group][comp]["down"] -= 1
        with open(output[0], "w") as out:
            out.write("group\tcomp\tdirection\tvalue\n")
            for group in order:
                for comp, cDat in num[group].items():
                    for direc in ["up", "down"]:
                        out.write("%s\t%s\t%s\t%i\n" % (group, comp, direc, cDat[direc]))

rule plotExpOverview_plot:
    input: "singleGenePlots/expOverview.tsv"
    output: "singleGenePlots/expOverview.svg"
    run:
        R("""
        library(ggplot2)
        d = read.table("{input}", header=T)
        d$group = factor(d$group, levels=unique(d$group))
        ggplot(d) + geom_bar(aes(x=group, y=value, fill=comp), stat="identity", position="dodge") + geom_hline(yintercept = 0) + coord_flip() + scale_fill_manual(values=c("#377a3e", "#c4b211"))
        ggsave("{output}")
        """)

rule geneOverlap:
    input: stma="expr/expr_straw-malt.tsv", alma="expr/expr_alder-malt.tsv"
    output: "overlap/genes.tsv"
    run:
        gNr = 0
        stma = set()
        with open(input.stma) as inStream:
            next(inStream) #header
            for line in inStream:
                gId, mean, change, lfc, start, p, padj = line.strip("\n").split("\t")
                if padj != "NA" and float(padj)<0.05 and abs(float(change))>1:
                    stma.add(gId)
                gNr += 1
        alma = set()
        with open(input.alma) as inStream:
            next(inStream) #header
            for line in inStream:
                gId, mean, change, lfc, start, p, padj = line.strip("\n").split("\t")
                if padj != "NA" and float(padj)<0.05 and abs(float(change))>1:
                    alma.add(gId)
        AandB = len(stma & alma)
        AbutNotB = len(stma - alma)
        BbutNotA = len(alma - stma)
        norAnorB = gNr - len(alma | stma)
        
        R("""
        d = c(%i, %i, %i, %i)
        cTab = matrix(d, nrow = 2)
        fRes=fisher.test(cTab, alternative = "greater")
        write.table(data.frame(A_and_B=d[1], A_but_not_B=d[2], B_but_not_A=d[3], not_A_not_B=d[4], pvalue=fRes$p.value), "{output}", row.names=F, sep="\t")
        """ % (AandB, AbutNotB, BbutNotA, norAnorB))

rule setOverlap:
    input: sm_go="activation/ca_active_go_straw-malt.tsv", sm_kegg="activation/ca_active_kegg_straw-malt.tsv", am_go="activation/ca_active_go_alder-malt.tsv", am_kegg="activation/ca_active_kegg_alder-malt.tsv"
    output: "overlap/geneSets.tsv"
    params: cutoff=0.6
    run:
        totalGO = 0
        totalKEGG = 0
        inFiles = [input.sm_go, input.sm_kegg, input.am_go, input.am_kegg]
        sets = [set(), set(), set(), set()]
        for i, inFile in enumerate(inFiles):
            with open(inFile) as inStream:
                header = next(inStream)
                for line in inStream:
                    if i==0:
                        totalGO += 1
                    elif i==1:
                        totalKEGG += 1
                    name, inP, inS, prop, se = line.split("\t")
                    if float(prop) > params.cutoff:
                        sets[i].add(name)
        sm_go, sm_kegg, am_go, am_kegg = sets
        
        go = (len(sm_go & am_go), #A and B
              len(sm_go - am_go), #A but not B
              len(am_go - sm_go), #B but not A
              totalGO - len(am_go| sm_go)#nor A nor B
             )
        kegg = (len(sm_kegg & am_kegg), #A and B
                len(sm_kegg - am_kegg), #A but not B
                len(am_kegg - sm_kegg), #B but not A
                totalKEGG - len(am_kegg| sm_kegg)#nor A nor B
               )
        arg = go+kegg
        R("""
        go = c(%i, %i, %i, %i)
        kegg = c(%i, %i, %i, %i)
        gTab = matrix(go, nrow = 2)
        gRes=fisher.test(gTab, alternative = "greater")
        outTab = data.frame(set="GO", A_and_B=go[1], A_but_not_B=go[2], B_but_not_A=go[3], not_A_not_B=go[4], pvalue=gRes$p.value)
        kTab = matrix(kegg, nrow = 2)
        kRes=fisher.test(kTab, alternative = "greater")
        outTab = rbind(outTab, data.frame(set="KEGG", A_and_B=kegg[1], A_but_not_B=kegg[2], B_but_not_A=kegg[3], not_A_not_B=kegg[4], pvalue=kRes$p.value))
        write.table(outTab, "{output}", row.names=F, sep="\t")
        """ % arg)


rule diffAnno:
    input: exp="expr/expr_{comp}.tsv", anno="anno.pic"
    output: "expr/diffAnno_{comp}.tsv"
    params: cdb = "/home/heeger/data/dbs/cazyNames.db", kdb="/home/heeger/data/dbs/KeggPathwayIdToNameMap.db", gUrl="http://purl.obolibrary.org/obo/go.obo"
    run:
        cName = CazyNameDict(params.cdb)
        kName = CachedKeggPathwayIdToNameMap(params.kdb)
        reader = codecs.getreader("utf-8")
        gName = parseGoOboDict(reader(urlopen(params.gUrl)))
        
        anno = pickle.load(open(input.anno, "rb"))
        exp = [line.strip("\n").split("\t") for line in open(input.exp).readlines()[1:]]
        exp.sort(key=lambda x: saveFloat(x[6]))
        with open(output[0], "w") as out:
            for line in exp[:5000]:
                gene, mean, change, lfc, stat, pvalue, padj = line
                if saveFloat(padj)>0.05:
                    break
                cIds = anno["cazy"][gene]
                cNames = ["%s:%s" % (cId, cName[cId]) for cId in cIds]
                kIds = anno["kegg"][gene]
                kNames = ["%s:%s" % (kId, kName[kId]) for kId in kIds]
                gIds = anno["go"][gene]
                gNames = ["%s:%s" % (gId, gName[gId].name) for gId in gIds]
                tAnno = list(zip_longest(gNames, kNames, cNames, fillvalue=""))
                if len(tAnno) == 0:
                    out.write("%s\t%s\t%s\t\t\t\n" % (gene, change, padj))
                else:
                    out.write("%s\t%s\t%s\t%s\t%s\t%s\n" % (gene, change, padj, tAnno[0][0], tAnno[0][1], tAnno[0][2]))
                if len(tAnno) > 1:
                    for g, k, c in tAnno[1:]:
                        out.write("\t\t\t%s\t%s\t%s\n" % (g, k, c))

rule transMem:
    input: stma="expr/diff_straw-malt.tsv", alma="expr/diff_alder-malt.tsv", alst="expr/diff_alder-straw.tsv", anno="anno.pic"
    output: gl="singleGenePlots/transMem_geneList.tsv", stma="transMem/straw-malt.txt", alst="transMem/alder-straw.txt", alma="transMem/alder-malt.txt"
    run:
        anno=pickle.load(open(input.anno, "rb"))
        stma = [line.strip("\n").split("\t")[0] for line in open(input.stma).readlines()]
        alma = [line.strip("\n").split("\t")[0] for line in open(input.alma).readlines()]
        alst = [line.strip("\n").split("\t")[0] for line in open(input.alst).readlines()]
        
        tm_stma = set([g for g in stma if "GO:0055085" in anno["go"][g]])
        tm_alma = set([g for g in alma if "GO:0055085" in anno["go"][g]])
        tm_alst = set([g for g in alst if "GO:0055085" in anno["go"][g]])
        
        with open(output.gl, "w") as out:
            for g in tm_stma.union(tm_alma).union(tm_alst):
                out.write("transMem\t%s\n" % g)
        
        open(output.stma, "w").write("\n".join(tm_stma))
        open(output.alma, "w").write("\n".join(tm_alma))
        open(output.alst, "w").write("\n".join(tm_alst))

rule gshTrans:
    input: anno="anno.pic", kegg="../kegg/user_ko.tsv"
    output: gl="singleGenePlots/gshTrans_geneList.tsv"
    run:
        iprs = ["IPR016639", "IPR014440", "IPR005442", "IPR003080", "IPR003081", "IPR003082", "IPR034339"]
        anno=pickle.load(open(input.anno, "rb"))
        gt = {}
        for g in anno["interpro"]:
            for a in anno["interpro"][g]:
                if a in iprs:
                    try:
                        gt[g].append(a)
                    except KeyError:
                        gt[g] = [a]
        for line in open(input.kegg):
            arr = line.strip("\n").split("\t")
            if len(arr) > 1 and arr[1] == "K00799":
                try:
                    gt[arr[0].split("-")[0]].append(arr[1])
                except KeyError:
                    gt[arr[0].split("-")[0]] = [arr[1]]
        
        with open(output.gl, "w") as out:
            for g, annos in gt.items():
                out.write("%s\t%s\n" % (",".join(annos + [anno["signalp"][g]]), g))

rule tyr:
    input: anno="anno.pic", kegg="../kegg/user_ko.tsv"
    output: gl="singleGenePlots/tyr_geneList.tsv"
    run:
        iprs = ["IPR002227", "IPR005203", "IPR005204", "IPR000896"]
        anno=pickle.load(open(input.anno, "rb"))
        gt = {}
        for g in anno["interpro"]:
            for a in anno["interpro"][g]:
                if a in iprs:
                    try:
                        gt[g].append(a)
                    except KeyError:
                        gt[g] = [a]
        for line in open(input.kegg):
            arr = line.strip("\n").split("\t")
            if len(arr) > 1 and arr[1] == "K00505":
                try:
                    gt[arr[0].split("-")[0]].append(arr[1])
                except KeyError:
                    gt[arr[0].split("-")[0]] = [arr[1]]
        
        with open(output.gl, "w") as out:
            for g, annos in gt.items():
                out.write("%s\t%s\n" % (",".join(annos + [anno["signalp"][g]]), g))

rule cyp450:
    input: anno="anno.pic"
    output: cl="singleGenePlots/cyp450_geneList.tsv"
    run:
        iprs = ["IPR001128"]
        anno=pickle.load(open(input.anno, "rb"))
        cl = {}
        for g in anno["interpro"]:
            for a in anno["interpro"][g]:
                if a in iprs:
                    try:
                        cl[g].append(a)
                    except KeyError:
                        cl[g] = [a]
        with open(output.cl, "w") as out:
            for g, iprs in cl.items():
                out.write("%s\t%s\n" % (";".join(iprs + [anno["signalp"][g]]), g))

rule peroxi:
    input: anno="anno.pic", fullGO="ca_proteins_goFull.tsv", kegg="../kegg/user_ko.tsv"
    output: gl="singleGenePlots/peroxi_intGeneList.tsv"
    run:
        kos= {"K20205": "manganese peroxidase",
              "K15733": "dye decolorizing peroxidase",
              "K00430": "peroxidase"
            }
        gos = {"GO:0004601": "peroxidase activity"
              }
        iprs = {"IPR001621": "Fungal ligninase",
                "IPR010255": "Haem peroxidase",
                "IPR019793": "Peroxidases heam-ligand binding site"
                }
        fGo = {}
        for line in open(input.fullGO):
            go, gStr = line.strip("\n").split("\t")
            if go in gos:
                for g in gStr.split(","):
                    try:
                        fGo[g].append(go)
                    except KeyError:
                        fGo[g] = [go]
        kegg = {}
        for line in open(input.kegg):
            arr = line.strip("\n").split("\t")
            if len(arr) > 1 and arr[1] in kos:
                try:
                    kegg[arr[0].split("-")[0]].append(arr[1])
                except KeyError:
                    kegg[arr[0].split("-")[0]] = [arr[1]]
        ipr = {}
        anno=pickle.load(open(input.anno, "rb"))
        with open(output.gl, "w") as out:
            for g in anno["interpro"]:
                tAnno = []
                for a in anno["interpro"][g]:
                    if a in iprs:
                        tAnno.append("%s|%s" % (a, iprs[a]))
                if g in fGo:
                    for a in fGo[g]:
                        tAnno.append("%s|%s" % (a, gos[a]))
                if g in kegg:
                    for a in kegg[g]:
                        tAnno.append("%s|%s" % (a, kos[a]))
                if tAnno:
                    out.write("%s\t%s\n" % (";".join(tAnno + [anno["signalp"][g]]), g))

rule extractPeroxiSeq:
    input: gl="singleGenePlots/peroxi_intGeneList.tsv", prot="../funannotate/caqua/predict_results/Clavariopsis_aquatica_WDA-00-1.proteins.fa"
    output: fasta="singleGenePlots/peroxi.fasta"
    run:
        gl = set()
        for line in open(input.gl):
            g = line.strip("\n").split("\t")[1]
            gl.add(g)
        with open(output.fasta, "w") as out:
            for rec in SeqIO.parse(input.prot, "fasta"):
                if rec.id.split("-")[0] in gl:
                    out.write(rec.format("fasta"))

rule extractCazySeq:
    input: anno="anno.pic", prot="../funannotate/caqua/predict_results/Clavariopsis_aquatica_WDA-00-1.proteins.fa"
    output: fasta="singleGenePlots/cazy-{fam}.fasta"
    run:
        anno=pickle.load(open(input.anno, "rb"))
        aa = {}
        for g in anno["cazy"]:
            for a in anno["cazy"][g]:
                if a.split("_")[0] == wildcards.fam:
                    try:
                        aa[g].append(a)
                    except KeyError:
                        aa[g] = [a]
        with open(output.fasta, "w") as out:
            for rec in SeqIO.parse(input.prot, "fasta"):
                if rec.id.split("-")[0] in aa:
                    rec.id = rec.id.split("-")[0]
                    rec.description = ""
                    out.write(rec.format("fasta"))

rule cazyFam:
    input: anno="anno.pic"
    output: aa="singleGenePlots/cazy-{fam}_geneList.tsv"
    run:
        anno=pickle.load(open(input.anno, "rb"))
        aa = {}
        for g in anno["cazy"]:
            for a in anno["cazy"][g]:
                if a.split("_")[0] == wildcards.fam:
                    try:
                        aa[g].append(a)
                    except KeyError:
                        aa[g] = [a]
        with open(output.aa, "w") as out:
            for g, annos in aa.items():
                if len(annos) < 5:
                    cString = ",".join(annos + [anno["signalp"][g]])
                else:
                    cString = ",".join(annos[:5] + ["...", anno["signalp"][g]])
                out.write("%s\t%s\n" % (cString, g))

rule keggPath:
    input: anno="anno.pic", kegg="../kegg/user_ko.tsv"
    output: aa="singleGenePlots/kegg-{path}_geneList.tsv"
    run:
        anno=pickle.load(open(input.anno, "rb"))
        ko = {}
        for line in open(input.kegg):
            arr = line.strip("\n").split("\t")
            if len(arr) > 1:
                ko[arr[0].split("-")[0]] = arr[1]
        aa = {}
        for g in anno["kegg"]:
            for a in anno["kegg"][g]:
                if a == wildcards.path:
                    try:
                        aa[g].append(a)
                    except KeyError:
                        aa[g] = [a]
        with open(output.aa, "w") as out:
            for g, annos in aa.items():
                out.write("%s,%s\t%s\n" % (ko[g], anno["signalp"][g], g))

###### HELPER FUNCTIONS ###########

def saveFloat(string):
    try:
        return float(string)
    except ValueError:
        if string == "NA":
            return 1
        else:
            raise

def dezToColor(value):
    if 0 < value < 255:
        return "#FF%02x%02x" % (int(value), int(value))
    if value <= 510:
        return "#%02xFF%02x" % (510-int(value), 510-int(value))
    raise ValueError("Can not convert to color. Value out of range: %i" % value)

class CazyHTMLParser(HTMLParser):
    target = "before"
    result = []

    def handle_starttag(self, tag, attrs):
        if tag=="td":
            for name, value in attrs:
                if name == "class" and value=="tdsum" and self.target == "before":
                    self.target = "in"

    def handle_endtag(self, tag):
        if tag=="td" and self.target == "in":
            self.target = "after"

    def handle_data(self, data):
        if self.target == "in":
            self.result.append(data.strip())

class CazyNameDict(object):
    def __init__(self, dbpath):
        self.cache = PersistantDict(dbpath=dbpath)
    
    def __getitem__(self, key):
        try:
            return self.cache[key]
        except KeyError:
            value = self.htmlReq(key)
            self.cache[key] = value
            return self.cache[key]
    
    def htmlReq(self, key):
        reader = codecs.getreader("utf-8")
        a = reader(urlopen("http://www.cazy.org/%s.html" % key))
        parser=CazyHTMLParser()
        parser.feed(a.read())
        return "".join(parser.result)

########### alternative go
from goViz import main as goViz
from GoTools import createGoTree

rule explicitGoParents:
    input: "ca_proteins_go.tsv"
    output: "ca_proteins_goFull.tsv"
    run:
        reader = codecs.getreader("utf-8")
        goDict=parseGoOboDict(reader(urlopen("http://purl.obolibrary.org/obo/go.obo")))
        goTree=createGoTree(goDict)
        go2gene = {}
        for line in open(input[0]):
            goTerm, geneList = line.strip().split("\t")
            genes = geneList.split(",")
            for tGo in goTree[goTerm].ancestors():
                for g in genes:
                    try:
                        go2gene[tGo].add(g)
                    except KeyError:
                        go2gene[tGo] = set([g])
        with open(output[0], "w") as out:
            for gId, geneList in go2gene.items():
                out.write("%s\t%s\n" % (gId, ",".join(geneList)))

rule getGoFrame:
    input: "ca_proteins_goFull.tsv"
    output: "ca_proteins_goFrame.tsv"
    run:
        with open(output[0], "w") as out:
            for line in open(input[0]):
                goId, geneStr = line.strip().split("\t")
                for gene in geneStr.split(","):
                    out.write("%s\tIEA\t%s\n" % (goId, gene))

rule goEnrichment:
    input: exprFile="expr/expr_{comp}.tsv", genesFile="../funannotate/caqua/predict_results/Clavariopsis_aquatica_WDA-00-1.proteins.fa", goFile="ca_proteins_goFrame.tsv"
    output: "gsea/enrich_{comp}.tsv", rscript="gsea/enrich_{comp}.R"
    run:
        genes = [rec.id.split("-")[0] for rec in SeqIO.parse(input.genesFile, "fasta")]
        diffGenes = set([])
        with open(input.exprFile) as inFile:
            next(inFile) #skip header
            for line in inFile:
                gene, mean, change, ltfc, stat, pval, padj = line.strip("\n").split("\t")
                if padj != "NA" and float(padj)<0.05 and abs(float(change))>1:
                    diffGenes.add(gene)
        
        rcmd = 'library("GSEABase")\n' \
               'library("GOstats")\n' \
               'library("GSEABase")\n' \
               'goFrameData = read.table("%s", sep="\\t", heade=F)\n' \
               'colnames(goFrameData) = c("frame.go_id", "frame.Evidence", "frame.gene_id")\n' \
               'goFrame = GOFrame(goFrameData, organism="Clavariopsis aquatica")\n' \
               'goAllFrame=GOAllFrame(goFrame)\n' % input.goFile 
        rcmd += 'gsc = GeneSetCollection(goAllFrame, setType = GOCollection())\n'
        rcmd += 'universe=c(%s)\n' % (",".join(['"%s"' % g for g in genes]))
        rcmd += 'diffGenes=c(%s)\n' % (",".join(['"%s"' % g for g in diffGenes]))
        rcmd += 'params = GSEAGOHyperGParams(name="My Custom GSEA based annot Params",\n' \
                '         geneSetCollection=gsc,\n' \
                '         geneIds = diffGenes,\n' \
                '         universeGeneIds = universe,\n' \
                '         ontology = "MF",\n' \
                '         pvalueCutoff = 0.05,\n' \
                '         conditional = T,\n' \
                '         testDirection = "over")\n' \
                'over = hyperGTest(params)\n' \
                'write.table(summary(over), "%s", sep="\\t", quote=F, row.names=F)\n' % output[0]
        with open(output.rscript, "w") as out:
            out.write(rcmd)
        shell("Rscript %s" % output.rscript)

rule plotGoEnrichment:
    input: "gsea/enrich_{comp}.tsv"
    output: "gsea/enrich_{comp}.goGraph.svg"
    params: oboFile="/home/heeger/data/go/go.obo", sigLvl="0.05"
    run:
        with open(input[0]) as inStream, open(output[0], "w") as outStream:
            goViz(inStream, outStream, None, float(params.sigLvl))
