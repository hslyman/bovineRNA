#two factor analysis (virgin/lactating AND individual)
#load DESeq2 and biomaRt
library(DESeq2) 
library(biomaRt)

###################################################
#   set up data structures for DESeq2 analysis    #
#   and import data                               #
###################################################

#set up data structures
project_directory="~/Documents/School/Davis_GGG/milk/RNAseq/differentialExpression/DESeq2/" 
data_directory=paste(project_directory, "data/", sep = "") 
output_directory=paste(project_directory, "results/", sep = "") 

#drops unpaired sample L495
sampleFiles=grep("L495", list.files(data_directory), value = TRUE, invert = TRUE)

#these will result in a data structure containing 15 paired samples
sampleCondition=substr(sampleFiles, 0, 1)
sampleCondition=sub("V", "Virgin", sampleCondition)
sampleCondition=sub("L", "Lactating", sampleCondition)
shortName=substr(sampleFiles, 0, 4)
sampleID=substr(sampleFiles, 2, 4)
sampleTable=data.frame(sampleName=shortName,
                       fileName=sampleFiles,
                       condition=sampleCondition,
                       individual=sampleID)
#import data
ddsHTSeq= DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                     directory = data_directory,
                                     design = ~ individual + condition)

# set control and treatment using colData
# otherwise DESeq2 will assume the first condition encountered,
# in this case, Lactating, is the control condition
colData(ddsHTSeq)$condition = factor(colData(ddsHTSeq)$condition, levels=c("Virgin", "Lactating"))

###################################################
#   Run DESeq2 Analysis                           #
###################################################

dds = DESeq(ddsHTSeq)
res = results(dds)
#add Ensembl ID column to res
res$ensemblID = rownames(res)
#make object for cow gene dataset
ensembl = useMart("ensembl", dataset = "btaurus_gene_ensembl")
#create data structure mapping Ensembl IDs to gene names
genemap = getBM(attributes= c("ensembl_gene_id", "external_gene_name", "refseq_mrna", "embl"),
                filters = "ensembl_gene_id",
                values = res$ensemblID,
                mart = ensembl)
#find matching indeces
indexID = match(res$ensemblID, genemap$ensembl_gene_id)
#add column containing gene name to res
res$geneName = genemap$external_gene_name[indexID]
res$refseq = genemap$refseq_mrna[indexID]

#add baseMean values for each condition
baseMeanConditions = sapply(levels(dds$condition), function(lvl) rowMeans(counts(dds,normalized=TRUE) [,dds$condition == lvl]))
bmc = as.data.frame(baseMeanConditions)
res$baseMeanVirgin = bmc$Virgin
res$baseMeanLactating = bmc$Lactating
#sort results on adjusted p-value
res = res[order(res$padj), ]


###################################################
#   Create plots from DESeq2 Analysis             #
###################################################

#create MA plot (log2 fold change vs mean normalized count, with 0-centered normal prior)
resMA = plotMA(dds, ylim=c(-4,4), main="MA plot (formula = ~ individual + condition)")

#plot counts of alpha casein subunit one in all samples
resCSN1S1 = plotCounts(dds, gene = "ENSBTAG00000007695", intgroup = "condition", main = "Counts of CSN1S1" )

#make heat map
resRLD = rlogTransformation(dds)
resDistsRL = dist(t(assay(resRLD)))
mat = as.matrix(resDistsRL)
rownames(mat) = colnames(mat) = shortName
library("gplots")
library("RColorBrewer")
colors=colorRampPalette(rev(brewer.pal(9, "Blues")))(255)
resHeat = heatmap.2(mat, trace="none", col=colors)

#make PCA
resPCA = plotPCA(resRLD, intgroup="condition")
print(resPCA)

###################################################
#   Write output (significant loci only)          #
###################################################

#isolate entries with adjusted p-values less than 0.05
resSig = res[which(res$padj < 0.05),c("geneName",
                                      "refseq",
                                      "genbank",
                                      "baseMeanVirgin",
                                      "baseMeanLactating",
                                      "log2FoldChange",
                                      "padj")]

#write output to tab-delimited file
write.table(resSig, file = "DESeq2_results.txt", quote = FALSE, sep = "\t")

