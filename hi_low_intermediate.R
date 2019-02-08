#compare lactating samples to each other
#three groups, based on their PCA clusters
#load DESeq2 and biomaRt
library(DESeq2)
library(biomaRt)

###################################################
#   set up data structures for DESeq2 analysis    #
#   and import data                               #
###################################################

#set up data structures
directory="~/Documents/School/Davis_GGG/milk/RNASeq/differentialExpression/DESeq2/data/"
sampleFiles=grep("L", list.files(directory), value = TRUE)
sampleCondition=c("high", "low", "high", "intermediate", "high", "high", "intermediate", "high", "high", "intermediate", "low", "intermediate", "high", "high", "high", "intermediate")
shortName=substr(sampleFiles, 0, 4)
sampleTable=data.frame(sampleName=shortName,
                       fileName=sampleFiles,
                       condition=sampleCondition)
#import data
ddsHTSeq= DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                     directory = directory,
                                     design = ~condition)

#set control and treatment 
colData(ddsHTSeq)$condition = factor(colData(ddsHTSeq)$condition, levels=c("high", "intermediate", "low"))

###################################################
#   Run DESeq2 Analysis                           #
###################################################

dds = DESeq(ddsHTSeq)
res_hi_low = results(dds, contrast = c("condition", "high", "low"))
res_hi_int = results(dds, contrast = c("condition", "high", "intermediate"))
res_int_low = results(dds, contrast = c("condition", "intermediate", "low"))
#add Ensembl ID column to res
res_hi_low$ensemblID = rownames(res_hi_low)
res_hi_int$ensemblID = rownames(res_hi_int)
res_int_low$ensemblID = rownames(res_int_low)
#make object for cow gene dataset
ensembl = useMart("ensembl", dataset = "btaurus_gene_ensembl")
#create data structure mapping Ensembl IDs to gene names
genemap = getBM(attributes= c("ensembl_gene_id", "external_gene_name", "refseq_mrna", "embl"),
                filters = "ensembl_gene_id",
                values = res_hi_low$ensemblID,
                mart = ensembl)
#find matching indeces
indexID = match(res_hi_low$ensemblID, genemap$ensembl_gene_id)
#add column containing gene name to res
res_hi_low$geneName = genemap$external_gene_name[indexID]
res_hi_low$refseq = genemap$refseq_mrna[indexID]
res_hi_int$geneName = genemap$external_gene_name[indexID]
res_hi_int$refseq = genemap$refseq_mrna[indexID]
res_int_low$geneName = genemap$external_gene_name[indexID]
res_int_low$refseq = genemap$refseq_mrna[indexID]
#add baseMean values for each condition
baseMeanConditions = sapply(levels(dds$condition), function(lvl) rowMeans(counts(dds,normalized=TRUE) [,dds$condition == lvl]))
bmc = as.data.frame(baseMeanConditions)
res_hi_low$baseMeanHigh = bmc$high
res_hi_low$baseMeanIntermediate = bmc$intermediate
res_hi_low$baseMeanLow = bmc$low
res_hi_int$baseMeanHigh = bmc$high
res_hi_int$baseMeanIntermediate = bmc$intermediate
res_hi_int$baseMeanLow = bmc$low
res_int_low$baseMeanHigh = bmc$high
res_int_low$baseMeanIntermediate = bmc$intermediate
res_int_low$baseMeanLow = bmc$low
#sort results on adjusted p-value
res_hi_low = res_hi_low[order(res_hi_low$padj), ]
res_hi_int = res_hi_int[order(res_hi_int$padj), ]
res_int_low = res_int_low[order(res_int_low$padj), ]

###################################################
#   Write output (significant loci only)          #
###################################################

#isolate entries with adjusted p-values less than 0.05
resSig_hi_low = res_hi_low[which(res_hi_low$padj < 0.05),c("geneName",
                                      "refseq",
                                      "baseMeanHigh",
                                      "baseMeanIntermediate",
                                      "baseMeanLow",
                                      "log2FoldChange",
                                      "padj")]
#isolate entries with adjusted p-values less than 0.05
resSig_hi_int = res_hi_int[which(res_hi_int$padj < 0.05),c("geneName",
                                                           "refseq",
                                                           "baseMeanHigh",
                                                           "baseMeanIntermediate",
                                                           "baseMeanLow",
                                                           "log2FoldChange",
                                                           "padj")]
#isolate entries with adjusted p-values less than 0.05
resSig_int_low = res_int_low[which(res_int_low$padj < 0.05),c("geneName",
                                                           "refseq",
                                                           "baseMeanHigh",
                                                           "baseMeanIntermediate",
                                                           "baseMeanLow",
                                                           "log2FoldChange",
                                                           "padj")]

#write output to tab-delimited file
write.table(resSig_hi_low, file = "~/Desktop/DESeq2_results_hi_low.txt", quote = FALSE, sep = "\t")
write.table(resSig_hi_int, file = "~/Desktop/DESeq2_results_hi_int.txt", quote = FALSE, sep = "\t")
write.table(resSig_int_low, file = "~/Desktop/DESeq2_results_int_low.txt", quote = FALSE, sep = "\t")
