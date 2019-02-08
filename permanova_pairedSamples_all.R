library(DESeq2)
library(reshape2)
library(vegan)

# set up data
data_directory = "~/Documents/School/Davis_GGG/milk/RNAseq/differentialExpression/DESeq2/data"

sampleFiles=grep("L495", list.files(data_directory), value = TRUE, invert = TRUE)
sampleCondition=substr(sampleFiles, 0, 1)
sampleCondition=sub("V", "Virgin", sampleCondition)
sampleCondition=sub("L", "Lactating", sampleCondition)
shortName=substr(sampleFiles, 0, 4)
sampleID=substr(sampleFiles, 2, 4)
sampleTable=data.frame(sampleName=shortName,
                       fileName=sampleFiles,
                       condition=sampleCondition,
                       animal=sampleID)
# import data
ddsHTSeq= DESeqDataSetFromHTSeqCount(sampleTable = sampleTable,
                                     directory = data_directory,
                                     design = ~ animal + condition)

# set control and treatment for experiment
colData(ddsHTSeq)$condition = factor(colData(ddsHTSeq)$condition, levels=c("Virgin", "Lactating"))

dds = DESeq(ddsHTSeq)

# regularized log transform to stabilize variance across mean
# this avoids skewing distances toward highly-expressed transcripts
rld = rlog(dds, blind=FALSE)

# get euclidian distances between log-transformed samples
sampleDists = dist(t(assay(rld)))
sampleDistMatrix = as.matrix(sampleDists)
#rownames(sampleDistMatrix) = paste( rld$condition, rld$animal, sep="" )
#colnames(sampleDistMatrix) = NULL

adonis(sampleDistMatrix ~ condition+animal+condition:animal, sampleTable, method = "euclidian")

adonis(sampleDistMatrix ~ condition, sampleTable, method = "euclidian")

