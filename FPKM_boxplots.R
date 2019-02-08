library(ggplot2)

# load in epithelial marker FPKM data from file
epithelial = read.table("~/Documents/School/Davis_GGG/milk/RNASeq/IMGC_2016_Abstract/fpkm_tracking/sig_med_hi_only_normalized.txt", header = TRUE)

epithelial$cluster = factor(epithelial$cluster, levels = c("low", "intermediate", "high"))

# create boxplot from epithelial markers
pdf(file = "~/Documents/School/Davis_GGG/milk/RNASeq/IMGC_2016_Abstract/epithelial_boxplot.pdf")
ggplot(epithelial, aes(x = cluster, y = fpkm, fill = cluster)) + geom_boxplot() + facet_wrap(~ gene, scales = "free") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("FPKM / internal control") + xlab("")
dev.off()

# load in casein FPKM data from file
caseins = read.table("~/Documents/School/Davis_GGG/milk/RNASeq/IMGC_2016_Abstract/fpkm_tracking/caseins_fpkm_normalized_to_internal_control.txt", header = TRUE)

caseins$cluster = factor(caseins$cluster, levels = c("low", "intermediate", "high"))

# create boxplot from caseins
pdf(file = "~/Documents/School/Davis_GGG/milk/RNASeq/IMGC_2016_Abstract/caseins_boxplot.pdf")
ggplot(caseins, aes(x = cluster, y = fpkm, fill = cluster)) + geom_boxplot() + facet_wrap(~ gene, scales = "free") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("FPKM / internal control") + xlab("")
dev.off()

# load in yield data from file
yield = read.table("~/Documents/School/Davis_GGG/milk/RNASeq/IMGC_2016_Abstract/three_clusters_lactation/cluster_milk_yield2.txt", header = TRUE)

yield$cluster = factor(yield$cluster, levels = c("low", "intermediate", "high"))

ggplot(yield, aes(x = cluster, y = volume, fill = cluster)) + geom_boxplot() + facet_wrap(~type, scales = "free") + theme(axis.text.x = element_blank(), axis.ticks.x = element_blank()) + ylab("yield") + xlab("")
