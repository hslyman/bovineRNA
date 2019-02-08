library(gplots)

# import data
dat <- read.table(file = "/share/lemaylab-backedup/tamu/Analysis/RNA-Seq_for_manuscript/tables_and_figures/normalized_log_RPKM_for_heatmap_short.txt", sep = "\t", header = TRUE, row.names = 1)

dat <- as.matrix(dat)

pdf("2018_07_07_bovine_RNASeq_heatmap_short.pdf")

heatmap.2(dat, col = bluered, Rowv = FALSE, dendrogram = "col", trace = "none")

dev.off()

