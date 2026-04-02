library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
setwd("~/LM/Figure3")
#read the expression value
dat <- read.xlsx("~/LM/data/dat_exp_all.xlsx", startRow = 2)
rownames(dat) <- dat$UniProt
dat <- dat[, -c(1,2)]
#read the sample info
sample_info <- read.xlsx("~/LM/data/patients_clinical_information.xlsx", startRow = 2)
#Discovery dat
dat <- dat[, sample_info$Patient[sample_info$Cohort=="Discovery"]]

#read the DEGs
degs <- read.csv("~/LM/dif/discovery_DEGs.csv", header= T)
rownames(degs) <- degs$UniProt
#set cutoff
logFC_cutoff <- 1
pValue_cutoff <- 0.05
#select up gene and down gene
up_gene <- subset(degs, logFC > logFC_cutoff  & adj.P.Val < pValue_cutoff)
down_gene <- subset(degs, logFC < -logFC_cutoff  & adj.P.Val < pValue_cutoff)
sigGenes <- rbind(down_gene, up_gene)

#ready for expression of DEGs
use <- dat[rownames(sigGenes), ]
rownames(use) <- sigGenes$Assay

#row annotation
annotation_row <- data.frame(direction = c(rep("DOWN", dim(down_gene)[1]), rep("UP", dim(up_gene)[1])))
rownames(annotation_row) <- rownames(use)

#col annotation
annotation_col <- data.frame( Cancer = ifelse(grepl("BC", colnames(use)), "BC", ifelse(grepl("LA", colnames(use)), "LA", ifelse(grepl("EC", colnames(use)), "EC", ifelse(grepl("MM", colnames(use)), "MM", "Control")))), case = ifelse(grepl("Control", colnames(use)), "Control", "LM"))
rownames(annotation_col) <- colnames(use)

#set the annotation colors, keep the colors same as pca plot.
anno_color <- list(Cancer = c(Control="#00A087", BC="#3C5488",LA="#E64B35", MM="#F39B7F",EC = "#4DBBD5"),direction = c(UP = "#F49C67",DOWN="#95D0A3"), condition = c(Control = "#00A087", LM = "#DC0000"))
						
#adjust the order of columns
use <- use[, c(colnames(use)[grepl("BC", colnames(use))], colnames(use)[grepl("LA", colnames(use))],  colnames(use)[grepl("EC", colnames(use))], colnames(use)[grepl("MM", colnames(use))],colnames(use)[grepl("Control", colnames(use))])]

#heatmap plot:scale
#adjust the limit of colors
my_breaks <- seq(-3, 3, length.out = 101)
hm_col <- colorRampPalette(
  c("#3C5488", "#F7F7F7", "#DC0000")
)(100)

pdf(paste0("DGE_pheatmap_discovery_scale.pdf"),width = 14,height = 6)
pheatmap(use,col = hm_col,cluster_row = TRUE, cluster_col = FALSE, border_color = NA, breaks = my_breaks,show_rownames= T, show_colnames = T, scale = 'row', annotation_col = annotation_col, annotation_colors = anno_color, annotation_row = annotation_row)
dev.off()