library(Seurat)
library(ggpubr)
setwd("~/LM/singleCell")
load("All_used_harmony.rds")
dat <- All
dat$Celltype_Harmony <- Idents(dat)
gene <- c("FASLG","GZMA","IL7","KDR","MMP7","MUC16","TNFSF14") 
gene <- gene[gene%in%rownames(dat)]
pdf("Dotheatmap_celltype.pdf", height = 4, width=7)
DotPlot(dat,feature=gene,group.by = "Celltype_Harmony",cols = c("#AFABAB", "#B2182B"))+RotatedAxis() + scale_size_continuous(range = c(1.5, 7))
dev.off()
#umap
pdf("gene_exp_umap.pdf")
for(i in gene){
p1 <- FeaturePlot(dat, i, reduction = "umap", raster = F)
print(p1)
}
dev.off()
#