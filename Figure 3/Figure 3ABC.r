setwd("~/LM/Figure3")
#
#read the DEGs
degs <- read.csv("~/LM/dif/discovery_DEGs.csv", header = T)
# degs <- read.csv("~/LM/dif/discovery_breast_DEGs.csv", header = T)
# degs <- read.csv("~/LM/dif/discovery_lung_DEGs.csv", header = T)
# degs <- read.csv("~/LM/dif/replication_DEGs.csv", header = T)

degs$direction <- factor(degs$direction, levels = c("NOT","Increased","Decreased"))
#calculating the number of up or down genes
up_num <- length(degs$direction[degs$direction=="Increased"])
down_num <- length(degs$direction[degs$direction=="Decreased"])

#signify the significant
degs$label <- ifelse(abs(degs$logFC) > 1.5 & -log10(degs$P.Value) > 3, degs$Assay, "")

#set the colors of directions
color_manual <- c("#7F8084","#D97C3F","#5FAF87") #Not，Increased，Decreased
directions <- names(table(degs$direction))
if("Increased" %in% directions & length(directions) == 2){manual_color_value <- color_manual[c(1,2)]}
if("Decreased" %in% directions & length(directions) == 2){manual_color_value <- color_manual[c(1,3)]}
if(length(directions) == 3){manual_color_value <-  color_manual[c(1, 2, 3)]}
if("NOT" %in% directions & length(directions) == 1){manual_color_value <-  color_manual[1]}

#set cutoff
logFC_cutoff <- 1
pValue_cutoff <- 0.05

#volcano plot
p1 <- ggplot(degs,aes(logFC,-1*log10(adj.P.Val))) +
		geom_point(aes(color = direction)) +   
		labs(title=paste0("Increased: ", up_num, "; Decreased: ", down_num),x="log2(FC)",y="-log10(padj)") +
		scale_color_manual(values = manual_color_value) +
		geom_hline(yintercept=-log10(pValue_cutoff),linetype=2)+
		geom_vline(xintercept=c(-logFC_cutoff,logFC_cutoff),linetype=2)+ 
		theme(legend.background=element_rect(fill="transparent"),axis.line = element_line(color = "black"),panel.background=element_rect(fill="transparent"),plot.title = element_text(hjust = 0.5))+
		geom_text_repel(aes(x = logFC, y = -1*log10(adj.P.Val), label=label), max.overlaps = 10000,size=3,box.padding=unit(0.5,'lines'),point.padding=unit(0.1, 'lines'),segment.color='black',show.legend=FALSE)
pdf("discovery_volcanoplot.pdf", width = 4.5, height = 4)
print(p1)
dev.off()