#load the necessary packages
library(ggplot2)
library(ggrepel)
library(pheatmap)
library(RColorBrewer)
library(dplyr)
library(pROC)

#set the work path
setwd("~/LM/biomarker_panel")
#read the expression value
dat <- read.xlsx("~/LM/data/dat_exp_all.xlsx", startRow = 2)
rownames(dat) <- dat$UniProt
dat <- dat[, -c(1,2)]
#read the sample info
sample_info <- read.xlsx("~/LM/data/patients_clinical_information.xlsx", startRow = 2)


#Discovery dat
dat <- dat[, sample_info$Patient[sample_info$cohort=="Discovery"]]

#read the DEGs
degs <- read.csv("~/LM/dif/discovery_DEGs.csv", header = T)
rownames(degs) <- degs$UniProt

#
sigGenes <- degs[c("P13232","P09237","Q8WXI7","P35968","P12544","P48023","O43557"), ]

#ready for expression of DEGs
use <- dat[rownames(sigGenes), ]
rownames(use) <- sigGenes$Assay

#预测概率
load("xgboost_model_top7.rds")
best_param <- xgboost_model$bestTune
res <- xgboost_model$pred %>%
filter(
nrounds == best_param$nrounds,
max_depth == best_param$max_depth,
eta == best_param$eta,
gamma == best_param$gamma,
colsample_bytree == best_param$colsample_bytree,
min_child_weight == best_param$min_child_weight,
subsample == best_param$subsample
)
#对每个样本预测的概率求均值
res <- res %>%
group_by(rowIndex, obs) %>%
summarise(
	LM = mean(LM, na.rm = TRUE)
) %>%
ungroup() %>% as.data.frame()
rownames(res) <- rownames(xgboost_model$trainingData)
res$sample <- rownames(res)
res$prediction <- ifelse(res$LM> 0.5, "LM", "Control")
res <- res[match(colnames(use), res$sample), ]
#ROC曲线
# 构建真实标签 & 预测概率
true_class <- ifelse(res$obs == "Control", 0, 1)
prob <- res$LM

# ROC & AUC
roc_obj <- roc(true_class, prob)
# AUC 及其 CI
ci_auc <- ci.auc(roc_obj, boot.n = 500)
auc_text <- sprintf(
  "%.3f (%.3f–%.3f)",
  auc(roc_obj),
  ci_auc[1],
  ci_auc[3]
)
ci_se <- ci.se(
  roc_obj,
  specificities = seq(0, 1, length.out = 20),
  boot.n = 500
)
pdf("train_roc.pdf", width = 4, height = 4)
plot(roc_obj,
		col = 'red',   
		legacy.axes = TRUE,  
		xlab = "1-Specificity",
		)
# CI 阴影（用基础颜色 gray）
plot(ci_se,
     type = "shape",
     col = "lightgray",
     border = NA)

# 再画一次 ROC 曲线
lines(roc_obj, col = "red", lwd = 2)

text(0.6, 0.2,  # 可以根据需要调整位置
     labels = paste0("AUC: ", auc_text),
     adj = 0, cex = 0.9)
dev.off()
#col annotation
annotation_col <- data.frame(group = ifelse(grepl("BC", colnames(use)), "BC", ifelse(grepl("LA", colnames(use)), "LA", ifelse(grepl("EC", colnames(use)), "EC", ifelse(grepl("MM", colnames(use)), "MM", "Control")))), 
                             Real = ifelse(grepl("Control", colnames(use)), "Control", "LM"),
							 prediction = res$prediction,
							 prob = res$LM)
rownames(annotation_col) <- colnames(use)

#set the annotation colors, keep the colors same as pca plot.
color_case <- brewer.pal(6, "Set1")[c(2,1)] 
anno_color <- list(group = c(Control=color_case[1], BC="#FF69B4",LA="#FFD700", MM="#8A2BE2",EC = "#20B2AA"), 
                   Real = c(Control = color_case[1], LM = color_case[2]),
				   prediction = c(Control = "#4DAF4A", LM = "#FF7F00"))
						
#adjust the order of columns
use <- use[, c(colnames(use)[grepl("BC", colnames(use))], colnames(use)[grepl("LA", colnames(use))],  colnames(use)[grepl("MM", colnames(use))], colnames(use)[grepl("EC", colnames(use))],colnames(use)[grepl("Control", colnames(use))])]

#heatmap plot:scale
#adjust the limit of colors
my_breaks <- seq(-3, 3, length.out = 101)
hm_col <- colorRampPalette(
  c("#3C5488", "#F7F7F7", "#DC0000")
)(100)
pdf(paste0("model_protein_pheatmap_discovery_scale.pdf"),width = 5,height = 2.5)
pheatmap(use,col = hm_col,cluster_row = FALSE,border_color = NA, breaks = my_breaks,cluster_col = FALSE, show_rownames= T, show_colnames = F, scale = 'row', annotation_col = annotation_col, annotation_colors = anno_color)
dev.off()
######################## Replication
##read the expression value
dat <- read.xlsx("~/LM/data/dat_exp_all.xlsx", startRow = 2)
rownames(dat) <- dat$UniProt
dat <- dat[, -c(1,2)]
#read the sample info
sample_info <- read.xlsx("~/LM/data/patients_clinical_information.xlsx", startRow = 2)

#Discovery dat
dat <- dat[, sample_info$Patient[sample_info$Cohort=="Replication"]]

#read the DEGs
degs <- read.csv("/data/home/houfei/project/Olink/combine/Correct/fourth_dif/discovery_DEGs.csv", header = T)
rownames(degs) <- degs$UniProt

#
sigGenes <- degs[c("P13232","P09237","Q8WXI7","P35968","P12544","P48023","O43557"), ]

#ready for expression of DEGs
use <- dat[rownames(sigGenes), ]
rownames(use) <- sigGenes$Assay 
#预测概率
res <- read.csv("xgboost_rep_result_top7.csv")
res$LM <- res$prediction
res$prediction <- ifelse(res$prediction > 0.5, "LM","Control")
res <- res[match(colnames(use), res$sample), ]
#ROC曲线
# 构建真实标签 & 预测概率
true_class <- ifelse(res$real_data == "Control", 0, 1)
prob <- res$LM
# ROC & AUC
roc_obj <- roc(true_class, prob)
# AUC 及其 CI
ci_auc <- ci.auc(roc_obj, boot.n = 500)
auc_text <- sprintf(
  "%.3f (%.3f–%.3f)",
  auc(roc_obj),
  ci_auc[1],
  ci_auc[3]
)
ci_se <- ci.se(
  roc_obj,
  specificities = seq(0, 1, length.out = 20),
  boot.n = 500
)
pdf("rep_roc.pdf", width = 4, height = 4)
plot(roc_obj,
		col = 'red',   
		legacy.axes = TRUE,  
		xlab = "1-Specificity",
		)
# CI 阴影（用基础颜色 gray）
plot(ci_se,
     type = "shape",
     col = "lightgray",
     border = NA)

# 再画一次 ROC 曲线
lines(roc_obj, col = "red", lwd = 2)

text(0.6, 0.2,  # 可以根据需要调整位置
     labels = paste0("AUC: ", auc_text),
     adj = 0, cex = 0.9)
dev.off()

#col annotation
annotation_col <- data.frame(group = ifelse(grepl("BC", colnames(use)), "BC", ifelse(grepl("LA", colnames(use)), "LA", ifelse(grepl("EC", colnames(use)), "EC", ifelse(grepl("MM", colnames(use)), "MM", "Control")))), 
                             Real = ifelse(grepl("Control", colnames(use)), "Control", "LM"),
							 prediction = res$prediction,
							 prob = res$LM)
rownames(annotation_col) <- colnames(use)

#set the annotation colors, keep the colors same as pca plot.
color_case <- brewer.pal(6, "Set1")[c(2,1)] 
anno_color <- list(group = c(Control=color_case[1], LA="#FFD700"), 
                   Real = c(Control = color_case[1], LM = color_case[2]),
				   prediction = c(Control = "#4DAF4A", LM = "#FF7F00"))
						

#heatmap plot:scale
#adjust the limit of colors
my_breaks <- seq(-3, 3, length.out = 101)
hm_col <- colorRampPalette(
  c("#3C5488", "#F7F7F7", "#DC0000")
)(100)
pdf(paste0("model_protein_pheatmap_rep_scale.pdf"),width = 4,height = 2.5)
pheatmap(use,col = hm_col,cluster_row = FALSE,border_color = NA, breaks = my_breaks,cluster_col = FALSE, show_rownames= T, show_colnames = F, scale = 'row', annotation_col = annotation_col, annotation_colors = anno_color)
dev.off()