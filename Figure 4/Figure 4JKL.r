library(pROC) #coords
library(caret) #confusionMatrix
library(ggplot2)
library(reshape2)  # 用于转换数据格式
library(dplyr)
library(tibble)
library(tidyr)
#
setwd("~/LM/biomarker_panel")
protein <- c("P35968","P13232","P12544","P09237","P48023","Q8WXI7","O43557")
dat <- read.xlsx("~/LM/data/dat_exp_all.xlsx", startRow = 2)
rownames(dat) <- dat$UniProt
dat <- dat[, -c(1,2)]
sample_info <- read.xlsx("~/LM/data/patients_clinical_information.xlsx", startRow = 2)
dat_dis <- dat[protein, sample_info$Patient[sample_info$Cohort=="Discovery"]]
dat_rep <- dat[protein, sample_info$Patient[sample_info$Cohort=="Replication"]]
##################################single gene of LM-IO7
#####################dis set
#our model
load("xgboost_model_top7.rds")
best_params <- xgboost_model$bestTune
pred <- xgboost_model$pred %>%
	filter(nrounds == best_params$nrounds,
	       max_depth == best_params$max_depth,
		   eta == best_params$eta,
		   gamma == best_params$gamma,
		   colsample_bytree == best_params$colsample_bytree,
		   min_child_weight == best_params$min_child_weight,
		   subsample == best_params$subsample)
res <- pred %>%
    group_by(rowIndex, obs) %>%
    summarise(
      LM = mean(LM, na.rm = TRUE),
      .groups = "drop"
    ) %>% as.data.frame()

true_class <- ifelse(res$obs=="Control", 0, 1)
prob <- res$LM
roc_myModel <- roc(true_class, prob)  # 计算 ROC 曲线
#VEGFR2
exp_value <- dat_dis["P35968", ] %>% as.numeric()
true_class <- ifelse(grepl("Control",colnames(dat_dis)), 0, 1)
roc_VEGFR2 <- roc(true_class, exp_value, direction = ">")  # 计算 ROC 曲线
#IL7
exp_value <- dat_dis["P13232", ] %>% as.numeric()
roc_IL7 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#GZMA
exp_value <- dat_dis["P12544", ] %>% as.numeric()
roc_GZMA <- roc(true_class, exp_value, direction = ">")  # 计算 ROC 曲线
#MMP7
exp_value <- dat_dis["P09237", ] %>% as.numeric()
roc_MMP7 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#FASLG
exp_value <- dat_dis["P48023", ] %>% as.numeric()
roc_FASLG <- roc(true_class, exp_value,direction = ">")  # 计算 ROC 曲线
#MUC16
exp_value <- dat_dis["Q8WXI7", ] %>% as.numeric()
roc_MUC16 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#TNFSF14
exp_value <- dat_dis["O43557", ] %>% as.numeric()
roc_TNFSF14 <- roc(true_class, exp_value, direction = ">")  # 计算 ROC 曲线
# 计算 AUC
auc_myModel <- auc(roc_myModel)
auc_VEGFR2 <- auc(roc_VEGFR2)
auc_IL7   <- auc(roc_IL7)
auc_GZMA  <- auc(roc_GZMA)
auc_MMP7    <- auc(roc_MMP7)
auc_FASLG   <- auc(roc_FASLG)
auc_MUC16   <- auc(roc_MUC16)
auc_TNFSF14  <- auc(roc_TNFSF14)

# 颜色（可按需调整）
cols <- c(
  myModel = "#1f77b4",
  VEGFR2   = "#ff7f0e",
  IL7  = "#2ca02c",
  GZMA    = "#d62728",
  MMP7   = "#9467bd",
  FASLG = "#8c564b",
  MUC16 = "#e377c2",
  TNFSF14 = "#7f7f7f"
)

# 绘图
pdf("roc_dis_single_gene.pdf", width = 4, height = 4)
op <- par(mar = c(4.2, 4.2, 3.2, 1.2))
plot(roc_myModel, col = cols["myModel"], lwd = 1.4, legacy.axes = TRUE,
     main = "dis Set",
     xlab = "1-specificity", ylab = "sensitivity")
plot(roc_VEGFR2, col = cols["VEGFR2"], lwd = 1.4, add = TRUE)
plot(roc_IL7,   col = cols["IL7"],   lwd = 1.4, add = TRUE)
plot(roc_GZMA,  col = cols["GZMA"],  lwd = 1.4, add = TRUE)
plot(roc_MMP7,  col = cols["MMP7"],  lwd = 1.4, add = TRUE)
plot(roc_FASLG,  col = cols["FASLG"],  lwd = 1.4, add = TRUE)
plot(roc_MUC16,  col = cols["MUC16"],  lwd = 1.4, add = TRUE)
plot(roc_TNFSF14,  col = cols["TNFSF14"],  lwd = 1.4, add = TRUE)
abline(0, 1, lty = 3, col = "gray50")

legend("bottomright",
       legend = c(
         sprintf("Our model (AUC = %.3f)", auc_myModel),
         sprintf("VEGFR2 (AUC = %.3f)",     auc_VEGFR2),
         sprintf("IL7 (AUC = %.3f)",    auc_IL7),
         sprintf("GZMA (AUC = %.3f)",      auc_GZMA),
         sprintf("MMP7 (AUC = %.3f)",     auc_MMP7),
		 sprintf("FASLG (AUC = %.3f)",     auc_FASLG),
		 sprintf("MUC16 (AUC = %.3f)",     auc_MUC16),
		 sprintf("TNFSF14 (AUC = %.3f)",     auc_TNFSF14)
       ),
       col = cols, lwd = 2.8, cex = 0.6, bty = "n")
par(op)
dev.off()
## ========== 1) 其他基因基于最优阈值(closest.topleft) + 指标汇总 ==========
#
get_metrics <- function(roc_obj, name, prob, true_class){
if(name == "myModel"){
	best_threshold = 0.5
	ACC = 0.972
	F1 = 0.982
	Sensitivity = 0.988
	Specificity = 0.926
	Precision = 0.976
	Recall = 0.988
	bACC = 0.957
} else {
co <- coords(roc_obj, x = "best", best.method = "closest.topleft",
               ret = c("threshold","sensitivity","specificity","ppv","npv"),
               transpose = FALSE)
best_threshold <- co["threshold"] %>% as.numeric()
if(name%in%c("VEGFR2","GZMA","FASLG","TNFSF14")){
	pred_class <- ifelse(prob < best_threshold, 1, 0)
} else {
	pred_class <- ifelse(prob > best_threshold, 1, 0)
}
pred_class <- factor(pred_class, levels = c(0, 1))
true_class <- factor(true_class, levels = c(0, 1))
cm <- confusionMatrix(pred_class, true_class, positive = "1")
ACC  <- as.numeric(signif(cm$overall["Accuracy"], 3))
F1   <- as.numeric(signif(cm$byClass["F1"], 3))
Sensitivity <- as.numeric(signif(cm$byClass["Sensitivity"], 3))
Specificity <- as.numeric(signif(cm$byClass["Specificity"], 3))
Precision <- as.numeric(signif(cm$byClass["Precision"], 3))
Recall <- as.numeric(signif(cm$byClass["Recall"], 3))
bACC <- as.numeric(signif(cm$byClass["Balanced Accuracy"], 3))
}
data.frame(
    Curve = name,
    AUC = as.numeric(auc(roc_obj)),
    threshold = best_threshold,
    sensitivity = Sensitivity,
    specificity = Specificity,
	ACC = ACC,
	F1 = F1,
	Precision = Precision,
	Recall = Recall,
	bACC = bACC,
    stringsAsFactors = FALSE
  )
}
metrics_tbl <- do.call(rbind, list(
  get_metrics(roc_myModel, "myModel", NA, NA),
  get_metrics(roc_VEGFR2,   "VEGFR2", dat_dis["P35968", ] %>% as.numeric(), true_class),
  get_metrics(roc_IL7,  "IL7", dat_dis["P13232", ] %>% as.numeric(), true_class),
  get_metrics(roc_GZMA,    "GZMA", dat_dis["P12544", ] %>% as.numeric(), true_class),
  get_metrics(roc_MMP7,   "MMP7", dat_dis["P09237", ] %>% as.numeric(), true_class),
  get_metrics(roc_FASLG,   "FASLG", dat_dis["P48023", ] %>% as.numeric(), true_class),
  get_metrics(roc_MUC16,   "MUC16", dat_dis["Q8WXI7", ] %>% as.numeric(), true_class),
  get_metrics(roc_TNFSF14,   "TNFSF14", dat_dis["O43557", ] %>% as.numeric(), true_class)
))

## -------- 使用 signif() 控制有效数字 --------
metrics_tbl$AUC         <- signif(metrics_tbl$AUC, 3)
metrics_tbl$threshold   <- signif(metrics_tbl$threshold, 3)
metrics_tbl$sensitivity <- signif(metrics_tbl$sensitivity, 3)
metrics_tbl$specificity <- signif(metrics_tbl$specificity, 3)
metrics_tbl$ACC <- signif(metrics_tbl$ACC, 3)
metrics_tbl$F1 <- signif(metrics_tbl$F1, 3)
metrics_tbl$Precision <- signif(metrics_tbl$Precision, 3)
metrics_tbl$Recall <- signif(metrics_tbl$Recall, 3)
metrics_tbl$bACC <- signif(metrics_tbl$bACC, 3)


print(metrics_tbl)

## 保存结果
write.csv(metrics_tbl, "evaluation_dis_single_gene.csv", row.names = FALSE)
##可视化
# 假设 metrics_tbl 中已有列：Curve, ACC, sensitivity, specificity
metrics_tbl_long <- metrics_tbl %>%
  select(Curve, ACC, sensitivity, specificity) %>%
  pivot_longer(cols = c(ACC, sensitivity, specificity),
               names_to = "Metric", values_to = "Value")

# 调整顺序（让ACC显示在最前）
metrics_tbl_long$Metric <- factor(metrics_tbl_long$Metric,
                                  levels = c("ACC", "sensitivity", "specificity"))

# 颜色可自定义
metric_colors <- c(
  "ACC" = "#1f77b4",
  "sensitivity" = "#ff7f0e",
  "specificity" = "#2ca02c"
)
# 绘图
metrics_tbl_long$Curve <- factor(metrics_tbl_long$Curve, levels = c("myModel", "VEGFR2", "IL7","GZMA","MMP7","FASLG","MUC16","TNFSF14"))
p <- ggplot(metrics_tbl_long,
            aes(x = Curve, y = Value, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = signif(Value, 3)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.3) +
  scale_fill_manual(values = metric_colors) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(title = "dis set",
       x = NULL, y = "Value") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top", axis.text.x = element_text(angle = 15, hjust = 1))
pdf("ACC_Spe_Sen_dis_single_gene.pdf")
print(p)
dev.off()
#####################Replication set
#our model
res <- read.csv("xgboost_rep_result_top7.csv")
true_class <- ifelse(res$real_data=="Control", 0, 1)
prob <- res$prediction
roc_myModel <- roc(true_class, prob)  # 计算 ROC 曲线
#VEGFR2
exp_value <- dat_rep["P35968", ] %>% as.numeric()
true_class <- ifelse(grepl("Control",colnames(dat_rep)), 0, 1)
roc_VEGFR2 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#IL7
exp_value <- dat_rep["P13232", ] %>% as.numeric()
roc_IL7 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#GZMA
exp_value <- dat_rep["P12544", ] %>% as.numeric()
roc_GZMA <- roc(true_class, exp_value)  # 计算 ROC 曲线
#MMP7
exp_value <- dat_rep["P09237", ] %>% as.numeric()
roc_MMP7 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#FASLG
exp_value <- dat_rep["P48023", ] %>% as.numeric()
roc_FASLG <- roc(true_class, exp_value)  # 计算 ROC 曲线
#MUC16
exp_value <- dat_rep["Q8WXI7", ] %>% as.numeric()
roc_MUC16 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#TNFSF14
exp_value <- dat_rep["O43557", ] %>% as.numeric()
roc_TNFSF14 <- roc(true_class, exp_value)  # 计算 ROC 曲线
# 计算 AUC
auc_myModel <- auc(roc_myModel)
auc_VEGFR2 <- auc(roc_VEGFR2)
auc_IL7   <- auc(roc_IL7)
auc_GZMA  <- auc(roc_GZMA)
auc_MMP7    <- auc(roc_MMP7)
auc_FASLG   <- auc(roc_FASLG)
auc_MUC16   <- auc(roc_MUC16)
auc_TNFSF14  <- auc(roc_TNFSF14)

# 颜色（可按需调整）
cols <- c(
  myModel = "#1f77b4",
  VEGFR2   = "#ff7f0e",
  IL7  = "#2ca02c",
  GZMA    = "#d62728",
  MMP7   = "#9467bd",
  FASLG = "#8c564b",
  MUC16 = "#e377c2",
  TNFSF14 = "#7f7f7f"
)

# 绘图
pdf("roc_rep_single_gene.pdf", width = 4, height = 4)
op <- par(mar = c(4.2, 4.2, 3.2, 1.2))
plot(roc_myModel, col = cols["myModel"], lwd = 1.4, legacy.axes = TRUE,
     main = "rep Set",
     xlab = "1-specificity", ylab = "sensitivity")
plot(roc_VEGFR2, col = cols["VEGFR2"], lwd = 1.4, add = TRUE)
plot(roc_IL7,   col = cols["IL7"],   lwd = 1.4, add = TRUE)
plot(roc_GZMA,  col = cols["GZMA"],  lwd = 1.4, add = TRUE)
plot(roc_MMP7,  col = cols["MMP7"],  lwd = 1.4, add = TRUE)
plot(roc_FASLG,  col = cols["FASLG"],  lwd = 1.4, add = TRUE)
plot(roc_MUC16,  col = cols["MUC16"],  lwd = 1.4, add = TRUE)
plot(roc_TNFSF14,  col = cols["TNFSF14"],  lwd = 1.4, add = TRUE)
abline(0, 1, lty = 3, col = "gray50")

legend("bottomright",
       legend = c(
         sprintf("Our model (AUC = %.3f)", auc_myModel),
         sprintf("VEGFR2 (AUC = %.3f)",     auc_VEGFR2),
         sprintf("IL7 (AUC = %.3f)",    auc_IL7),
         sprintf("GZMA (AUC = %.3f)",      auc_GZMA),
         sprintf("MMP7 (AUC = %.3f)",     auc_MMP7),
		 sprintf("FASLG (AUC = %.3f)",     auc_FASLG),
		 sprintf("MUC16 (AUC = %.3f)",     auc_MUC16),
		 sprintf("TNFSF14 (AUC = %.3f)",     auc_TNFSF14)
       ),
       col = cols, lwd = 2.8, cex = 0.6, bty = "n")
par(op)
dev.off()
## 
#
get_metrics <- function(roc_obj, name, prob, true_class){
if(name == "myModel"){
	best_threshold = 0.5
	ACC = 0.839
	F1 = 0.783
	Sensitivity = 0.818
	Specificity = 0.850
	Precision = 0.750
	Recall = 0.818
	bACC = 0.834
} else {
best_threshold <- metrics_tbl[metrics_tbl$Curve==name, "threshold"] %>% as.numeric()
if(name%in%c("VEGFR2","GZMA","FASLG","TNFSF14")){
	pred_class <- ifelse(prob < best_threshold, 1, 0)
} else {
	pred_class <- ifelse(prob > best_threshold, 1, 0)
}
pred_class <- factor(pred_class, levels = c(0, 1))
true_class <- factor(true_class, levels = c(0, 1))
cm <- confusionMatrix(pred_class, true_class, positive = "1")
ACC  <- as.numeric(signif(cm$overall["Accuracy"], 3))
F1   <- as.numeric(signif(cm$byClass["F1"], 3))
Sensitivity <- as.numeric(signif(cm$byClass["Sensitivity"], 3))
Specificity <- as.numeric(signif(cm$byClass["Specificity"], 3))
Precision <- as.numeric(signif(cm$byClass["Precision"], 3))
Recall <- as.numeric(signif(cm$byClass["Recall"], 3))
bACC <- as.numeric(signif(cm$byClass["Balanced Accuracy"], 3))
}
data.frame(
    Curve = name,
    AUC = as.numeric(auc(roc_obj)),
    threshold = best_threshold,
    sensitivity = Sensitivity,
    specificity = Specificity,
	ACC = ACC,
	F1 = F1,
	Precision = Precision,
	Recall = Recall,
	bACC = bACC,
    stringsAsFactors = FALSE
  )
}
metrics_tbl <- do.call(rbind, list(
  get_metrics(roc_myModel, "myModel", NA, NA),
  get_metrics(roc_VEGFR2,   "VEGFR2", dat_rep["P35968", ] %>% as.numeric(), true_class),
  get_metrics(roc_IL7,  "IL7", dat_rep["P13232", ] %>% as.numeric(), true_class),
  get_metrics(roc_GZMA,    "GZMA", dat_rep["P12544", ] %>% as.numeric(), true_class),
  get_metrics(roc_MMP7,   "MMP7", dat_rep["P09237", ] %>% as.numeric(), true_class),
  get_metrics(roc_FASLG,   "FASLG", dat_rep["P48023", ] %>% as.numeric(), true_class),
  get_metrics(roc_MUC16,   "MUC16", dat_rep["Q8WXI7", ] %>% as.numeric(), true_class),
  get_metrics(roc_TNFSF14,   "TNFSF14", dat_rep["O43557", ] %>% as.numeric(), true_class)
))

## -------- 使用 signif() 控制有效数字 --------
metrics_tbl$AUC         <- signif(metrics_tbl$AUC, 3)
metrics_tbl$threshold   <- signif(metrics_tbl$threshold, 3)
metrics_tbl$sensitivity <- signif(metrics_tbl$sensitivity, 3)
metrics_tbl$specificity <- signif(metrics_tbl$specificity, 3)
metrics_tbl$ACC <- signif(metrics_tbl$ACC, 3)
metrics_tbl$F1 <- signif(metrics_tbl$F1, 3)
metrics_tbl$Precision <- signif(metrics_tbl$Precision, 3)
metrics_tbl$Recall <- signif(metrics_tbl$Recall, 3)
metrics_tbl$bACC <- signif(metrics_tbl$bACC, 3)


print(metrics_tbl)

## 保存结果
write.csv(metrics_tbl, "evaluation_rep_single_gene.csv", row.names = FALSE)
##可视化
# 假设 metrics_tbl 中已有列：Curve, ACC, sensitivity, specificity
metrics_tbl_long <- metrics_tbl %>%
  select(Curve, ACC, sensitivity, specificity) %>%
  pivot_longer(cols = c(ACC, sensitivity, specificity),
               names_to = "Metric", values_to = "Value")

# 调整顺序（让ACC显示在最前）
metrics_tbl_long$Metric <- factor(metrics_tbl_long$Metric,
                                  levels = c("ACC", "sensitivity", "specificity"))

# 颜色可自定义
metric_colors <- c(
  "ACC" = "#1f77b4",
  "sensitivity" = "#ff7f0e",
  "specificity" = "#2ca02c"
)
# 绘图
metrics_tbl_long$Curve <- factor(metrics_tbl_long$Curve, levels = c("myModel", "VEGFR2", "IL7","GZMA","MMP7","FASLG","MUC16","TNFSF14"))
p <- ggplot(metrics_tbl_long,
            aes(x = Curve, y = Value, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = signif(Value, 3)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.3) +
  scale_fill_manual(values = metric_colors) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(title = "rep set",
       x = NULL, y = "Value") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 15, hjust = 1))
pdf("ACC_Spe_Sen_rep_single_gene.pdf")
print(p)
dev.off()

#####################################established biomarkers
#基于已知的单个基因标志物的表达，直接预测LM
library(pROC) #coords
library(caret) #confusionMatrix
library(ggplot2)
library(reshape2)  # 用于转换数据格式
library(dplyr)
library(tibble)
library(tidyr)
protein <- c("P15692","P02778","Q14116","P01137")

#####################dis set
#our model
load("xgboost_model_top7.rds")
best_params <- xgboost_model$bestTune
pred <- xgboost_model$pred %>%
	filter(nrounds == best_params$nrounds,
	       max_depth == best_params$max_depth,
		   eta == best_params$eta,
		   gamma == best_params$gamma,
		   colsample_bytree == best_params$colsample_bytree,
		   min_child_weight == best_params$min_child_weight,
		   subsample == best_params$subsample)
res <- pred %>%
    group_by(rowIndex, obs) %>%
    summarise(
      LM = mean(LM, na.rm = TRUE),
      .groups = "drop"
    ) %>% as.data.frame()

true_class <- ifelse(res$obs=="Control", 0, 1)
prob <- res$LM
roc_myModel <- roc(true_class, prob)  # 计算 ROC 曲线
#VEGFA
exp_value <- dat_dis["P15692", ] %>% as.numeric()
true_class <- ifelse(grepl("Control",colnames(dat_dis)), 0, 1)
roc_VEGFA <- roc(true_class, exp_value)  # 计算 ROC 曲线
#CXCL10
exp_value <- dat_dis["P02778", ] %>% as.numeric()
true_class <- ifelse(grepl("Control",colnames(dat_dis)), 0, 1)
roc_CXCL10 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#IL18，这个蛋白表达越小，预测为case，需要根据文献核实是否对应。
exp_value <- dat_dis["Q14116", ] %>% as.numeric()
true_class <- ifelse(grepl("Control",colnames(dat_dis)), 0, 1)
roc_IL18 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#TGFB1
exp_value <- dat_dis["P01137", ] %>% as.numeric()
true_class <- ifelse(grepl("Control",colnames(dat_dis)), 0, 1)
roc_TGFB1 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#

# 计算 AUC
auc_myModel <- auc(roc_myModel)
auc_VEGFA   <- auc(roc_VEGFA)
auc_CXCL10  <- auc(roc_CXCL10)
auc_IL18    <- auc(roc_IL18)
auc_TGFB1   <- auc(roc_TGFB1)

# 颜色（可按需调整）
cols <- c(
  myModel = "#1f77b4",
  VEGFA   = "#ff7f0e",
  CXCL10  = "#2ca02c",
  IL18    = "#d62728",
  TGFB1   = "#9467bd"
)

# 绘图
pdf("roc_dis.pdf", width = 4, height = 4)
op <- par(mar = c(4.2, 4.2, 3.2, 1.2))
plot(roc_myModel, col = cols["myModel"], lwd = 1.4, legacy.axes = TRUE,
     main = "dis Set",
     xlab = "1-specificity", ylab = "sensitivity")
plot(roc_VEGFA,  col = cols["VEGFA"],  lwd = 1.4, add = TRUE)
plot(roc_CXCL10, col = cols["CXCL10"], lwd = 1.4, add = TRUE)
plot(roc_IL18,   col = cols["IL18"],   lwd = 1.4, add = TRUE)
plot(roc_TGFB1,  col = cols["TGFB1"],  lwd = 1.4, add = TRUE)
abline(0, 1, lty = 3, col = "gray50")

legend("bottomright",
       legend = c(
         sprintf("Our model (AUC = %.3f)", auc_myModel),
         sprintf("VEGFA (AUC = %.3f)",     auc_VEGFA),
         sprintf("CXCL10 (AUC = %.3f)",    auc_CXCL10),
         sprintf("IL18 (AUC = %.3f)",      auc_IL18),
         sprintf("TGFB1 (AUC = %.3f)",     auc_TGFB1)
       ),
       col = cols, lwd = 2.8, cex = 0.6, bty = "n")
par(op)
dev.off()
## ========== 1) 其他基因基于最优阈值(closest.topleft) + 指标汇总 ==========
#
get_metrics <- function(roc_obj, name, prob, true_class){
if(name == "myModel"){
	best_threshold = 0.5
	ACC = 0.972
	F1 = 0.982
	Sensitivity = 0.988
	Specificity = 0.926
	Precision = 0.976
	Recall = 0.988
	bACC = 0.957
} else {
co <- coords(roc_obj, x = "best", best.method = "closest.topleft",
               ret = c("threshold","sensitivity","specificity","ppv","npv"),
               transpose = FALSE)
best_threshold <- co["threshold"] %>% as.numeric()
pred_class <- ifelse(prob > best_threshold, 1, 0)
pred_class <- factor(pred_class, levels = c(0, 1))
true_class <- factor(true_class, levels = c(0, 1))
cm <- confusionMatrix(pred_class, true_class, positive = "1")
ACC  <- as.numeric(signif(cm$overall["Accuracy"], 3))
F1   <- as.numeric(signif(cm$byClass["F1"], 3))
Sensitivity <- as.numeric(signif(cm$byClass["Sensitivity"], 3))
Specificity <- as.numeric(signif(cm$byClass["Specificity"], 3))
Precision <- as.numeric(signif(cm$byClass["Precision"], 3))
Recall <- as.numeric(signif(cm$byClass["Recall"], 3))
bACC <- as.numeric(signif(cm$byClass["Balanced Accuracy"], 3))
}
data.frame(
    Curve = name,
    AUC = as.numeric(auc(roc_obj)),
    threshold = best_threshold,
    sensitivity = Sensitivity,
    specificity = Specificity,
	ACC = ACC,
	F1 = F1,
	Precision = Precision,
	Recall = Recall,
	bACC = bACC,
    stringsAsFactors = FALSE
  )
}

metrics_tbl <- do.call(rbind, list(
  get_metrics(roc_myModel, "myModel", NA, NA),
  get_metrics(roc_VEGFA,   "VEGFA", dat_dis["P15692", ] %>% as.numeric(), true_class),
  get_metrics(roc_CXCL10,  "CXCL10", dat_dis["P02778", ] %>% as.numeric(), true_class),
  get_metrics(roc_IL18,    "IL18", dat_dis["Q14116", ] %>% as.numeric(), true_class),
  get_metrics(roc_TGFB1,   "TGFB1", dat_dis["P01137", ] %>% as.numeric(), true_class)
))

## -------- 使用 signif() 控制有效数字 --------
metrics_tbl$AUC         <- signif(metrics_tbl$AUC, 3)
metrics_tbl$threshold   <- signif(metrics_tbl$threshold, 3)
metrics_tbl$sensitivity <- signif(metrics_tbl$sensitivity, 3)
metrics_tbl$specificity <- signif(metrics_tbl$specificity, 3)
metrics_tbl$ACC <- signif(metrics_tbl$ACC, 3)
metrics_tbl$F1 <- signif(metrics_tbl$F1, 3)
metrics_tbl$Precision <- signif(metrics_tbl$Precision, 3)
metrics_tbl$Recall <- signif(metrics_tbl$Recall, 3)
metrics_tbl$bACC <- signif(metrics_tbl$bACC, 3)


print(metrics_tbl)

## 保存结果
write.csv(metrics_tbl, "evaluation_dis.csv", row.names = FALSE)
##可视化
# 假设 metrics_tbl 中已有列：Curve, ACC, sensitivity, specificity
metrics_tbl_long <- metrics_tbl %>%
  select(Curve, ACC, sensitivity, specificity) %>%
  pivot_longer(cols = c(ACC, sensitivity, specificity),
               names_to = "Metric", values_to = "Value")

# 调整顺序（让ACC显示在最前）
metrics_tbl_long$Metric <- factor(metrics_tbl_long$Metric,
                                  levels = c("ACC", "sensitivity", "specificity"))

# 颜色可自定义
metric_colors <- c(
  "ACC" = "#1f77b4",
  "sensitivity" = "#ff7f0e",
  "specificity" = "#2ca02c"
)

# 绘图
metrics_tbl_long$Curve <- factor(metrics_tbl_long$Curve, levels = c("myModel", "VEGFA", "CXCL10","TGFB1","IL18"))
p <- ggplot(metrics_tbl_long,
            aes(x = Curve, y = Value, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = signif(Value, 3)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.3) +
  scale_fill_manual(values = metric_colors) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(title = "dis set",
       x = NULL, y = "Value") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 15, hjust = 1))
pdf("ACC_Spe_Sen_dis.pdf")
print(p)
dev.off()
#####################replication set
#our model
res <- read.csv("xgboost_rep_result_top7.csv")
true_class <- ifelse(res$real_data=="Control", 0, 1)
prob <- res$prediction
roc_myModel <- roc(true_class, prob)  # 计算 ROC 曲线
#VEGFA
exp_value <- dat_rep["P15692", ] %>% as.numeric()
true_class <- ifelse(grepl("Control",colnames(dat_rep)), 0, 1)
roc_VEGFA <- roc(true_class, exp_value)  # 计算 ROC 曲线
#CXCL10
exp_value <- dat_rep["P02778", ] %>% as.numeric()
true_class <- ifelse(grepl("Control",colnames(dat_rep)), 0, 1)
roc_CXCL10 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#IL18
exp_value <- dat_rep["Q14116", ] %>% as.numeric()
true_class <- ifelse(grepl("Control",colnames(dat_rep)), 0, 1)
roc_IL18 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#TGFB1
exp_value <- dat_rep["P01137", ] %>% as.numeric()
true_class <- ifelse(grepl("Control",colnames(dat_rep)), 0, 1)
roc_TGFB1 <- roc(true_class, exp_value)  # 计算 ROC 曲线
#

# 计算 AUC
auc_myModel <- auc(roc_myModel)
auc_VEGFA   <- auc(roc_VEGFA)
auc_CXCL10  <- auc(roc_CXCL10)
auc_IL18    <- auc(roc_IL18)
auc_TGFB1   <- auc(roc_TGFB1)

# 颜色（可按需调整）
cols <- c(
  myModel = "#1f77b4",
  VEGFA   = "#ff7f0e",
  CXCL10  = "#2ca02c",
  IL18    = "#d62728",
  TGFB1   = "#9467bd"
)

# 绘图
pdf("roc_rep.pdf", width = 4, height = 4)
op <- par(mar = c(4.2, 4.2, 3.2, 1.2))
plot(roc_myModel, col = cols["myModel"], lwd = 1.4, legacy.axes = TRUE,
	 main = "rep Set",
     xlab = "1-specificity", ylab = "sensitivity")
plot(roc_VEGFA,  col = cols["VEGFA"],  lwd = 1.4, add = TRUE)
plot(roc_CXCL10, col = cols["CXCL10"], lwd = 1.4, add = TRUE)
plot(roc_IL18,   col = cols["IL18"],   lwd = 1.4, add = TRUE)
plot(roc_TGFB1,  col = cols["TGFB1"],  lwd = 1.4, add = TRUE)
abline(0, 1, lty = 3, col = "gray50")

legend("bottomright",
       legend = c(
         sprintf("Our model (AUC = %.3f)", auc_myModel),
         sprintf("VEGFA (AUC = %.3f)",     auc_VEGFA),
         sprintf("CXCL10 (AUC = %.3f)",    auc_CXCL10),
         sprintf("IL18 (AUC = %.3f)",      auc_IL18),
         sprintf("TGFB1 (AUC = %.3f)",     auc_TGFB1)
       ),
       col = cols, lwd = 2.8, cex = 0.6, bty = "n")
par(op)
dev.off()
## ========== 1) 根据训练的阈值确定rep集的阈值
get_metrics <- function(roc_obj, name, prob, true_class){
if(name == "myModel"){
	best_threshold = 0.5
	ACC = 0.839
	F1 = 0.783
	Sensitivity = 0.818
	Specificity = 0.850
	Precision = 0.750
	Recall = 0.818
	bACC = 0.834
} else {
best_threshold <- metrics_tbl[metrics_tbl$Curve==name, "threshold"] %>% as.numeric()
pred_class <- ifelse(prob > best_threshold, 1, 0)
pred_class <- factor(pred_class, levels = c(0, 1))
true_class <- factor(true_class, levels = c(0, 1))
cm <- confusionMatrix(pred_class, true_class, positive = "1")
ACC  <- as.numeric(signif(cm$overall["Accuracy"], 3))
F1   <- as.numeric(signif(cm$byClass["F1"], 3))
Sensitivity <- as.numeric(signif(cm$byClass["Sensitivity"], 3))
Specificity <- as.numeric(signif(cm$byClass["Specificity"], 3))
Precision <- as.numeric(signif(cm$byClass["Precision"], 3))
Recall <- as.numeric(signif(cm$byClass["Recall"], 3))
bACC <- as.numeric(signif(cm$byClass["Balanced Accuracy"], 3))
}
data.frame(
    Curve = name,
    AUC = as.numeric(auc(roc_obj)),
    threshold = best_threshold,
    sensitivity = Sensitivity,
    specificity = Specificity,
	ACC = ACC,
	F1 = F1,
	Precision = Precision,
	Recall = Recall,
	bACC = bACC,
    stringsAsFactors = FALSE
  )
}
metrics_tbl <- do.call(rbind, list(
  get_metrics(roc_myModel, "myModel", NA, NA),
  get_metrics(roc_VEGFA,   "VEGFA", dat_rep["P15692", ] %>% as.numeric(), true_class),
  get_metrics(roc_CXCL10,  "CXCL10", dat_rep["P02778", ] %>% as.numeric(), true_class),
  get_metrics(roc_IL18,    "IL18", dat_rep["Q14116", ] %>% as.numeric(), true_class),
  get_metrics(roc_TGFB1,   "TGFB1", dat_rep["P01137", ] %>% as.numeric(), true_class)
))

## -------- 使用 signif() 控制有效数字 --------
metrics_tbl$AUC         <- signif(metrics_tbl$AUC, 3)
metrics_tbl$threshold   <- signif(metrics_tbl$threshold, 3)
metrics_tbl$sensitivity <- signif(metrics_tbl$sensitivity, 3)
metrics_tbl$specificity <- signif(metrics_tbl$specificity, 3)
metrics_tbl$ACC <- signif(metrics_tbl$ACC, 3)
metrics_tbl$F1 <- signif(metrics_tbl$F1, 3)
metrics_tbl$Precision <- signif(metrics_tbl$Precision, 3)
metrics_tbl$Recall <- signif(metrics_tbl$Recall, 3)
metrics_tbl$bACC <- signif(metrics_tbl$bACC, 3)

print(metrics_tbl)

## 保存结果
write.csv(metrics_tbl, "evaluation_rep.csv", row.names = FALSE)
##可视化
# 假设 metrics_tbl 中已有列：Curve, ACC, sensitivity, specificity
metrics_tbl_long <- metrics_tbl %>%
  select(Curve, ACC, sensitivity, specificity) %>%
  pivot_longer(cols = c(ACC, sensitivity, specificity),
               names_to = "Metric", values_to = "Value")

# 调整顺序（让ACC显示在最前）
metrics_tbl_long$Metric <- factor(metrics_tbl_long$Metric,
                                  levels = c("ACC", "sensitivity", "specificity"))

# 颜色可自定义
metric_colors <- c(
  "ACC" = "#1f77b4",
  "sensitivity" = "#ff7f0e",
  "specificity" = "#2ca02c"
)

# 绘图
metrics_tbl_long$Curve <- factor(metrics_tbl_long$Curve, levels = c("myModel", "VEGFA", "CXCL10","TGFB1","IL18"))
p <- ggplot(metrics_tbl_long,
            aes(x = Curve, y = Value, fill = Metric)) +
  geom_col(position = position_dodge(width = 0.7), width = 0.6) +
  geom_text(aes(label = signif(Value, 3)),
            position = position_dodge(width = 0.7),
            vjust = -0.5, size = 3.3) +
  scale_fill_manual(values = metric_colors) +
  coord_cartesian(ylim = c(0, 1.05)) +
  labs(title = "rep set",
       x = NULL, y = "Value") +
  theme_minimal(base_size = 14) +
  theme(legend.position = "top",
        axis.text.x = element_text(angle = 15, hjust = 1))
pdf("ACC_Spe_Sen_rep.pdf")
print(p)
dev.off()