library(openxlsx)
library(survival)
library(survminer)
library(forestplot)
library(dplyr)
library(grid)
library(pROC)
library(timeROC)
setwd("~/LM/application")

#read the survival info
patient_details <-read.xlsx("~/LM/data/patients_clinical_information.xlsx", startRow = 2)
survival_info <- patient_details[, c("Patient", "iPFS_time","iPFS_status")]
#survival_info <- patient_details[, c("Patient", "OS_time","OS_status")]
survival_info <- na.omit(survival_info)
print(dim(survival_info))

#read the expression value
dat_overall <- read.table("~/LM/data/dat_all_exp.txt", sep = "\t", header = T, as.is = T, check.names = F)
select_gene <- read.table("~/LM/biomarker_panel/SHAP_importance.txt", header = T)
select_gene <- select_gene$Feature
dat_overall <- dat_overall[select_gene, ]

dat <- as.data.frame(t(dat_overall))
dat$sample <- rownames(dat)

#######overall cohort
#covariates
covariates <- colnames(dat)[-which(colnames(dat)%in%"sample")]
print(length(covariates))
survival_pfs <- merge(survival_info, dat, by.x = "Patient", by.y = "sample", all = F)

surv_formulas <- sapply(covariates, function(x){as.formula(paste0("Surv(iPFS_time, iPFS_status) ~ ", x))}) 
res <- lapply(surv_formulas, function(x){coxph(x, dat = survival_pfs)})
res_frame <- lapply(res, function(x){
	tmp <- summary(x)
	#p value
	p_value <- signif(tmp$wald['pvalue'], 2)
	#HR
	HR <- signif(tmp$coef[2], 2)
	#CI
	HR_confint_lower <- signif(tmp$conf.int[, 'lower .95'], 2)
	HR_confint_upper <- signif(tmp$conf.int[, 'upper .95'], 2)
	res <- c(p_value, HR, HR_confint_lower, HR_confint_upper)
	names(res) <- c("p_value", "HR", "CI_lower", "CI_upper")
	return(res)
})
data <- t(as.data.frame(res_frame)) %>% as.data.frame()
#
data$protein <- rownames(data)
print(dim(data))
protein_id_name <- read.table("~/LM/data/protein_id_name.txt", header= T, sep ="\t")
data <- merge(data, protein_id_name, by.x = "protein",by.y="UniProt",all.x = T, all.y = F)
print(head(data))
#
data <- data %>%
  mutate(p_value = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)),
         p_value = ifelse(p_value < 0.05, paste0(p_value, "*"), p_value))
print(head(data))
labeltext <- cbind(
  c("Variable", data$Assay),  
  c("HR", data$HR),             
  c("95% CI", paste0("(", data$CI_lower, "-", data$CI_upper, ")")),
  c("p-value", data$p_value)    
)
# 
pdf("xgboost_ipfs_forestplot_overall.pdf", width = 6, height = 5)

forestplot(
  labeltext = labeltext,
  mean  = c(NA, data$HR),
  lower = c(NA, data$CI_lower),
  upper = c(NA, data$CI_upper),

  zero = 1,
  xlog = TRUE,
  clip = c(0.5, 2.0),

  is.summary = c(TRUE, rep(FALSE, nrow(data))),
  graph.pos = 2,

  col = fpColors(
    box = "#2C7BB6",        
    line = "#2C7BB6",       
    summary = "#D7191C",    
    zero = "black"          
  ),

  # 
  boxsize = 0.18,
  lwd.ci = 2,              
  lwd.zero = 1.5,           
  ci.vertices = TRUE,       
  ci.vertices.height = 0.1,

  # 

  xticks = c(0.5, 0.75, 1, 1.5, 2.0),
  xticks.cex = 0.9,

  
  txt_gp = fpTxtGp(
    label = gpar(fontsize = 11),
    ticks = gpar(fontsize = 10),
    xlab  = gpar(fontsize = 12, fontface = "bold"),
    title = gpar(fontsize = 13, fontface = "bold")
  ),

  
  lineheight = unit(0.7, "cm"),
  colgap = unit(6, "mm"),

  
  graphwidth = unit(60, "mm"),

 
  title = "Hazard Ratio for iPFS (overall)",

  xlab = "Hazard Ratio (log scale)"
)
dev.off()
#KM curve for IL7 (within overall cohorts)
protein <- "P13232"
protein_value <- as.numeric(survival_pfs[, protein]) 
protein_name <- protein_id_name$Assay[protein_id_name$UniProt==protein]
group <- ifelse(protein_value>median(protein_value), "high_risk","low_risk")
group <- factor(group, levels = c("low_risk","high_risk"))
risk_score <- data.frame(Sample = survival_info$Patient, group)
surv_info <- merge(survival_info, risk_score, by.x = "Patient", by.y = "Sample", all = F)
#
fit <- survfit(Surv(iPFS_time, iPFS_status) ~ group, data = surv_info)
colour_man <- c("#1d52a1","#f3878d")
title_name <- paste0("All: ", protein_name)
pdf(paste0(protein_name, "_ipfs_KM_overall.pdf"),width = 4, height = 5)
p1 <-  ggsurvplot(
           data = surv_info,
           fit,
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           surv.median.line = "hv",
		   ggtheme = theme_bw() + theme(panel.grid = element_blank()), 
           palette = colour_man,
           title = title_name,
           font.main = 10,     
           font.x = 10,        
           font.y = 10,        
           font.tickslab = 10,
           risk.table.height = 0.3,
           pval.size = 3,
           risk.table.fontsize = 4,
           font.legend = 8,
		   size = 0.5
           )
print(p1)
dev.off()
#ROC curve for IL7 (within overall cohorts)
time_roc_res <- timeROC(
  T = survival_pfs$iPFS_time,
  delta = survival_pfs$iPFS_status,
  marker = as.numeric(survival_pfs[,protein]),
  cause = 1,
  weighting = "marginal",
  times = c(0.5 * 12, 1 * 12, 1.5 * 12),
  ROC = TRUE,
  iid = TRUE
)
roc_df <- data.frame(
  FPR_1 = time_roc_res$FP[,1],
  TPR_1 = time_roc_res$TP[,1],
  FPR_2 = time_roc_res$FP[,2],
  TPR_2 = time_roc_res$TP[,2],
  FPR_3 = time_roc_res$FP[,3],
  TPR_3 = time_roc_res$TP[,3]
)

# AUC
auc_vals <- round(time_roc_res$AUC[1:3], 3)
labels <- c(
  paste0("0.5 Years (AUC=", auc_vals[1], ")"),
  paste0("1 Years (AUC=", auc_vals[2], ")"),
  paste0("1.5 Years (AUC=", auc_vals[3], ")")
)

# 
roc_long <- dplyr::bind_rows(
  data.frame(FPR = time_roc_res$FP[,1], TPR = time_roc_res$TP[,1], Time = labels[1]),
  data.frame(FPR = time_roc_res$FP[,2], TPR = time_roc_res$TP[,2], Time = labels[2]),
  data.frame(FPR = time_roc_res$FP[,3], TPR = time_roc_res$TP[,3], Time = labels[3])
)

# 
roc_long$Time <- factor(roc_long$Time, levels = labels)

# 
cols <- setNames(c("#716db2", "#72c15a", "#f3793b"), labels)

# ggplot
title_name <- protein
pdf(paste0(protein,"_roc_overall_ipfs.pdf"), width = 5, height = 5)
p<-ggplot(roc_long, aes(x = FPR, y = TPR, color = Time)) +
  geom_line(size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = cols) +
  labs(
    title = title_name,
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw(base_size = 14)+
  theme(
    panel.grid = element_blank(), 
    legend.position =  c(0.7, 0.3),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.2)))+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.02)))+
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02)))
print(p)
dev.off()


#######Primary breast cancer or Lung cancer
##################
condition <- c("BC","LA")

for(i in 1:length(condition)){
dat <- dat_overall[, grepl(condition[i], colnames(dat_overall))]
dat <- as.data.frame(t(dat))
dat$sample <- rownames(dat)

#covariates
covariates <- colnames(dat)[-which(colnames(dat)%in%"sample")]
print(length(covariates))
survival_pfs <- merge(survival_info, dat, by.x = "Patient", by.y = "sample", all = F)

surv_formulas <- sapply(covariates, function(x){as.formula(paste0("Surv(iPFS_time, iPFS_status) ~ ", x))})

res <- lapply(surv_formulas, function(x){coxph(x, dat = survival_pfs)})

res_frame <- lapply(res, function(x){
	tmp <- summary(x)
	
	p_value <- signif(tmp$wald['pvalue'], 2)
	#HR
	HR <- signif(tmp$coef[2], 2)
	
	HR_confint_lower <- signif(tmp$conf.int[, 'lower .95'], 2)
	HR_confint_upper <- signif(tmp$conf.int[, 'upper .95'], 2)
	res <- c(p_value, HR, HR_confint_lower, HR_confint_upper)
	names(res) <- c("p_value", "HR", "CI_lower", "CI_upper")
	return(res)
})
data <- t(as.data.frame(res_frame)) %>% as.data.frame()
data$protein <- rownames(data)
data <- merge(data, protein_id_name, by.x = "protein",by.y="UniProt",all.x = T, all.y = F)

data <- data %>%
  mutate(p_value = ifelse(p_value < 0.001, "<0.001", sprintf("%.3f", p_value)),
         p_value = ifelse(p_value < 0.05, paste0(p_value, "*"), p_value))
print(head(data))
labeltext <- cbind(
  c("Variable", data$Assay),  
  c("HR", data$HR),              
  c("95% CI", paste0("(", data$CI_lower, "-", data$CI_upper, ")")), 
  c("p-value", data$p_value)    
)
# 
pdf(paste0("xgboost_ipfs_forestplot_",condition[i],".pdf"),width = 6, height = 5)
if(condition[i]=="BC")
  {xticks_value = c(0.2, 0.5, 1, 1.5, 2,2.5,3.0)}
else {
  xticks_value = c(0.5, 0.75, 1, 1.5, 2)
  }
p1<-
forestplot(
  labeltext = labeltext,
  mean  = c(NA, data$HR),
  lower = c(NA, data$CI_lower),
  upper = c(NA, data$CI_upper),

  zero = 1,
  xlog = TRUE,
  clip = c(0.5, 2.0),

  is.summary = c(TRUE, rep(FALSE, nrow(data))),
  graph.pos = 2,

 
  col = fpColors(
    box = "#2C7BB6",        
    line = "#2C7BB6",       
    summary = "#D7191C",   
    zero = "black"          
  ),

  # 
  boxsize = 0.18,
  lwd.ci = 2,              
  lwd.zero = 1.5,           
  ci.vertices = TRUE,      
  ci.vertices.height = 0.1,

  #
  xticks = xticks_value,
  xticks.cex = 0.9,

  # 
  txt_gp = fpTxtGp(
    label = gpar(fontsize = 11),
    ticks = gpar(fontsize = 10),
    xlab  = gpar(fontsize = 12, fontface = "bold"),
    title = gpar(fontsize = 13, fontface = "bold")
  ),

  
  lineheight = unit(0.7, "cm"),
  colgap = unit(6, "mm"),

  
  graphwidth = unit(60, "mm"),

 
  title = "Hazard Ratio for iPFS",

  xlab = "Hazard Ratio (log scale)"
)
print(p1)
dev.off()

#KM curve for GZMA/IL7 (within specific cancer type)
protein <- "P12544" #P13232
protein_value <- as.numeric(survival_pfs[, protein]) 
protein_name <- protein_id_name$Assay[protein_id_name$UniProt==protein]
group <- ifelse(protein_value>median(protein_value), "high_risk","low_risk")
group <- factor(group, levels = c("low_risk","high_risk"))
risk_score <- data.frame(Sample = survival_pfs$Patient, group)
surv_info <- merge(survival_pfs, risk_score, by.x = "Patient", by.y = "Sample", all = F)
#
fit <- survfit(Surv(iPFS_time, iPFS_status) ~ group, data = surv_info)
colour_man <- c("#1d52a1","#f3878d")
title_name <- paste0(condition[i],":", protein_name)
pdf(paste0(protein_name, "_ipfs_KM_",condition[i],".pdf"),width = 4, height = 5)
p1 <-  ggsurvplot(
           data = surv_info,
           fit,
           pval = TRUE,
           risk.table = TRUE,
           risk.table.col = "strata",
           surv.median.line = "hv",
		   ggtheme = theme_bw() + theme(panel.grid = element_blank()), 
           palette = colour_man,
           title = title_name,
           font.main = 10,     
           font.x = 10,        
           font.y = 10,        
           font.tickslab = 10,
           risk.table.height = 0.3,
           pval.size = 3,
           risk.table.fontsize = 4,
           font.legend = 8,
		   size = 0.5
           )
print(p1)
dev.off()

#ROC curve for GZMA/IL7 (within specific cancer type)
time_roc_res <- timeROC(
  T = survival_pfs$iPFS_time,
  delta = survival_pfs$iPFS_status,
  marker = as.numeric(survival_pfs[,protein]),
  cause = 1,
  weighting = "marginal",
  times = c(0.5 * 12, 1 * 12, 1.5 * 12),
  ROC = TRUE,
  iid = TRUE
)
roc_df <- data.frame(
  FPR_1 = time_roc_res$FP[,1],
  TPR_1 = time_roc_res$TP[,1],
  FPR_2 = time_roc_res$FP[,2],
  TPR_2 = time_roc_res$TP[,2],
  FPR_3 = time_roc_res$FP[,3],
  TPR_3 = time_roc_res$TP[,3]
)

# AUC
auc_vals <- round(time_roc_res$AUC[1:3], 3)
labels <- c(
  paste0("0.5 Years (AUC=", auc_vals[1], ")"),
  paste0("1 Years (AUC=", auc_vals[2], ")"),
  paste0("1.5 Years (AUC=", auc_vals[3], ")")
)

# 
roc_long <- dplyr::bind_rows(
  data.frame(FPR = time_roc_res$FP[,1], TPR = time_roc_res$TP[,1], Time = labels[1]),
  data.frame(FPR = time_roc_res$FP[,2], TPR = time_roc_res$TP[,2], Time = labels[2]),
  data.frame(FPR = time_roc_res$FP[,3], TPR = time_roc_res$TP[,3], Time = labels[3])
)

# 
roc_long$Time <- factor(roc_long$Time, levels = labels)

# 
cols <- setNames(c("#716db2", "#72c15a", "#f3793b"), labels)

# ggplot
title_name <- protein
pdf(paste0(protein,"_roc_",condition[i],"_ipfs.pdf"), width = 5, height = 5)
p<-ggplot(roc_long, aes(x = FPR, y = TPR, color = Time)) +
  geom_line(size = 0.8) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "gray") +
  scale_color_manual(values = cols) +
  labs(
    title = title_name,
    x = "1 - Specificity",
    y = "Sensitivity"
  ) +
  theme_bw(base_size = 14)+
  theme(
    panel.grid = element_blank(), 
    legend.position =  c(0.7, 0.3),
    legend.text = element_text(size = 12)
  ) +
  guides(color = guide_legend(override.aes = list(size = 1.2)))+
  scale_y_continuous(expand = expansion(mult = c(0.01, 0.02)))+
  scale_x_continuous(expand = expansion(mult = c(0.01, 0.02)))
print(p)
dev.off()
}
############drug response comparison
response <- patient_details[, c("Patient", "Response")]
print(dim(response))

#read the expression value
dat_overall <- read.table("~/LM/data/dat_all_exp.txt", sep = "\t", header = T, as.is = T, check.names = F)
select_gene <- read.table("~/LM/biomarker_panel/SHAP_importance.txt", header = T)# select_gene <- select_gene$Feature[select_gene$SHAP_Mean_Scaled!=0]
select_gene <- select_gene$Feature
dat_overall <- dat_overall[select_gene, ]
#All
dat <- as.data.frame(t(dat_overall))
dat$sample <- rownames(dat)
dat <- merge(dat, response, by.x = "sample", by.y = "Patient")
pdf("All_exp_compare.pdf", width = 5, height = 5)
for(protein_id in select_gene){
	pdat <- dat[, c(protein_id,"Response")]
	# Response 顺序
pdat$Response <- factor(
  pdat$Response,
  levels = c("PR", "SD", "PD")
)
comparisons <- list(
  c("PR", "SD"),
  c("PR", "PD"),
  c("SD", "PD")
)


# 颜色（高区分度，论文常用）
resp_cols <- c(
  PR = "#1b9e77",
  SD = "#7570b3",
  PD = "#d95f02"
)
p <- ggplot(pdat, aes(x = Response, y = .data[[protein_id]], color = Response)) +
  geom_boxplot(
    width = 0.6,
    outlier.shape = NA
  ) +
  geom_point(
    size = 2,
    alpha = 0.8
  ) +
  scale_color_manual(values = resp_cols) +
  stat_compare_means(
    comparisons = comparisons,
    method = "wilcox.test",
    label = "p.format"   # 显示 *, **, ***
  ) +
  labs(
    x = "Response",
    y = paste0(protein_id, " expression"),
    color = "Response"
  ) +
  theme_classic()
print(p)
}
dev.off()
#subgroup with primary breast cancer / lung cancer 