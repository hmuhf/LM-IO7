library(openxlsx) #read.xlsx
library(ggplot2)
library(dplyr)
library(tidyr)
#
setwd("~/LM/biomarker_panel")

#read the patient info
sample_info <- read.xlsx("~/LM/data/patients_clinical_information.xlsx", startRow = 2)
case <- sample_info[, c("Patient","Cytology","MRI")]
#train set
sample_train <- patient_info[patient_info$Cohort=="Discovery" & patient_info$Group=="LM", ]
sample_train <- merge(sample_train, case, by = "Patient", all = F)
train_cytology_sensitivity <- mean(sample_train[,"Cytology"] == "Positive")
train_MRI_sensitivity <- mean(sample_train[,"MRI"] == "Positive")

#replication
sample_rep <- patient_info[patient_info$Cohort=="Relication" & patient_info$Group=="LM", ]
sample_rep <- merge(sample_rep, case, by = "Patient", all = F)
rep_cytology_sensitivity <- mean(sample_rep[,"Cytology"] == "Positive")
rep_MRI_sensitivity <- mean(sample_rep[,"MRI"] == "Positive")

# 创建数据框
metrics_data <- data.frame(
  Dataset = rep(c("Training", "Replication"), each = 3),
  Method = rep(c("LM-IO7","Cytology", "MRI"), times = 2),
  Sensitivity = c(0.988,train_cytology_sensitivity, train_MRI_sensitivity, 0.818,
                 rep_cytology_sensitivity, rep_MRI_sensitivity)
)
metrics_data$Dataset <- factor(metrics_data$Dataset, levels = c("Training", "Replication"))
metrics_data$Method <- factor(metrics_data$Method, levels = c("LM-IO7", "Cytology","MRI"))
# 自定义颜色
colors <- c(
  "Training" = "#1f77b4",  # 训练集用蓝色
  "Replication" = "#ff7f0e"    # 测试集用橙色
)

# 绘制直方图
pdf("Cytology_MRI_sensitivity.pdf")
p1 <- ggplot(metrics_data, aes(x = Method, y = Sensitivity, fill = Dataset)) +
  geom_col(position = position_dodge(width = 0.8), width = 0.7) +
  geom_text(aes(label = sprintf("%.2f", Sensitivity)), 
            position = position_dodge(width = 0.8),
            vjust = -0.5, size = 4) +
  scale_fill_manual(values = colors) +
  scale_y_continuous(limits = c(0, 1.05), expand = c(0, 0)) +
  labs(
    title = "Sensitivity: Cytology and MRI",
    x = "Method",
    y = "Sensitivity",
    fill = "Dataset"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(size = 12),
    panel.grid.major.x = element_blank()
  )
print(p1)
dev.off()