setwd("~/LM/Figure3")
library(NCmisc) #p.to.Z
#z-score calculation
degs_dis <- read.csv("~/LM/dif/discovery_DEGs.csv", header= T)
degs_dis <- degs_dis[degs_dis$direction!="NOT", ]
degs_BC <- read.csv("/LM/dif/discovery_breast_DEGs.csv", header = T)
degs_BC <- degs_BC[degs_BC$UniProt%in%degs_dis$UniProt, ]
degs_BC <- degs_BC[match(degs_dis$UniProt, degs_BC$UniProt), ]
degs_LA <- read.csv("/LM/dif/discovery_lung_DEGs.csv", header = T)
degs_LA <- degs_LA[degs_LA$UniProt%in%degs_dis$UniProt, ]
degs_LA <- degs_LA[match(degs_dis$UniProt, degs_LA$UniProt), ]

#
degs_dis$z_score <- p.to.Z(degs_dis$adj.P.Val) * (ifelse(degs_dis$logFC>0, 1, -1))
degs_BC$z_score <-  p.to.Z(degs_BC$adj.P.Val) * (ifelse(degs_BC$logFC>0, 1, -1))
degs_LA$z_score <-  p.to.Z(degs_LA$adj.P.Val) * (ifelse(degs_LA$logFC>0, 1, -1))
z_score <- data.frame(ID = degs_dis$UniProt, Name = degs_dis$Assay, dis = degs_dis$z_score, BC = degs_BC$z_score, LA = degs_LA$z_score)

#All和breast
# 计算相关性
cor_test <- cor.test(z_score$dis, z_score$BC)
r_value <- signif(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
label_text <- paste0("r = ", r_value, ", p = ", p_value)
#
pdf("All_Breast_Zscore_compare.pdf", width = 4, height = 4)
# 绘图
ggplot(z_score, aes(x = dis, y = BC)) +
  geom_point(color = "#1f77b4", alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3, color = "black") +
  theme_classic() +
  xlab("LM vs Control") +
  ylab("LM (BC) vs Control") +
  ggtitle("Comparison of Z-score") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  annotate("text", x = min(z_score$dis), y = max(z_score$BC),
           label = label_text, hjust = 0, vjust = 1, size = 5, color = "black")
dev.off()
#Lung和Breast
# 计算相关性
cor_test <- cor.test(z_score$BC, z_score$LA)
r_value <- signif(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
label_text <- paste0("r = ", r_value, ", p = ", p_value)
#
pdf("Breast_Lung_Zscore_compare.pdf", width = 4, height = 4)
# 绘图
ggplot(z_score, aes(x = BC, y = LA)) +
  geom_point(color = "#1f77b4", alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3, color = "black") +
  theme_classic() +
  xlab("LM (BC) vs Control") +
  ylab("LM (LA) vs Control") +
  ggtitle("Comparison of Z-score") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  annotate("text", x = min(z_score$BC), y = max(z_score$LA),
           label = label_text, hjust = 0, vjust = 1, size = 5, color = "black")
dev.off()
#All和Lung
# 计算相关性
cor_test <- cor.test(z_score$dis, z_score$LA)
r_value <- signif(cor_test$estimate, 3)
p_value <- signif(cor_test$p.value, 3)
label_text <- paste0("r = ", r_value, ", p = ", p_value)
#
pdf("All_Lung_Zscore_compare.pdf", width = 4, height = 4)
# 绘图
ggplot(z_score, aes(x = dis, y = LA)) +
  geom_point(color = "#1f77b4", alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = TRUE) +
  geom_vline(xintercept = 0, linetype = "dashed", size = 0.3, color = "black") +
  geom_hline(yintercept = 0, linetype = "dashed", size = 0.3, color = "black") +
  theme_classic() +
  xlab("LM vs Control") +
  ylab("LM (LA) vs Control") +
  ggtitle("Comparison of Z-score") +
  theme(
    plot.title = element_text(hjust = 0.5, size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  ) +
  annotate("text", x = min(z_score$dis), y = max(z_score$LA),
           label = label_text, hjust = 0, vjust = 1, size = 5, color = "black")
dev.off()

