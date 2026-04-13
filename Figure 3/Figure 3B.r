library(pROC) #coords
library(caret) #confusionMatrix
library(ggplot2)
library(reshape2)  
library(dplyr)
library(tidyr)
library(pROC)
#
setwd("~/LM/biomarker_panel")
method <- c("enet","lasso","ridge","xgboost","gbm","rf","knn","svm","nb")

index_train <- data.frame()
best_threshold_train <- c()

# color
method_colors <- setNames(
  rainbow(length(method)),
  method
)



pdf("roc_train.pdf", width = 4, height = 4)

first_plot <- TRUE  # 
auc_list <- numeric(length(method))
names(auc_list) <- method

for (i in method) {

  res <- read.csv(paste0(i, "_train_cv_result.csv"))

  # 
  res <- res %>%
    group_by(rowIndex, obs) %>%
    summarise(
      LM = mean(LM, na.rm = TRUE),
      .groups = "drop"
    )

  # 
  true_class <- ifelse(res$obs == "Control", 0, 1)
  prob <- res$LM

  # ROC & AUC
  roc_obj <- roc(true_class, prob)
  auc_val <- as.numeric(signif(auc(roc_obj), 3))
  auc_list[i] <- auc_val

  # ===============================
  # threshold
  # ===============================
  best_threshold <- 0.5
  pred_class <- ifelse(prob > best_threshold, 1, 0)

  pred_class <- factor(pred_class, levels = c(0, 1))
  true_class <- factor(true_class, levels = c(0, 1))

  cm <- confusionMatrix(pred_class, true_class, positive = "1")

  ACC  <- as.numeric(signif(cm$overall["Accuracy"], 3))
  F1   <- as.numeric(signif(cm$byClass["F1"], 3))
  Sens <- as.numeric(signif(cm$byClass["Sensitivity"], 3))
  Spec <- as.numeric(signif(cm$byClass["Specificity"], 3))
  # ===============================
  # save
  # ===============================
  index_train <- rbind(
    index_train,
    data.frame(
      Method = i,
      AUC = auc_val,
      ACC = ACC,
      F1 = F1,
      Sensitivity = Sens,
      Specificity = Spec
    )
  )

  # ===============================
  # plot
  # ===============================
  if (first_plot) {
    plot(
      roc_obj,
      col = method_colors[i],
      legacy.axes = TRUE,
      xlab = "1 - Specificity",
      ylab = "Sensitivity",
      main = "Training ROC Curves",
      lwd = 2
    )
    first_plot <- FALSE
  } else {
    lines(
      roc_obj,
      col = method_colors[i],
      lwd = 2
    )
  }
}

## ====== sort AUC======
ord <- order(index_train$AUC, decreasing = TRUE)

legend(
  "bottomright",
  legend = paste0(
    index_train$Method[ord],
    " (AUC = ",
    index_train$AUC[ord],
    ")"
  ),
  col = method_colors[index_train$Method[ord]],
  lwd = 2,
  cex = 0.8,
  bty = "n"
)

dev.off()

colnames(index_train) <- paste0("train_", colnames(index_train))
###柱状图可视化所有指标
# =========================
# 1. 整理数据为长表格
# =========================
train_df <- index_train %>%
  select(train_Method,train_AUC, train_ACC) %>%
  rename(Method = train_Method) %>%
  pivot_longer(
    cols = -Method,
    names_to = "Metric",
    values_to = "Value"
  )

# =========================
# 2. 设置 Method 顺序（xgboost 第一个）
# =========================
method_order <- c("xgboost","gbm","nb","enet","rf","svm","knn","lasso","ridge")
train_df$Method <- factor(train_df$Method, levels = method_order)
method_cols <- c(
  "xgboost" = "#DC0000",  # 主模型用红色
  "gbm"     = "#E64B35",
  "nb"      = "#F39B7F",
  "enet"    = "#00A087",
  "rf"      = "#4DBBD5",
  "svm"     = "#4575B4",
  "knn"     = "#3C5488",
  "lasso"   = "#762A83",
  "ridge"   = "#7E6148"
)

# =========================
# 3. 画柱状图
# =========================
pdf("All_metric_train.pdf", width = 8, height = 4)
ggplot(train_df, aes(x = Metric, y = Value, fill = Method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.85), width = 0.7) +
  geom_text(aes(label = Value),
            position = position_dodge(width = 0.85),
            vjust = -0.3,
            size = 3) +
  scale_y_continuous(limits = c(0,1), expand = expansion(mult = c(0, 0.05))) +
  coord_cartesian(ylim = c(0.75, 1)) +
  labs(x = "Method", y = "Metric Value", fill = "Metric") +
  scale_fill_manual(values = method_cols)+
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 0, hjust = 0.5)
  )
dev.off()