#直观展示基于前n个基因，xgboost模型的性能。
library(pROC)
library(dplyr)
library(ggplot2)
library(caret)
setwd("~/LM/biomarker_panel")

#计算ACC,F1,敏感性和特异性指标的函数
calc_metrics <- function(true, prob, positive = "LM", threshold = 0.5) {
  
  pred_class <- factor(
    ifelse(prob >= threshold, positive, setdiff(levels(true), positive)),
    levels = levels(true)
  )
  
  cm <- confusionMatrix(pred_class, true, positive = positive)
  
  data.frame(
    ACC  = cm$overall["Accuracy"],
    Sens = cm$byClass["Sensitivity"],
    Spec = cm$byClass["Specificity"],
    F1   = cm$byClass["F1"]
  )
}
#计算AUC和95% CI

auc_ci_df <- data.frame()
auc_ci_df_rep <- data.frame()
n_seq <- seq(1,22,1)
for(n in n_seq)
{	
	load(paste0("xgboost_model_top",n,".rds"))
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
	ungroup()
	#train set
	roc_n <- roc(
    response  = res$obs,
    predictor = res$LM,
    levels = c("Control", "LM"),
    direction = "<"
    )

    ci_n <- ci.auc(roc_n, method = "delong")
	# 其他指标（train）
    metrics_train <- calc_metrics(
    true = res$obs,
    prob = res$LM,
    positive = "LM",
    threshold = 0.5
    )
    auc_ci_df <- rbind(
    auc_ci_df,
    data.frame(
      n = n,
      auc = as.numeric(auc(roc_n)),
      ci_low = ci_n[1],
      ci_high = ci_n[3],
	  ACC  = metrics_train$ACC,
    F1   = metrics_train$F1,
    Sens = metrics_train$Sens,
    Spec = metrics_train$Spec
    )
    )
	##rep
	res <- read.csv(paste0("xgboost_rep_result_top",n,".csv"))
	roc_n <- roc(
    response  = res$real_data,
    predictor = res$prediction,
    levels = c("Control", "LM"),
    direction = "<"
    )

    ci_n <- ci.auc(roc_n, method = "delong")
	# 其他指标（rep）
  metrics_rep <- calc_metrics(
  true = factor(res$real_data, levels = c("Control", "LM")),
  prob = res$prediction,
  positive = "LM",
  threshold = 0.5
  )
    auc_ci_df_rep <- rbind(
    auc_ci_df_rep,
    data.frame(
      n = n,
      auc = as.numeric(auc(roc_n)),
      ci_low = ci_n[1],
      ci_high = ci_n[3],
	   ACC  = metrics_rep$ACC,
    F1   = metrics_rep$F1,
    Sens = metrics_rep$Sens,
    Spec = metrics_rep$Spec
    )
    )
}
#
auc_ci_df$auc_rep <- auc_ci_df_rep$auc
auc_ci_df$ci_low_rep <- auc_ci_df_rep$ci_low
auc_ci_df$ci_high_rep <- auc_ci_df_rep$ci_high
auc_ci_df$ACC_rep <- auc_ci_df_rep$ACC
auc_ci_df$F1_rep <- auc_ci_df_rep$F1
auc_ci_df$Sens_rep <- auc_ci_df_rep$Sens
auc_ci_df$Spec_rep <- auc_ci_df_rep$Spec
plot_df <- auc_ci_df


write.csv(plot_df, file = "all_metric.csv",row.names = F)
library(ggplot2)


pdf("auc_train_only.pdf", width = 8, height = 5)

# 左轴范围（只用训练集数据）
auc_min <- min(plot_df$auc)
auc_max <- max(plot_df$auc)
ACC_min <- min(plot_df$ACC)
ACC_max <- max(plot_df$ACC)
scale_factor <- (auc_max - auc_min) / (ACC_max - ACC_min)

# 转长表格，只保留训练集
library(tidyr)
plot_long <- plot_df %>%
  pivot_longer(
    cols = c(auc, ACC),
    names_to = "metric",
    values_to = "value"
  ) %>%
  mutate(
    type = ifelse(metric == "auc", "AUC", "ACC"),
    dataset = "Training",  # 全部是训练集
    value_scaled = ifelse(type == "ACC", (value - ACC_min) * scale_factor + auc_min, value)
  )


ggplot(plot_long, aes(x = n, y = value, color = type)) +
  geom_point(aes(shape = type), size = 2) +
  geom_line(linewidth = 1) +
  geom_text(
    aes(label = sprintf("%.3f", value)),
    vjust = ifelse(plot_long$type == "AUC", -0.8, 1.5),
    size = 2.5,
    show.legend = FALSE
  ) +
  theme_classic() +
  scale_color_manual(values = c("AUC" = "#E64B35","ACC" = "#3C5488")) +
  scale_linetype_manual(values = c("AUC" = "solid", "ACC" = "dashed")) +
  scale_shape_manual(values = c("AUC" = 16, "ACC" = 17)) +
  labs(x = "Number of features (n)", color = "Metric", linetype = "Metric", shape = "Metric") +
  scale_x_continuous(breaks = plot_df$n) +
  theme_classic() +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )
dev.off()