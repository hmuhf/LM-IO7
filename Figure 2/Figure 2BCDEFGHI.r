setwd("~/LM/Figure2")
library(openxlsx)
library(dplyr)
###pca
pca_plot <- function(title_name, pdat, group, color_man, sample_name){
        library(ggplot2)
        library(ggrepel)
		pdat <- t(pdat)
        pca <- prcomp(pdat, center = TRUE, scale. = TRUE)
        dat <- as.data.frame(pca$x)
        dat$group <- group
        #
        summ <- summary(pca)
        xlab <- paste0("PCA Dim1(",round(summ$importance[2,1]*100,2),"%)")
        ylab <- paste0("PCA Dim2(",round(summ$importance[2,2]*100,2),"%)")
        

        p1 <- ggplot(data = dat, aes(x = PC1, y = PC2, color = group)) +
			geom_point(size = 2) +
			labs(x = xlab, y = ylab, title = title_name) +
			scale_colour_manual(values = color_man) +
			theme_bw() +
	theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13),
    plot.margin = unit(c(0.4,0.4,0.4,0.4), 'cm'),
    panel.border = element_blank(),          # 
    panel.grid = element_blank(),            # 
    axis.line = element_line(color = "black") # 
  )


        return(p1)
}
#
plsda_plot <- function(title_name, pdat, group, color_man, sample_name){
        library(ggplot2)
        #PLS-DA
		library(mixOmics)
		pdat <- t(pdat)
		pls_da <- plsda(pdat, group, ncomp = 2)
        dat <- as.data.frame(pls_da$variates$X)
        dat$group <- group
        #
        xlab <- paste0("PLS-DA Dim1(",round(pls_da$prop_expl_var$X[1]*100,2),"%)")
        ylab <- paste0("PLS-DA Dim2(",round(pls_da$prop_expl_var$X[2]*100,2),"%)")

        p1 <- ggplot(data = dat, aes(x = comp1, y = comp2, color = group)) +
			geom_point(size = 2) +
			labs(x = xlab, y = ylab, title = title_name) +
			scale_colour_manual(values = color_man) +
			theme_bw() +
    theme(
    plot.title = element_text(hjust = 0.5, size = 15),
    axis.text = element_text(size = 11),
    axis.title = element_text(size = 13),
    legend.text = element_text(size = 11),
    legend.title = element_text(size = 13),
    plot.margin = unit(c(0.4,0.4,0.4,0.4), 'cm'),
    panel.border = element_blank(),          # 
    panel.grid = element_blank(),            # 
    axis.line = element_line(color = "black") # 
  )
        return(p1)
}
#
library(RColorBrewer)
dat <- read.xlsx("~/LM/data/dat_exp_all.xlsx", startRow = 2)
rownames(dat) <- dat$UniProt
dat <- dat[, -c(1,2)]
sample_info <- read.xlsx("~/LM/data/patients_clinical_information.xlsx", startRow = 2)
#Discovery
pdat <- dat[, colnames(dat)%in%sample_info$Patient[sample_info$Cohort=="Discovery"]]
title_name <- "Discovery"
sample <- colnames(pdat)
#
group_info <- data.frame(sample=colnames(pdat), group=ifelse(grepl("BC", sample), "BC", ifelse(grepl("LA", sample), "LA", ifelse(grepl("EC", sample), "EC", ifelse(grepl("MM", sample), "MM", "Control")))))
group <- factor(group_info$group, levels = c("Control", "BC", "LA","EC","MM"))
color_man <- c("#00A087", "#3C5488", "#E64B35", "#4DBBD5", "#F39B7F")
pdf(paste0(title_name, "_pca.pdf"), width = 4.5, height = 4)
p1 <- pca_plot(title_name, pdat, group,color_man, sample)
print(p1)
dev.off()
#pda-da
pdf(paste0(title_name, "_pls_da.pdf"), width = 4.5, height = 4)
p1 <- plsda_plot(title_name, pdat, group,color_man, sample)
print(p1)
dev.off()


#Replication
pdat <- dat[, colnames(dat)%in%sample_info$Patient[sample_info$Cohort=="Replication"]]
title_name <- "Replication"
sample <- colnames(pdat)
#
group_info <- data.frame(sample=colnames(pdat), group=ifelse(grepl("BC", sample), "BC", ifelse(grepl("LA", sample), "LA", ifelse(grepl("EC", sample), "EC", ifelse(grepl("MM", sample), "MM", "Control")))))
group <- factor(group_info$group, levels = c("Control",  "LA"))
library(RColorBrewer)
color_man <- c("#00A087",  "#E64B35")
pdf(paste0(title_name, "_pca.pdf"), width = 4.5, height = 4)
p1 <- pca_plot(title_name, pdat, group,color_man, sample)
print(p1)
dev.off()
#pda-da
pdf(paste0(title_name, "_pls_da.pdf"), width = 4.5, height = 4)
p1 <- plsda_plot(title_name, pdat, group,color_man, sample)
print(p1)
dev.off()

###Correlation analysis between PCA/PLS-DA components and sex,age,group
dat <- read.xlsx("~/LM/data/dat_exp_all.xlsx", startRow = 2)
rownames(dat) <- dat$UniProt
dat <- dat[, -c(1,2)]
sample_info <- read.xlsx("~/LM/data/patients_clinical_information.xlsx", startRow = 2)

#Discovery
pdat <- dat[, colnames(dat)%in%sample_info$Patient[sample_info$Cohort=="Discovery"]]
title_name <- "Discovery"
group_info <- data.frame(sample=colnames(pdat), group=ifelse(grepl("Control", colnames(pdat)), "Control","LM"))
#pca
pdat <- t(pdat)
pca <- prcomp(pdat, center = TRUE, scale. = TRUE)
res <- as.data.frame(pca$x)
res$group <- group_info$group
res$sample <- group_info$sample
res <- merge(res, pheno, by.x = "sample",by.y = "Patient", all.x = T, all.y = F)

result <- c()
for(i in 1:5){
	value <- res[, paste0("PC", i)] %>% as.numeric()
	value_control <- value[which(res$group=="Control")]
	value_case <- value[which(res$group=="LM")]
	p <- wilcox.test(value_control, value_case)
	p_group <- p$p.value
	#
	value_control <- value[which(res$Sex=="Female")]
	value_case <- value[which(res$Sex=="Male")]
	p <- wilcox.test(value_control, value_case)
	p_sex <- p$p.value
	#
	value_x <- res[, paste0("PC", i)] %>% as.numeric()
	value_y <- res[, "Age"] %>% as.numeric()
	p <- cor.test(value_x, value_y, method = c("pearson"))
	p_age <- p$p.value
	tmp <- c(p_group, p_sex, p_age)
	result <- rbind(result, tmp)
}
colnames(result) <- c("condition","sex","age")
rownames(result) <- paste0("PC",1:5)
summ <- summary(pca)
result <- as.data.frame(result)
result$contribution <- paste0(round(summ$importance[2,1:5] *100, 1),"%")
result_dis_pca <- result
#
#Replication
pdat <- dat[, colnames(dat)%in%sample_info$Patient[sample_info$Cohort=="Replication"]]
title_name <- "Replication"
group_info <- data.frame(sample=colnames(pdat), group=ifelse(grepl("Control", colnames(pdat)), "Control","LM"))
#pca
pdat <- t(pdat)
pca <- prcomp(pdat, center = TRUE, scale. = TRUE)
res <- as.data.frame(pca$x)
res$group <- group_info$group
res$sample <- group_info$sample
res <- merge(res, pheno, by.x = "sample",by.y = "Patient", all.x = T, all.y = F)

result <- c()
for(i in 1:5){
	value <- res[, paste0("PC", i)] %>% as.numeric()
	value_control <- value[which(res$group=="Control")]
	value_case <- value[which(res$group=="LM")]
	p <- wilcox.test(value_control, value_case)
	p_group <- p$p.value
	#
	value_control <- value[which(res$Sex=="Female")]
	value_case <- value[which(res$Sex=="Male")]
	p <- wilcox.test(value_control, value_case)
	p_sex <- p$p.value
	#
	value_x <- res[, paste0("PC", i)] %>% as.numeric()
	value_y <- res[, "Age"] %>% as.numeric()
	p <- cor.test(value_x, value_y, method = c("pearson"))
	p_age <- p$p.value
	tmp <- c(p_group, p_sex, p_age)
	result <- rbind(result, tmp)
}
colnames(result) <- c("condition","sex","age")
rownames(result) <- paste0("PC",1:5)
summ <- summary(pca)
result <- as.data.frame(result)
result$contribution <- paste0(round(summ$importance[2,1:5] *100, 1),"%")
result_rep_pca <- result


#PLS-DA
library(mixOmics)
#Discovery
pdat <- dat[, colnames(dat)%in%sample_info$Patient[sample_info$Cohort=="Discovery"]]
title_name <- "Discovery"
group_info <- data.frame(sample=colnames(pdat), group=ifelse(grepl("Control", colnames(pdat)), "Control","LM"))
#pls-da
pdat <- t(pdat)
pls_da <- plsda(pdat, group_info$group, ncomp = 5)
res <- as.data.frame(pls_da$variates$X)
res$group <- group_info$group
res$sample <- group_info$sample
res <- merge(res, pheno, by.x = "sample",by.y = "Patient", all.x = T, all.y = F)

result <- c()
for(i in 1:5){
	value <- res[, paste0("comp", i)] %>% as.numeric()
	value_control <- value[which(res$group=="Control")]
	value_case <- value[which(res$group=="LM")]
	p <- wilcox.test(value_control, value_case)
	p_group <- p$p.value
	#
	value_control <- value[which(res$Sex=="Female")]
	value_case <- value[which(res$Sex=="Male")]
	p <- wilcox.test(value_control, value_case)
	p_sex <- p$p.value
	#
	value_x <- res[, paste0("comp", i)] %>% as.numeric()
	value_y <- res[, "Age"] %>% as.numeric()
	p <- cor.test(value_x, value_y, method = c("pearson"))
	p_age <- p$p.value
	tmp <- c(p_group, p_sex, p_age)
	result <- rbind(result, tmp)
}
colnames(result) <- c("condition","sex","age")
rownames(result) <- paste0("comp",1:5)
summ <- summary(pca)
result <- as.data.frame(result)
result$contribution <- paste0(round(pls_da$prop_expl_var$X[1:5]*100,2),"%")
result_dis_plsda <- result
#
#Replication
pdat <- dat[, colnames(dat)%in%sample_info$Patient[sample_info$Cohort=="Replication"]]
title_name <- "Replication"
group_info <- data.frame(sample=colnames(pdat), group=ifelse(grepl("Control", colnames(pdat)), "Control","LM"))
#pls-da
pdat <- t(pdat)
pls_da <- plsda(pdat, group_info$group, ncomp = 5)
res <- as.data.frame(pls_da$variates$X)
res$group <- group_info$group
res$sample <- group_info$sample
res <- merge(res, pheno, by.x = "sample",by.y = "Patient", all.x = T, all.y = F)

result <- c()
for(i in 1:5){
	value <- res[, paste0("comp", i)] %>% as.numeric()
	value_control <- value[which(res$group=="Control")]
	value_case <- value[which(res$group=="LM")]
	p <- wilcox.test(value_control, value_case)
	p_group <- p$p.value
	#
	value_control <- value[which(res$Sex=="Female")]
	value_case <- value[which(res$Sex=="Male")]
	p <- wilcox.test(value_control, value_case)
	p_sex <- p$p.value
	#
	value_x <- res[, paste0("comp", i)] %>% as.numeric()
	value_y <- res[, "Age"] %>% as.numeric()
	p <- cor.test(value_x, value_y, method = c("pearson"))
	p_age <- p$p.value
	tmp <- c(p_group, p_sex, p_age)
	result <- rbind(result, tmp)
}
colnames(result) <- c("condition","sex","age")
rownames(result) <- paste0("comp",1:5)
summ <- summary(pca)
result <- as.data.frame(result)
result$contribution <- paste0(round(pls_da$prop_expl_var$X[1:5]*100,2),"%")
result_rep_plsda <- result
#
write.table(result_dis_pca, file = "result_dis_pca.txt", row.names = T, col.names =T, sep = "\t")
write.table(result_dis_plsda, file = "result_dis_plsda.txt", row.names = T, col.names =T, sep = "\t")
write.table(result_rep_pca, file = "result_rep_pca.txt", row.names = T, col.names =T, sep = "\t")
write.table(result_rep_plsda, file = "result_rep_plsda.txt", row.names = T, col.names =T, sep = "\t")
