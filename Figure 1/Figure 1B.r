library(dplyr)
library(ComplexHeatmap)
library(circlize)
library(grid)
library(openxlsx)
#
setwd("~/LM/Figure1")
data <- read.xlsx("~/LM/data/patients_clinical_information.xlsx", startRow = 2)
data[is.na(data)] <- "N"
NA_COL <- "#EDEDED"
data2 <- data %>%
  mutate(
    Cohort    = factor(Cohort, levels = c("Discovery", "Replication")),
    Group = factor(Group, levels = c("LM", "Control")),
	Response  = factor(Response, levels = c("PR","SD","PD","N"))
  ) %>%
  arrange(Cohort, Group, Source, Cancer)
sample_order <- data2$Patient
data2 <- data2 %>%
  mutate(
    iPFS_time = as.numeric(ifelse(iPFS_time == "N", NA, iPFS_time)),
    OS_time   = as.numeric(ifelse(OS_time   == "N", NA, OS_time))
  )
clinical_vars <- c(
  "Cohort", "Group", "Source", 
  "Age", "Sex", "Cancer", "Subtype",
  "iPFS_status", "iPFS_time",
  "OS_status", "OS_time",
  "Response", "MRI", "Cytology"
)

anno_df <- data2 %>%
  select(all_of(clinical_vars))
col_fun_age <- colorRamp2(
  c(min(data2$Age, na.rm = TRUE), max(data2$Age, na.rm = TRUE)),
  c("#4575B4", "#D73027")
)

col_fun_iPFS <- colorRamp2(
  c(0, max(data2$iPFS_time, na.rm = TRUE)),
  c("white", "#762A83")
)

col_fun_OS <- colorRamp2(
  c(0, max(data2$OS_time, na.rm = TRUE)),
  c("white", "#A50F15")
)
col_list <- list(

  Cohort = c(
    Discovery   = "#4DBBD5",
    Replication = "#E64B35"
  ),

  Group = c(
    LM = "#00A087",
    Control   = "#DC0000"
  ),

  Sex = c(
    Male   = "#3C5488",
    Female = "#E64B35"
  ),

  Response = c(
    PR = "#00A087",
    SD = "#F39B7F",
    PD = "#DC0000",
    N  = NA_COL
  ),

  MRI = c(
    Positive = "#DC0000",
	Negative = "#4DBBD5",
    N = NA_COL
  ),

  Cytology = c(
    Positive = "#DC0000",
	Negative = "#4DBBD5",
    N = NA_COL
  ),

  Cancer = c(
    BC = "#3C5488",
    LA = "#E64B35",
    EC = "#00A087",
	MM = "#F39B7F",
    N  = NA_COL
  ),

  Subtype = c(
    Adenocarcinoma = "#4DBBD5",
    TNBC           = "#DC0000",
    `HER2+`        = "#00A087",
	Other          = "#7E6148",  
    `Luminal A/B`  = "#F39B7F",  
    N              = NA_COL
  ),

  iPFS_status = c(
    "0" = "#4DBBD5",
    "1" = "#DC0000",
    "N" = NA_COL
  ),

  OS_status = c(
    "0" = "#4DBBD5",
    "1" = "#DC0000",
    "N" = NA_COL
  )
)

rownames(anno_df) <- data2$Patient

ha <- HeatmapAnnotation(
  df = anno_df,
  col = c(
    list(
      Age       = col_fun_age,
      iPFS_time = col_fun_iPFS,
      OS_time   = col_fun_OS
    ),
    col_list
  ),
  na_col = NA_COL,
  annotation_name_side = "left",
  gp = gpar(col = "gray", lwd = 0.5)  
)

cohort_split <- anno_df$Cohort   
pdf("sample_info.pdf", width = 12)
p <- Heatmap(
  matrix(0, nrow = 1, ncol = nrow(anno_df)),
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = FALSE,

  column_split = cohort_split,         
  column_gap = unit(3, "mm"),     

  top_annotation = ha,
  height = unit(3, "mm"),
  rect_gp = gpar(col = NA)
)

print(p)
dev.off()