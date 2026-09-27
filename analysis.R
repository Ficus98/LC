##########################Figure1C
library(ggplot2)
library(reshape2)
library(RColorBrewer)

# 假设 mat 是原始数据，症状列在 13:ncol(mat)
mat <- read.csv("heatmap.csv", stringsAsFactors = FALSE)
symptom_cols <- 13:ncol(mat)
symptom_names <- colnames(mat)[13:ncol(mat)]

# 生成统一颜色，21种症状
colors <- colorRampPalette(brewer.pal(12, "Set3"))(length(symptom_names))

# 定义分组列表
group_list <- list(
  Age = list(colname = "agegroup", levels = c(2,3,4), labels = c("31–40", "41–50", "≥51")),
  CD4 = list(colname = "CD4group", levels = c(1,2,3), labels = c("<350","350–499","≥500")),
  VL = list(colname = "VLgroup", levels = c(1,2), labels = c("<50", "≥50")),
  Treatment = list(colname = "treatgroup", levels = c(1,3,4), labels = c("2NRTIs+NNRTI","2NRTIs+INSTI","3TC+DTG")),
  Vaccination = list(colname = "vaccgroup", levels = c(1,2), labels = c("Unvaccinated","Vaccinated")),
  Severity = list(colname = "severegroup", levels = c(1,2), labels = c("Non-severe","Severe"))
)

# 初始化空的dataframe
plot_df <- data.frame()

# 循环处理每个分组
for (grp_name in names(group_list)) {
  grp_info <- group_list[[grp_name]]
  col <- grp_info$colname
  mat[, col] <- factor(mat[, col], levels = grp_info$levels, labels = grp_info$labels)
  
  # 转成长格式
  temp_long <- melt(mat[, c(col, symptom_names)],
                    id.vars = col,
                    variable.name = "Symptom",
                    value.name = "Present")
  
  temp_long$GroupType <- grp_name  # 标记是哪类分组
  colnames(temp_long)[1] <- "Group" # 统一列名
  plot_df <- rbind(plot_df, temp_long)
}

# 绘图：横坐标为组别+分组类型，堆叠显示症状
ph <- ggplot(plot_df, aes(x = Group, y = Present, fill = Symptom)) +
  geom_bar(stat = "identity") +
  facet_wrap(~GroupType, scales = "free_x", nrow = 2) +
  scale_fill_manual(values = colors) +
  theme_bw(base_family = "Arial") + # 整体字体为 Arial
  ylab("Number of subjects with symptom") +
  xlab("Group") +
  theme(
    axis.text.x = element_text(angle = 40, hjust = 1, color="black",face = "bold",size = 10),
    axis.text.y = element_text(face = "bold"),
    axis.title.x = element_text(face = "bold",size = 12),
    axis.title.y = element_text(face = "bold", margin = margin(r = 10),size = 12),
    strip.text = element_text(face = "bold",size = 10), # facet 标签
    strip.background = element_blank(),
    legend.text = element_text(face = "bold",size = 11),
    legend.title = element_text(face = "bold"),
    legend.position = "right",
    legend.key.size = unit(0.5, "cm"),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  ) +
  guides(fill = guide_legend(ncol = 1))
ph
ggsave("症状热图600.tiff", plot = ph, width = 6.5, height = 5.5,dpi=600)
##############################################################Figure2
###########################火山图 500X400
# 加载必要包
library(dplyr)
library(ggplot2)
library(ggrepel)

# 读取数据
df <- read.csv("volcano.csv")
group <- df[[2]]  # 第二列为分组 LC（0/1）
expr <- df[, 3:45]  # 第3至第49列为细胞因子

# 差异分析 + log2FC 计算
res <- lapply(colnames(expr), function(gene) {
  group0 <- expr[group == 0, gene]
  group1 <- expr[group == 1, gene]
  test <- wilcox.test(group1, group0)
  log2fc <- log2(median(group1, na.rm = TRUE) / median(group0, na.rm = TRUE))
  data.frame(
    marker = gene,
    log2FC = log2fc,
    pvalue = test$p.value
  )
}) %>% bind_rows()

# 根据 p 值和 log2FC 标记显著性
res$threshold <- with(res, ifelse(pvalue < 0.05 & abs(log2FC) > 0.25,
                                  ifelse(log2FC > 0.25, "Up", "Down"), "NS"))

# 火山图绘制（基于 p 值）
ggplot(res, aes(x = log2FC, y = -log10(pvalue), color = threshold, label = marker)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(data = subset(res, threshold != "NS"), size = 4, max.overlaps = 15) +
  scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey60")) +
  geom_vline(xintercept =0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlim(-1.5, 1.5) +  # 限制横坐标范围
  theme_bw(base_size = 14) +
  labs(title = NULL,
       x = "log2(Fold Change)",
       y = "-log10(p-value)",
       color = "")
############################################### PCA 500X400
library(ggplot2)
library(dplyr)

# 读数据
df <- read.csv("volcano.csv")

# 设置分组
group <- factor(df[[2]], levels = c(0, 1), labels = c("NLC", "LC"))

# 选择表达矩阵并处理NA
selected_markers <- c("CCL27","IL1A", "IL2", "IL3", "CXCL8", "IL17A","IFNA2")
expr <- df[, selected_markers]

#expr <- df[, 3:45]
expr[is.na(expr)] <- 0

# 标准化
expr_scaled <- scale(expr)

# PCA
pca_res <- prcomp(expr_scaled, center = TRUE, scale. = TRUE)
explained_var <- summary(pca_res)$importance[2, 1:2] * 100

# 构建数据框
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  group = group
)

# 绘图：带填充的置信椭圆
ggplot(pca_df, aes(x = PC1, y = PC2, color = group, fill = group)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +  # 填充椭圆
  geom_point(size = 2, alpha = 1) +
  scale_color_manual(values = c("NLC" = "#4575B4", "LC" = "#F8766D")) + # 自定义点边框颜色
  scale_fill_manual(values = c("NLC" = "#4575B4", "LC" = "#F8766D")) +   # 自定义椭圆填充颜色
  theme_bw(base_size = 14) +
  labs(
    title = NULL,
    x = paste0("PC1 (", round(explained_var[1], 1), "% variance)"),
    y = paste0("PC2 (", round(explained_var[2], 1), "% variance)"),
    color = "Group",
    fill = "Group"
  )
################################################热图
library(ComplexHeatmap)
library(circlize)
library(dplyr)

# 读入数据
df <- read.csv("volcano.csv")

# 提取分组信息
group <- factor(df[[2]], levels = c(0, 1), labels = c("NLC", "LC"))
sample_names <- df[[1]]

# 提取表达矩阵
selected_markers <- c("CCL27", "IL1A", "IL2", "IL3", "CXCL8", "IL17A","IFNA2")
expr <- df[, selected_markers]
rownames(expr) <- sample_names
expr[is.na(expr)] <- 0

# 按照分组排序样本
df$GroupLabel <- group
df_sorted <- df[order(df$GroupLabel), ]
group_sorted <- df_sorted$GroupLabel
sample_sorted <- df_sorted[[1]]
expr_sorted <- df_sorted[, selected_markers]
rownames(expr_sorted) <- sample_sorted

# 转置表达矩阵（行为marker，列为样本）
expr_t <- t(expr_sorted)

# 标准化（按marker标准化）
expr_scaled <- t(scale(t(expr_t)))

# 构建上方注释
ha_col <- HeatmapAnnotation(
  Group = group_sorted,
  show_annotation_name = FALSE,
  col = list(Group = c("LC" = "#81C784", "NLC" = "#CE93D8"))
)

# 绘制热图
Heatmap(
  expr_scaled,
  name = "Z-score",
  col = colorRamp2(c(-3, 0, 3), c("#2166ACFF", "white", "#B2182BFF")),
  top_annotation = ha_col,
  cluster_columns = FALSE,   # ❗不对样本聚类
  cluster_rows = TRUE,       # 对marker聚类
  show_column_names = FALSE,
  show_row_names = TRUE,
  column_title = NULL,
  row_title = NULL,
  column_names_rot = 45,
  heatmap_legend_param = list(at = c(-3, 0, 3), labels = c("-3", "0", "3"))
)

#################################富集分析
# 1. 加载必要的包
# install.packages("BiocManager")
# BiocManager::install(c("clusterProfiler", "org.Hs.eg.db", "enrichplot"))

library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 2. 准备你的差异基因列表（Symbol）
gene_list <- c("IL1A", "IL2", "IL3", "CXCL8", "IL17A", "IFNA2")

# 3. ID 转换：将基因 Symbol 转换为 Entrez Gene ID
# org.Hs.eg.db 是人类基因注释数据库
# 从 'SYMBOL' 列映射到 'ENTREZID' 列
entrez_ids <- bitr(geneID = gene_list,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# 4. 进行 GO 和 KEGG 富集分析

# GO - 生物过程 (Biological Process, BP)
ego_bp <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

# GO - 细胞组分 (Cellular Component, CC)
ego_cc <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 1
)

# GO - 分子功能 (Molecular Function, MF)
ego_mf <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05,
                   qvalueCutoff  = 0.05)

# KEGG 富集分析
ekegg <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                    organism     = 'hsa', # 人类 KEGG 库
                    pvalueCutoff = 0.05)
#options(timeout = 300)
# 5. 导出所有 GO 和 KEGG 富集结果为 CSV 文件

# 导出 GO-BP 结果
if (!is.null(ego_bp)) {
  ego_bp_df <- as.data.frame(ego_bp)
  write.csv(ego_bp_df, "GO_BP_enrichment_results.csv", row.names = FALSE)
  print("GO生物过程富集结果已导出至 'GO_BP_enrichment_results.csv'")
} else {
  print("没有显著的GO生物过程富集通路，未导出文件。")
}

# 导出 GO-CC 结果
if (!is.null(ego_cc)) {
  ego_cc_df <- as.data.frame(ego_cc)
  write.csv(ego_cc_df, "GO_CC_enrichment_results.csv", row.names = FALSE)
  print("GO细胞组分富集结果已导出至 'GO_CC_enrichment_results.csv'")
} else {
  print("没有显著的GO细胞组分富集通路，未导出文件。")
}

# 导出 GO-MF 结果
if (!is.null(ego_mf)) {
  ego_mf_df <- as.data.frame(ego_mf)
  write.csv(ego_mf_df, "GO_MF_enrichment_results.csv", row.names = FALSE)
  print("GO分子功能富集结果已导出至 'GO_MF_enrichment_results.csv'")
} else {
  print("没有显著的GO分子功能富集通路，未导出文件。")
}

# 导出 KEGG 结果
if (!is.null(ekegg)) {
  ekegg_df <- as.data.frame(ekegg)
  write.csv(ekegg_df, "KEGG_enrichment_results.csv", row.names = FALSE)
  print("KEGG富集结果已导出至 'KEGG_enrichment_results.csv'")
} else {
  print("没有显著的KEGG富集通路，未导出文件。")
}
############################################气泡GO
#########################################################Go
library(ggplot2)
library(dplyr)
library(forcats)

# 读取数据
go_df <- read.csv("GO1.csv", stringsAsFactors = FALSE)

# 按显著性排序 GO 条目
go_df <- go_df %>%
  group_by(Category) %>%
  arrange(p.adjust, .by_group = TRUE) %>%
  ungroup()

# 将 GO 条目按照 -log10(pvalue) 排序
go_df <- go_df %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# 定义颜色
colors <- c(
  "Biological Process" = "#1f77b4",
  "Molecular Function" = "#ff7f0e",
  "Cellular Component" = "#2ca02c"
)

# 创建纵坐标对应颜色向量
y_colors <- colors[go_df$Category]
names(y_colors) <- levels(go_df$Description)  # 对应 Description

# 绘图
p <- ggplot(go_df, aes(x = -log10(p.adjust), y = Description)) +
  geom_point(aes(size = Count, fill = Category), shape = 21, color = "black") +
  scale_fill_manual(values = colors) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw() +   # 四边框
  theme(
    axis.text.y = element_text(face = "bold", size = 11, color = y_colors),
    axis.text.x = element_text(face = "bold", size = 11, color = "black"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 10, face = "bold"),
    legend.text = element_text(size = 10, face = "bold")
  ) +
  ylab("") +
  xlab("-Log10 adjust.p") +
  guides(
    fill = guide_legend(title = "Category", override.aes = list(size = 5))
  )

# 保存
ggsave("GO_triplet_bubbleplot_colored_y.tiff", plot = p, width = 10.5, height = 8)
########################################################KEGG
library(ggplot2)
library(readr)
library(dplyr)
library(forcats)
library(ggtext)  # 用于Markdown文本渲染
#install.packages("ggtext")
# 读取数据
df <- read_csv("kegg.csv")

# 设定颜色
category_colors <- c(
  "Metabolism" = "#3C5488",
  "Genetic Information Processing" = "#8491B4",
  "Environmental Information Processing" = "#D64B39",
  "Cellular Processes" = "#00A087",
  "Organismal Systems" = "#E69F00",
  "Human Diseases" = "#16B4E9"
)

# 排序并重排因子
df <- df %>%
  arrange(Category, desc(adjust.p)) %>%
  mutate(Description = fct_reorder(Description, adjust.p, .desc = TRUE),
         # 给每条通路名称加上颜色的HTML标签
         Description_colored = paste0("<span style='color:", category_colors[Category], "'>", Description, "</span>")
  )

# 缩放系数
scaleFactor <- max(-log10(df$adjust.p)) / max(df$Count)

ggplot(df, aes(x = Description, y = -log10(adjust.p), fill = Category)) +
  geom_col() +
  geom_point(aes(y = Count * scaleFactor), color = "black", size = 2) +
  geom_line(aes(y = Count * scaleFactor, group = 1), color = "black") +
  scale_y_continuous(
    name = expression(-log[10](adjust.p)),
    sec.axis = sec_axis(~ . / scaleFactor, name = "Gene Count")
  ) +
  coord_flip() +
  scale_fill_manual(values = category_colors) +
  theme_bw() +
  theme(
    axis.text.y = element_markdown(size = 11, face = "bold"),  # 使用Markdown渲染带颜色文本
    axis.title.y = element_blank(),
    legend.position = c(0.51, 0.45),  # 图内位置，x=0.05靠左，y=0.85靠上
    legend.justification = c(0, 1),   # 图例左上角对齐
    legend.background = element_rect(fill = alpha("white", 0.7), color = NA), # 半透明背景
    legend.title = element_text(size = 8),
    legend.text = element_text(size = 8)
  ) +
  aes(x = Description_colored)+  # 替换为带颜色的纵坐标
  guides(fill = guide_legend(nrow = 5, byrow = TRUE))
##########################一致性聚类pam+pearson
library(ConsensusClusterPlus)

# 输入文件路径和输出目录，请自行修改
inputFile <- "LC90_scaled.csv"
outputDir <- "consensus_results"

# 创建输出目录（如果不存在）
if (!dir.exists(outputDir)) dir.create(outputDir)

# 读取临床数据（第一列为样本名，已设为行名）
clinicalData <- read.csv(inputFile, header = TRUE, row.names = 1, na.strings = "NA")

# 查看数据结构（可选）
str(clinicalData)

# 删除有缺失值的样本（行）
completeData <- na.omit(clinicalData)

# 转换为数值矩阵
#dataMatrix <- as.matrix(completeData)
dataMatrix <- t(as.matrix(completeData))
# 设置聚类参数
maxK <- 6      # 最大聚类数
reps <- 50      # 重抽样次数

# 执行共识聚类
results <- ConsensusClusterPlus(
  d = dataMatrix,
  maxK = maxK,
  reps = reps,
  pItem = 0.8,           # 每次抽80%样本
  pFeature = 1,          # 每次用100%特征
  clusterAlg = "pam",    # PAM算法
  distance = "pearson",# 距离：欧氏距离euclidean
  seed = 1234,
  title = outputDir,
  plot = "png"           # 输出png格式图形
)

# 计算并绘制聚类一致性指标（可选）
icl <- calcICL(results, title = outputDir, plot = "png")

# 根据结果选择最佳聚类数，比如这里选2
clusterNum <- 2

# 提取最佳聚类的分类结果
cluster <- results[[clusterNum]][["consensusClass"]]
cluster_df <- data.frame(Cluster = cluster)

# 给每个聚类编号分配字母标签，方便后续展示
letters_vec <- c("A","B","C","D","E","F")
uniqClu <- sort(unique(cluster_df$Cluster))
cluster_df$Cluster <- letters_vec[match(cluster_df$Cluster, uniqClu)]

# 合并聚类结果和完整数据
finalData <- cbind(completeData, Cluster = cluster_df$Cluster)

# 保存结果到CSV文件
write.csv(finalData, file = file.path(outputDir, "clinical_clustering_results.csv"), row.names = TRUE)
write.csv(cluster_df, file = file.path(outputDir, "cluster_assignments.csv"), row.names = TRUE)
#################################################火山图A
# 加载必要包
library(dplyr)
library(ggplot2)
library(ggrepel)

# 读取数据
df <- read.csv("LC84AAb90-clusterA.csv")
group <- df[[4]]  # 第二列为分组 LC（0/1）
expr <- df[, 5:47]  # 第3至第49列为细胞因子

# 差异分析 + log2FC 计算
res <- lapply(colnames(expr), function(gene) {
  group0 <- expr[group == 0, gene]
  group1 <- expr[group == 1, gene]
  test <- wilcox.test(group1, group0)
  log2fc <- log2(median(group1, na.rm = TRUE) / median(group0, na.rm = TRUE))
  data.frame(
    marker = gene,
    log2FC = log2fc,
    pvalue = test$p.value
  )
}) %>% bind_rows()

# 根据 p 值和 log2FC 标记显著性
res$threshold <- with(res, ifelse(pvalue < 0.05 & abs(log2FC) > 0.15,
                                  ifelse(log2FC > 0.15, "Up", "Down"), "NS"))

# 火山图绘制（基于 p 值）
pv <- ggplot(res, aes(x = log2FC, y = -log10(pvalue), color = threshold, label = marker)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text_repel(data = subset(res, threshold != "NS"), size = 4.5, fontface = "bold",max.overlaps = 15) +
  scale_color_manual(values = c("Up" = "#fF4F39", "Down" = "#2CA02C", "NS" = "grey60")) +
  geom_vline(xintercept =0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlim(-1.5, 1.5) +  # 限制横坐标范围
  theme_bw(base_size = 14) +
  labs(title = NULL,
       x = "log2(Fold Change)",
       y = "-log10(p-value)",
       color = "")
ggsave("vocanol_plotA.tiff", plot = pv, width = 5.5, height = 4.5, dpi = 600)

########################################LASSO
# -----------------------------
# Repeated 5-fold LASSO (no nested CV)
# -----------------------------
library(data.table)
library(caret)
library(glmnet)
library(dplyr)
library(ggplot2)

# 读取数据
data <- fread("LC84AAb90-clusterA.csv")

# 特征与标签
feature_names <- c(
  "CXCL10", "CCL5","IL2RA","IL1RN",
  "CXCL9","TNFSF10","IL18","CXCL8","CD4CD8"
)

X <- as.matrix(data[, ..feature_names])
y <- as.numeric(data[[4]])  # 0 / 1

# -----------------------------
# 参数设置
# -----------------------------
set.seed(1111)
n_repeat <- 50
k <- 5

coef_list <- list()
cnt <- 1

# -----------------------------
# Repeated 5-fold CV
# -----------------------------
for (r in 1:n_repeat) {
  
  folds <- createFolds(y, k = k)
  
  for (i in seq_along(folds)) {
    
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)
    
    # 仅在训练集标准化（防止数据泄露）
    x_train <- scale(X[train_idx, ])
    y_train <- y[train_idx]
    
    cvfit <- cv.glmnet(
      x = x_train,
      y = y_train,
      family = "binomial",
      alpha = 1,
      type.measure = "deviance"
    )
    
    coef_vec <- as.numeric(coef(cvfit, s = "lambda.min"))[-1]
    names(coef_vec) <- feature_names
    
    coef_list[[cnt]] <- data.frame(
      Feature = feature_names,
      Coef = coef_vec
    )
    cnt <- cnt + 1
  }
}

# -----------------------------
# 汇总结果
# -----------------------------
coef_df <- bind_rows(coef_list) %>%
  mutate(
    Selected = Coef != 0,
    AbsCoef = abs(Coef)
  ) %>%
  group_by(Feature) %>%
  summarise(
    SelectionFreq = mean(Selected),     # 被选中的比例
    MeanCoef = mean(Coef),
    MeanAbsCoef = mean(AbsCoef),
    .groups = "drop"
  ) %>%
  arrange(desc(SelectionFreq), desc(MeanAbsCoef))

print(coef_df)

p_lasso <- ggplot(
  coef_df,
  aes(
    x = reorder(Feature, MeanAbsCoef),
    y = MeanAbsCoef,
    fill = SelectionFreq        # 👈 核心：颜色 = 选择频率
  )
) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  scale_fill_gradient(
    low = "#D9D9D9",
    high = "#1ABAB3",
    name = "Selection\nfrequency"
  ) +
  theme_minimal(base_size = 13) +
  labs(
    x = "Feature",
    y = "Mean |Coefficient|",
    title = "Bootstrap-based LASSO Variable Importance"
  ) +
  theme(
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(color = "black", face = "bold"),
    plot.title = element_text(size = 10, face = "bold"),
    axis.title.x = element_text(color = "black", face = "bold"),
    axis.title.y = element_text(color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      linewidth = 0.5
    ),
    legend.title = element_text(face = "bold",size = 8),
    legend.text  = element_text(size = 8)
  )
ggsave("LASSO_importanceA.tiff",
       plot = p_lasso, width = 4.3, height = 2.8, dpi = 600)
library(pROC)
ranked_features <- feature_names  # 如果 coef_df 还没生成，可以用全部

# -----------------------------
# 2. 定义 incremental ROC 函数
# -----------------------------
library(pROC)

calc_auc_boot <- function(X, y, features, n_boot = 2000) {
  
  df <- data.frame(y = y, X[, features, drop = FALSE])
  
  fit <- glm(y ~ ., data = df, family = binomial)
  prob <- predict(fit, type = "response")
  
  roc_obj <- roc(y, prob, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  # bootstrap CI
  ci_auc <- ci.auc(roc_obj, method = "bootstrap", boot.n = n_boot)
  
  data.frame(
    AUC = auc_val,
    LowerCI = ci_auc[1],
    UpperCI = ci_auc[3]
  )
}
auc_results <- list()

for(k in seq_along(ranked_features)) {
  
  current_features <- ranked_features[1:k]
  
  res <- calc_auc_boot(
    X = X,
    y = y,
    features = current_features,
    n_boot = 2000
  )
  
  auc_results[[k]] <- data.frame(
    Model = paste0("Top ", k),
    Features = paste(current_features, collapse = ", "),
    MeanAUC = res$AUC,
    LowerCI = res$LowerCI,
    UpperCI = res$UpperCI
  )
}

auc_df <- bind_rows(auc_results)

# -----------------------------
# 4. 绘制森林图（带 AUC ± CI 标签）
# -----------------------------
p_auc <- ggplot(auc_df,
                aes(y = Model, x = MeanAUC)) +
  geom_point(size = 3, color = "#1F77B4") +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.2, color = "#1F77B4") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_text(
    aes(
      x = UpperCI + 0.03,
      label = sprintf("%.2f (%.2f–%.2f)", MeanAUC, LowerCI, UpperCI)
    ),
    fontface = "bold",
    size = 3.8,
    color = "black",
    hjust = 0
  )+  # 标签在点右侧
  scale_x_continuous(
    limits = c(0.4, 1.1),
    expand = c(0, 0)
  )+  # 给标签留空间
  theme_minimal(base_size = 13) +
  labs(
    x = "AUC",
    y = "Cumulative feature set",
    title = "Performance Based on LASSO Ranking"
  )
p_auc <- p_auc +
  theme(
    # 坐标轴文字
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(color = "black", face = "bold"),
    
    # 坐标轴标题
    axis.title.x = element_text(color = "black", face = "bold"),
    axis.title.y = element_text(color = "black", face = "bold"),
    
    # 外框（重点）
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    # 去掉网格线（森林图更干净）
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("AUC_forestA.tiff", plot = p_auc, width = 6.5, height = 3, dpi = 600)
######################################################PCA
library(ggplot2)
library(dplyr)

# 读数据
df <- read.csv("LC84AAb90-clusterA.csv")

# 设置分组
group <- factor(df[[4]], levels = c(0, 1), labels = c("NLC", "LC"))

# 选择表达矩阵并处理NA
selected_markers <- c(
  "CXCL10","CD4CD8"
)#, "CCL5","IL2RA","IL1RN","CXCL9","TNFSF10","IL18","CXCL8"
expr <- df[, selected_markers] 

#expr <- df[, 3:45]
expr[is.na(expr)] <- 0

# 标准化
expr_scaled <- scale(expr)

# PCA
pca_res <- prcomp(expr_scaled, center = TRUE, scale. = TRUE)
explained_var <- summary(pca_res)$importance[2, 1:2] * 100

# 构建数据框
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  group = group
)

# 绘图：带填充的置信椭圆
pp <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, fill = group)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +  # 填充椭圆
  geom_point(size = 2, alpha = 1) +
  scale_color_manual(values = c("NLC" = "#fF4F39","LC" = "#2CA02C" )) + # 自定义点边框颜色
  scale_fill_manual(values = c("NLC" = "#fF4F39","LC" = "#2CA02C")) +   # 自定义椭圆填充颜色
  theme_bw(base_size = 14) +
  labs(
    title = NULL,
    x = paste0("PC1 (", round(explained_var[1], 1), "% variance)"),
    y = paste0("PC2 (", round(explained_var[2], 1), "% variance)"),
    color = "Group",
    fill = "Group"
  )
ggsave("pca_plotA.tiff", plot = pp, width = 4, height = 3, dpi = 600)

#################################################火山图B
# 加载必要包
library(dplyr)
library(ggplot2)
library(ggrepel)

# 读取数据
df <- read.csv("LC84AAb90-clusterB.csv")
group <- df[[4]]  # 第二列为分组 LC（0/1）
expr <- df[, 5:47]  # 第3至第49列为细胞因子

# 差异分析 + log2FC 计算
res <- lapply(colnames(expr), function(gene) {
  group0 <- expr[group == 0, gene]
  group1 <- expr[group == 1, gene]
  test <- wilcox.test(group1, group0)
  log2fc <- log2(median(group1, na.rm = TRUE) / median(group0, na.rm = TRUE))
  data.frame(
    marker = gene,
    log2FC = log2fc,
    pvalue = test$p.value
  )
}) %>% bind_rows()

# 根据 p 值和 log2FC 标记显著性
res$threshold <- with(res, ifelse(pvalue < 0.05 & abs(log2FC) > 0.15,
                                  ifelse(log2FC > 0.15, "Up", "Down"), "NS"))

# 火山图绘制（基于 p 值）
pv <- ggplot(res, aes(x = log2FC, y = -log10(pvalue), color = threshold, label = marker)) +
  geom_point(size = 4, alpha = 0.7) +
  geom_text_repel(data = subset(res, threshold != "NS"), size = 4.5, fontface = "bold",max.overlaps = 15) +
  scale_color_manual(values = c("Up" = "#fF4F39", "Down" = "#2CA02C", "NS" = "grey60")) +
  geom_vline(xintercept =0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlim(-1.5, 1.5) +  # 限制横坐标范围
  theme_bw(base_size = 14) +
  labs(title = NULL,
       x = "log2(Fold Change)",
       y = "-log10(p-value)",
       color = "")
ggsave("vocanol_plotB.tiff", plot = pv, width = 5.5, height = 4.5, dpi = 600)
########################################LASSO
# -----------------------------
# Repeated 5-fold LASSO (no nested CV)
# -----------------------------
library(data.table)
library(caret)
library(glmnet)
library(dplyr)
library(ggplot2)

# 读取数据
data <- fread("LC84AAb90-clusterB.csv")

# 特征与标签
feature_names <- c("CCL27","IL3","IL2","IL1A","IL17A","Treattype","VD")

X <- as.matrix(data[, ..feature_names])
y <- as.numeric(data[[4]])  # 0 / 1

# -----------------------------
# 参数设置
# -----------------------------
set.seed(1111)
n_repeat <- 50
k <- 5

coef_list <- list()
cnt <- 1

# -----------------------------
# Repeated 5-fold CV
# -----------------------------
for (r in 1:n_repeat) {
  
  folds <- createFolds(y, k = k)
  
  for (i in seq_along(folds)) {
    
    test_idx  <- folds[[i]]
    train_idx <- setdiff(seq_len(nrow(X)), test_idx)
    
    # 仅在训练集标准化（防止数据泄露）
    x_train <- scale(X[train_idx, ])
    y_train <- y[train_idx]
    
    cvfit <- cv.glmnet(
      x = x_train,
      y = y_train,
      family = "binomial",
      alpha = 1,
      type.measure = "deviance"
    )
    
    coef_vec <- as.numeric(coef(cvfit, s = "lambda.min"))[-1]
    names(coef_vec) <- feature_names
    
    coef_list[[cnt]] <- data.frame(
      Feature = feature_names,
      Coef = coef_vec
    )
    cnt <- cnt + 1
  }
}

# -----------------------------
# 汇总结果
# -----------------------------
coef_df <- bind_rows(coef_list) %>%
  mutate(
    Selected = Coef != 0,
    AbsCoef = abs(Coef)
  ) %>%
  group_by(Feature) %>%
  summarise(
    SelectionFreq = mean(Selected),     # 被选中的比例
    MeanCoef = mean(Coef),
    MeanAbsCoef = mean(AbsCoef),
    .groups = "drop"
  ) %>%
  arrange(desc(SelectionFreq), desc(MeanAbsCoef))

print(coef_df)

p_lasso <- ggplot(
  coef_df,
  aes(
    x = reorder(Feature, MeanAbsCoef),
    y = MeanAbsCoef,
    fill = SelectionFreq        # 👈 核心：颜色 = 选择频率
  )
) +
  geom_bar(stat = "identity", width = 0.75) +
  coord_flip() +
  scale_fill_gradient(
    low = "#D9D9D9",
    high = "#1ABAB3",
    name = "Selection\nfrequency"
  ) +
  theme_minimal(base_size = 13) +
  labs(
    x = "Feature",
    y = "Mean |Coefficient|",
    title = "Bootstrap-based LASSO Variable Importance"
  ) +
  theme(
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(color = "black", face = "bold"),
    plot.title = element_text(size = 8.5, face = "bold"),
    axis.title.x = element_text(color = "black", face = "bold"),
    axis.title.y = element_text(color = "black", face = "bold"),
    panel.grid = element_blank(),
    panel.border = element_rect(
      color = "black",
      fill  = NA,
      linewidth = 0.5
    ),
    legend.title = element_text(face = "bold",size = 8),
    legend.text  = element_text(size = 8)
  )
ggsave("LASSO_importanceB.tiff",
       plot = p_lasso, width = 3.6, height = 2.5, dpi = 600)
library(pROC)
ranked_features <- feature_names  # 如果 coef_df 还没生成，可以用全部

# -----------------------------
# 2. 定义 incremental ROC 函数
# -----------------------------
library(pROC)

calc_auc_boot <- function(X, y, features, n_boot = 2000) {
  
  df <- data.frame(y = y, X[, features, drop = FALSE])
  
  fit <- glm(y ~ ., data = df, family = binomial)
  prob <- predict(fit, type = "response")
  
  roc_obj <- roc(y, prob, quiet = TRUE)
  auc_val <- as.numeric(auc(roc_obj))
  
  # bootstrap CI
  ci_auc <- ci.auc(roc_obj, method = "bootstrap", boot.n = n_boot)
  
  data.frame(
    AUC = auc_val,
    LowerCI = ci_auc[1],
    UpperCI = ci_auc[3]
  )
}
auc_results <- list()

for(k in seq_along(ranked_features)) {
  
  current_features <- ranked_features[1:k]
  
  res <- calc_auc_boot(
    X = X,
    y = y,
    features = current_features,
    n_boot = 2000
  )
  
  auc_results[[k]] <- data.frame(
    Model = paste0("Top ", k),
    Features = paste(current_features, collapse = ", "),
    MeanAUC = res$AUC,
    LowerCI = res$LowerCI,
    UpperCI = res$UpperCI
  )
}

auc_df <- bind_rows(auc_results)

# -----------------------------
# 4. 绘制森林图（带 AUC ± CI 标签）
# -----------------------------
p_auc <- ggplot(auc_df,
                aes(y = Model, x = MeanAUC)) +
  geom_point(size = 3, color = "#1F77B4") +
  geom_errorbarh(aes(xmin = LowerCI, xmax = UpperCI), height = 0.2, color = "#1F77B4") +
  geom_vline(xintercept = 0.5, linetype = "dashed", color = "grey50") +
  geom_text(
    aes(
      x = LowerCI - 0.15,
      label = sprintf("%.2f (%.2f–%.2f)", MeanAUC, LowerCI, UpperCI)
    ),
    fontface = "bold",
    size = 3.8,
    color = "black",
    hjust = 0
  )+  # 标签在点右侧
  scale_x_continuous(
    limits = c(0.4, 1.10),
    expand = c(0, 0)
  )+  # 给标签留空间
  theme_minimal(base_size = 13) +
  labs(
    x = "AUC",
    y = "Cumulative feature set",
    title = "Performance Based on LASSO Ranking"
  )
p_auc <- p_auc +
  theme(
    # 坐标轴文字
    axis.text.y = element_text(color = "black", face = "bold"),
    axis.text.x = element_text(color = "black", face = "bold"),
    
    # 坐标轴标题
    axis.title.x = element_text(color = "black", face = "bold"),
    axis.title.y = element_text(color = "black", face = "bold"),
    
    # 外框（重点）
    panel.border = element_rect(
      color = "black",
      fill = NA,
      linewidth = 0.8
    ),
    
    # 去掉网格线（森林图更干净）
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank()
  )

ggsave("AUC_forestB.tiff", plot = p_auc, width = 6.5, height = 3, dpi = 600)
######################################################PCA
library(ggplot2)
library(dplyr)

# 读数据
df <- read.csv("LC84AAb90-clusterB.csv")

# 设置分组
group <- factor(df[[4]], levels = c(0, 1), labels = c("NLC", "LC"))

# 选择表达矩阵并处理NA
selected_markers <- c(
  "CCL27","IL17A"
)
expr <- df[, selected_markers] 

#expr <- df[, 3:45]
expr[is.na(expr)] <- 0

# 标准化
expr_scaled <- scale(expr)

# PCA
pca_res <- prcomp(expr_scaled, center = TRUE, scale. = TRUE)
explained_var <- summary(pca_res)$importance[2, 1:2] * 100

# 构建数据框
pca_df <- data.frame(
  PC1 = pca_res$x[, 1],
  PC2 = pca_res$x[, 2],
  group = group
)

# 绘图：带填充的置信椭圆
pp <- ggplot(pca_df, aes(x = PC1, y = PC2, color = group, fill = group)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +  # 填充椭圆
  geom_point(size = 2, alpha = 1) +
  scale_color_manual(values = c("NLC" = "#fF4F39","LC" = "#2CA02C" )) + # 自定义点边框颜色
  scale_fill_manual(values = c("NLC" = "#fF4F39","LC" = "#2CA02C")) +   # 自定义椭圆填充颜色
  theme_bw(base_size = 14) +
  labs(
    title = NULL,
    x = paste0("PC1 (", round(explained_var[1], 1), "% variance)"),
    y = paste0("PC2 (", round(explained_var[2], 1), "% variance)"),
    color = "Group",
    fill = "Group"
  )
ggsave("pca_plotB.tiff", plot = pp, width = 4, height = 3, dpi = 600)
###################################################Figure3-6,FS1-2
#########################################################首先确定候选肽段及总体概览
########################################################
# 0️⃣ 读取数据
reads <- read.csv("reads.csv", header = TRUE, row.names = 1, check.names = FALSE)
mlxp <- read.csv("mlxp.csv", header = TRUE, row.names = 1, check.names = FALSE)

# 确保行列一致
stopifnot(all(rownames(reads) == rownames(mlxp)))
stopifnot(all(colnames(reads) == colnames(mlxp)))

# 定义阳性矩阵 (1=阳性, 0=阴性)
positive_matrix <- (reads > 100 & mlxp > 1.3) * 1
########################################################
# 1️⃣ 阳性率计算
library(reshape2)
library(dplyr)

positive_df <- as.data.frame(positive_matrix)
positive_df$Peptide <- rownames(positive_df)
positive_long <- melt(positive_df, id.vars = "Peptide",
                      variable.name = "Sample", value.name = "Positive")

# 样本分组（保留C/H/L/N）
sample_groups <- substr(colnames(positive_matrix), 1, 1)
names(sample_groups) <- colnames(positive_matrix)
positive_long$Group <- sample_groups[positive_long$Sample]

# 计算每个肽段阳性率
positive_rate <- positive_long %>%
  group_by(Peptide, Group) %>%
  summarise(
    PositiveRate = mean(Positive),
    PositiveN = sum(Positive),
    TotalN = n(),
    .groups = "drop"
  )

positive_rate_wide <- positive_rate %>%
  select(Peptide, Group, PositiveRate) %>%
  tidyr::pivot_wider(names_from = Group, values_from = PositiveRate, values_fill = 0)
########################################################总体与HC比较
# 2️⃣ 初筛：L > 0.1 或 N > 0.1或C>0.1
candidate_peptides <- positive_rate_wide %>%
  filter(L > 0.1 | N > 0.1| C > 0.1) %>%
  pull(Peptide)
######################################################
#####################################################zscore PCA 三个组健康对照 
library(ggplot2)
library(ggrepel)

# -------------------------------
# 1️⃣ PCA
# 三组：L/N/C
zscore <- read.csv("zscore.csv", row.names = 1, check.names = FALSE)
keep_samples <- substr(colnames(zscore),1,1) %in% c("L","N","C")
zscore_sub <- zscore[rownames(zscore) %in% candidate_peptides, keep_samples]

pca_res <- prcomp(t(zscore_sub), scale. = TRUE)

# 方差解释比例
explained_var <- summary(pca_res)$importance[2, 1:2] * 100

# -------------------------------
# 2️⃣ 构建 PCA 数据框
pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Group = substr(rownames(pca_res$x),1,1)
)
pca_df$Group[pca_df$Group=="L"] <- "LC"
pca_df$Group[pca_df$Group=="N"] <- "NLC"
pca_df$Group[pca_df$Group=="C"] <- "HC"

# -------------------------------
# 3️⃣ 绘图
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 1, linewidth = 1) +
  scale_color_manual(values = c("LC" = "#9f0117", "NLC" = "#377EB8", "HC" = "#fdae61")) +
  scale_fill_manual(values = c("LC" = "#9f0117", "NLC" = "#377EB8", "HC" = "#fdae61")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(
    title = "PCA ",
    x = paste0("PC1 (", round(explained_var[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_var[2], 1), "%)"),
    color = "Group",
    fill = "Group"
  )
#############################################阳性率t-SNE图
library(Rtsne)
library(vegan)
library(ggplot2)
library(dplyr)

# -----------------------------
# 1️⃣ 保留三组样本：L/N/C
keep_samples <- colnames(positive_matrix)[substr(colnames(positive_matrix),1,1) %in% c("L","N","C")]

# -----------------------------
# 2️⃣ 只保留 candidate_peptides
positive_sub <- positive_matrix[candidate_peptides, keep_samples]

# -----------------------------
# 3️⃣ 计算样本间 Jaccard 距离
dist_jaccard <- vegdist(t(positive_sub), method = "jaccard")

# -----------------------------
# 4️⃣ t-SNE 降维
set.seed(123)
tsne_res <- Rtsne(as.matrix(dist_jaccard), is_distance = TRUE, perplexity = 10)

# -----------------------------
# 5️⃣ 构建可视化数据框
tsne_df <- data.frame(
  Dim1 = tsne_res$Y[,1],
  Dim2 = tsne_res$Y[,2],
  Group = substr(colnames(positive_sub), 1, 1)
)
tsne_df$Group <- recode(tsne_df$Group,
                        L = "LC",
                        N = "NLC",
                        C = "HC")

# -----------------------------
# 6️⃣ 绘图（只画点，不要圈）
ggplot(tsne_df, aes(x=Dim1, y=Dim2, color=Group)) +
  geom_point(size=3, alpha=0.8) +
  scale_color_manual(values = c("LC"="#9f0117", "NLC"="#377EB8", "HC"="#fdae61")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(title="t-SNE (Jaccard distance, candidate peptides only)")
#######################################################差异蛋白热图
library(ComplexHeatmap)
library(circlize)

# -------------------------------
# 取出候选肽段对应数据
z_mat <- zscore[rownames(zscore) %in% candidate_peptides, ]
pos_mat <- positive_matrix[rownames(positive_matrix) %in% candidate_peptides, ]

# 根据列名获取样本组（首字母表示组）
sample_groups_orig <- substr(colnames(z_mat),1,1)

# 只保留 L/N/C 样本
keep_samples <- sample_groups_orig %in% c("L","N","C")
z_mat <- z_mat[, keep_samples]
pos_mat <- pos_mat[, keep_samples]
sample_groups_orig <- sample_groups_orig[keep_samples]

# 按组顺序排序：L/N/C
desired_order <- order(factor(sample_groups_orig, levels = c("L", "N", "C")))
z_mat <- z_mat[, desired_order]
pos_mat <- pos_mat[, desired_order]
sample_groups_orig <- sample_groups_orig[desired_order]

# 生成标签
sample_groups_label <- sample_groups_orig
sample_groups_label[sample_groups_label=="L"] <- "LC"
sample_groups_label[sample_groups_label=="N"] <- "NLC"
sample_groups_label[sample_groups_label=="C"] <- "HC"

# 设置颜色和注释
ha <- HeatmapAnnotation(
  Group = sample_groups_label,
  show_annotation_name = FALSE,
  col = list(Group = c(
    "LC" = "#4daf4a",
    "NLC" = "#984ea3",
    "HC" = "#fdae61"
  ))
)

# Z-score 热图
ht1 <- Heatmap(
  as.matrix(z_mat),
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#61AACF", "white", "#D6604D")),
  top_annotation = ha,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE
)

# Positive/Negative 热图
pos_mat_chr <- apply(pos_mat, 2, as.character)
ht2 <- Heatmap(
  pos_mat,#如果列想不聚类，改成pos_mat_chr
  name = "Positive",
  col = c("0"="#e2f6ff","1"="#57adde"),
  top_annotation = ha,
  cluster_columns = TRUE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = FALSE,
  heatmap_legend_param = list(at=c("0","1"), labels=c("Negative","Positive"))
)

# 拼接绘制
draw(ht1 + ht2, heatmap_legend_side="right", annotation_legend_side="right")
######################################################通过多轮筛选确定差异肽
########################################################
# 2️⃣ 初筛：L > 0.1 或 N > 0.1
candidate_peptides2 <- positive_rate_wide %>%
  filter(L > 0.1 | N > 0.1) %>%
  pull(Peptide)
################################################
########################################################
# 2.5️⃣ 二次筛选：L vs C 和 N vs C 的阳性率比较（Fisher 检验 + FDR）
compare_positive_FDR <- function(group1, group2, positive_long, candidate_peptides2) {
  
  df <- positive_long %>%
    filter(Peptide %in% candidate_peptides2 & Group %in% c(group1, group2)) %>%
    group_by(Peptide, Group) %>%
    summarise(
      PositiveN = sum(Positive),
      TotalN = n(),
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(
      names_from = Group,
      values_from = c(PositiveN, TotalN),
      names_sep = "_",
      values_fill = 0
    )
  
  results <- df %>%
    rowwise() %>%
    mutate(
      p_value = {
        pos1 <- get(paste0("PositiveN_", group1))
        pos2 <- get(paste0("PositiveN_", group2))
        neg1 <- get(paste0("TotalN_", group1)) - pos1
        neg2 <- get(paste0("TotalN_", group2)) - pos2
        
        # 判断表格是否有效（不能全0或全正）
        if ((pos1 + pos2) > 0 && (neg1 + neg2) > 0) {
          fisher.test(matrix(c(pos1, pos2, neg1, neg2), nrow = 2))$p.value
        } else {
          1
        }
      }
    ) %>%
    ungroup() %>%
    mutate(FDR = p.adjust(p_value, method = "fdr"),
           Comparison = paste0(group1, "_vs_", group2)) %>%
    select(Peptide, PositiveN_1 = paste0("PositiveN_", group1),
           TotalN_1 = paste0("TotalN_", group1),
           PositiveN_2 = paste0("PositiveN_", group2),
           TotalN_2 = paste0("TotalN_", group2),
           p_value, FDR, Comparison)
  
  return(results)
}

# 分别计算 L vs C 和 N vs C
res_LC <- compare_positive_FDR("L", "C", positive_long, candidate_peptides2)
res_NC <- compare_positive_FDR("N", "C", positive_long, candidate_peptides2)

########################################################
# 3️⃣ 输出 CSV
#write.csv(res_LC, "PositiveRate_L_vs_C.csv", row.names = FALSE)
#write.csv(res_NC, "PositiveRate_N_vs_C.csv", row.names = FALSE)

#cat("完成：输出 L_vs_C 和 N_vs_C 的阳性率比较表格。\n")
########################################################
# 4️⃣ 基于 FDR < 0.05 的二次筛选，得到新的候选肽段
candidate_peptides_LC <- res_LC$Peptide[res_LC$FDR < 0.05]
candidate_peptides_NC <- res_NC$Peptide[res_NC$FDR < 0.05]

# 合并两组候选肽段，去重
candidate_peptides<- unique(c(candidate_peptides_LC, candidate_peptides_NC))

cat("新的候选肽段数量:", length(candidate_peptides), "\n")

# 可选：保存新的候选肽段列表
#write.csv(data.frame(Peptide = candidate_peptides),
#          "Candidate_Peptides_LC_NC_FDR05.csv", row.names = FALSE)
################################################
# 3️⃣ 组间阳性率比较 (L vs N)
positive_rate_LN <- positive_long %>%
  filter(Peptide %in% candidate_peptides & Group %in% c("L","N")) %>%
  group_by(Peptide, Group) %>%
  summarise(
    PositiveRate = mean(Positive),
    PositiveN = sum(Positive),
    TotalN = n(),
    .groups = "drop"
  ) %>%
  tidyr::pivot_wider(names_from = Group, values_from = c(PositiveRate, PositiveN, TotalN),
                     names_sep = "_")

# Fisher检验
test_results <- positive_rate_LN %>%
  rowwise() %>%
  do({
    pos_L <- .$PositiveN_L; pos_N <- .$PositiveN_N
    neg_L <- .$TotalN_L - pos_L; neg_N <- .$TotalN_N - pos_N
    tab <- rbind(c(pos_L, pos_N), c(neg_L, neg_N))
    test <- fisher.test(tab)
    data.frame(Peptide = .$Peptide, p_value = test$p.value)
  }) %>%
  ungroup()

# 多重校正
test_results$FDR_rate <- p.adjust(test_results$p_value, method = "fdr")

########################################################
# 4️⃣ Z-score 差异 (L vs N)
zscore <- read.csv("zscore.csv", row.names = 1, check.names = FALSE)
zscore_positive <- zscore[rownames(zscore) %in% candidate_peptides, ]

# 样本分组（四组）
sample_groups_z <- substr(colnames(zscore_positive), 1, 1)
names(sample_groups_z) <- colnames(zscore_positive)

library(tibble)
zscore_long <- zscore_positive %>%
  as.data.frame() %>%
  rownames_to_column("Peptide") %>%
  tidyr::pivot_longer(-Peptide, names_to = "Sample", values_to = "Zscore") %>%
  mutate(Group = sample_groups_z[Sample]) %>%
  filter(Group %in% c("L","N"))

# t-test (Welch)
zscore_test <- zscore_long %>%
  group_by(Peptide) %>%
  summarise(
    p_value_zscore = t.test(Zscore ~ Group, var.equal = FALSE)$p.value,
    .groups = "drop"
  )

zscore_test$FDR_zscore <- p.adjust(zscore_test$p_value_zscore, method = "fdr")

########################################################
# 5️⃣ 取交集最终差异肽
sig_peptides_rate <- test_results$Peptide[test_results$p_value < 0.05]
sig_peptides_zscore <- zscore_test$Peptide[zscore_test$p_value_zscore < 0.05]
final_diff_peptides <- intersect(sig_peptides_rate, sig_peptides_zscore)
#final_diff_peptides <- sig_peptides_rate
########################################################
# 5️⃣ 计算Zscore平均值并得到FC（取绝对值）
zscore_summary <- data.frame(
  Peptide = rownames(zscore_positive),
  Z_mean_L = apply(zscore_positive[, sample_groups_z=="L"], 1, mean, na.rm=TRUE),
  Z_mean_N = apply(zscore_positive[, sample_groups_z=="N"], 1, mean, na.rm=TRUE)
) %>%
  mutate(
    # 计算FC（取绝对值，防止负号或除0）
    FC = ifelse(is.na(Z_mean_N) | Z_mean_N == 0, NA, abs(Z_mean_L / Z_mean_N)),
    log2FC = log2(FC),
    abs_log2FC = abs(log2FC)
  )

########################################################
# 6️⃣ 综合筛选（保持与手动标准一致）
# 条件：
#   (1) Fisher检验 p<0.05
#   (2) t检验 p<0.05
#   (3) FC > 1.5 或 FC < 0.667（即 |log2FC| > 0.585）
final_diff_peptides <- test_results %>%
  left_join(zscore_test, by = "Peptide") %>%
  left_join(zscore_summary, by = "Peptide") %>%
  filter(
    !is.na(FC),
    p_value < 0.05,
    p_value_zscore < 0.05,
    (FC > 1.5 | FC < 0.667)
  ) %>%
  pull(Peptide)
#write.csv(final_diff_peptides, "diff_peptide.csv", row.names = FALSE)

########################################################
# 输出完整结果表
final_table <- test_results %>%
  left_join(zscore_test, by = "Peptide") %>%
  left_join(positive_rate_wide, by = "Peptide") %>%
  left_join(zscore_summary, by = "Peptide") %>%
  filter(Peptide %in% final_diff_peptides & !is.na(FC)) %>%
  mutate(Direction = ifelse(Z_mean_L > Z_mean_N, "Up_in_L", "Down_in_L"))

#write.csv(final_table, "final_diff_peptides_full_info_LN_FC_abs.csv", row.names = FALSE)
#######################################################
#####################################################差异蛋白PCA 三个组健康对照 
library(ggplot2)
library(ggrepel)

# -------------------------------
# 1️⃣ PCA
# 三组：L/N/C
zscore <- read.csv("zscore.csv", row.names = 1, check.names = FALSE)
keep_samples <- substr(colnames(zscore),1,1) %in% c("L","N","C")
zscore_sub <- zscore[rownames(zscore) %in% final_diff_peptides, keep_samples]

pca_res <- prcomp(t(zscore_sub), scale. = TRUE)

# 方差解释比例
explained_var <- summary(pca_res)$importance[2, 1:2] * 100

# -------------------------------
# 2️⃣ 构建 PCA 数据框
pca_df <- data.frame(
  Sample = rownames(pca_res$x),
  PC1 = pca_res$x[,1],
  PC2 = pca_res$x[,2],
  Group = substr(rownames(pca_res$x),1,1)
)
pca_df$Group[pca_df$Group=="L"] <- "LC"
pca_df$Group[pca_df$Group=="N"] <- "NLC"
pca_df$Group[pca_df$Group=="C"] <- "HC"

# -------------------------------
# 3️⃣ 绘图
ggplot(pca_df, aes(x = PC1, y = PC2, color = Group, fill = Group)) +
  stat_ellipse(geom = "polygon", alpha = 0.2, level = 0.95) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, linetype = 1, size = 1) +
  scale_color_manual(values = c("LC" = "#9f0117", "NLC" = "#377EB8", "HC" = "#fdae61")) +
  scale_fill_manual(values = c("LC" = "#9f0117", "NLC" = "#377EB8", "HC" = "#fdae61")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(
    title = "PCA of final differential peptides",
    x = paste0("PC1 (", round(explained_var[1], 1), "%)"),
    y = paste0("PC2 (", round(explained_var[2], 1), "%)"),
    color = "Group",
    fill = "Group"
  )
####################################################tsne 阳性率
library(Rtsne)
library(vegan)
library(ggplot2)
library(dplyr)

# -----------------------------
# 0️⃣ 保留三组样本：L/N/C
keep_samples <- colnames(positive_matrix)[substr(colnames(positive_matrix), 1, 1) %in% c("L","N","C")]
positive_sub <- positive_matrix[final_diff_peptides, keep_samples, drop=FALSE]

# -----------------------------
# 1️⃣ 去掉全阴性肽（行）和全阴性样本（列）
positive_sub <- positive_sub[rowSums(positive_sub) > 0, ]
positive_sub <- positive_sub[, colSums(positive_sub) > 0]

# -----------------------------
# 2️⃣ 计算样本间 Jaccard 距离
dist_jaccard <- vegdist(t(positive_sub), method = "jaccard")

# -----------------------------
# 3️⃣ t-SNE 降维
set.seed(666)#1111
tsne_res <- Rtsne(as.matrix(dist_jaccard), is_distance = TRUE, perplexity = 10)

# -----------------------------
# 4️⃣ 构建可视化数据框
tsne_df <- data.frame(
  Dim1 = tsne_res$Y[,1],
  Dim2 = tsne_res$Y[,2],
  Sample = colnames(positive_sub)
)
tsne_df$Group <- substr(tsne_df$Sample, 1, 1)
tsne_df$Group <- recode(tsne_df$Group,
                        L = "LC",
                        N = "NLC",
                        C = "HC")

# -----------------------------
# 5️⃣ 绘图（只画点，不要圈）
ggplot(tsne_df, aes(x = Dim1, y = Dim2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  scale_color_manual(values = c("LC"="#9f0117", "NLC"="#377EB8", "HC"="#fdae61")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(title = "t-SNE (Jaccard distance)")
#带圈
ggplot(tsne_df, aes(x = Dim1, y = Dim2, color = Group)) +
  geom_point(size = 3, alpha = 0.8) +
  stat_ellipse(level = 0.95, size = 1, linetype = 2, alpha = 0.6) +
  scale_color_manual(values = c("LC"="#9f0117", "NLC"="#377EB8", "HC"="#fdae61")) +
  theme_bw(base_size = 14) +
  theme(panel.grid = element_blank()) +
  labs(title = "t-SNE (Jaccard distance)")
#######################################################差异蛋白热图
library(ComplexHeatmap)
library(circlize)

# -------------------------------
# 取出候选肽段对应数据
z_mat <- zscore[rownames(zscore) %in% final_diff_peptides, ]
pos_mat <- positive_matrix[rownames(positive_matrix) %in% final_diff_peptides, ]

# 根据列名获取样本组（首字母表示组）
sample_groups_orig <- substr(colnames(z_mat),1,1)

# 只保留 L/N/C 样本
keep_samples <- sample_groups_orig %in% c("L","N","C")
z_mat <- z_mat[, keep_samples]
pos_mat <- pos_mat[, keep_samples]
sample_groups_orig <- sample_groups_orig[keep_samples]

# 按组顺序排序：L/N/C
desired_order <- order(factor(sample_groups_orig, levels = c("L", "N", "C")))
z_mat <- z_mat[, desired_order]
pos_mat <- pos_mat[, desired_order]
sample_groups_orig <- sample_groups_orig[desired_order]

# 生成标签
sample_groups_label <- sample_groups_orig
sample_groups_label[sample_groups_label=="L"] <- "LC"
sample_groups_label[sample_groups_label=="N"] <- "NLC"
sample_groups_label[sample_groups_label=="C"] <- "HC"

# 设置颜色和注释
ha <- HeatmapAnnotation(
  Group = sample_groups_label,
  show_annotation_name = FALSE,
  col = list(Group = c(
    "LC" = "#4daf4a",
    "NLC" = "#984ea3",
    "HC" = "#fdae61"
  ))
)

# Z-score 热图
ht1 <- Heatmap(
  as.matrix(z_mat),
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#61AACF", "white", "#D6604D")),
  top_annotation = ha,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE
)

# Positive/Negative 热图
pos_mat_chr <- apply(pos_mat, 2, as.character)
ht2 <- Heatmap(
  pos_mat_chr,
  name = "Positive",
  col = c("0"="#e2f6ff","1"="#57adde"),
  top_annotation = ha,
  cluster_columns = FALSE,
  cluster_rows = TRUE,
  show_column_names = FALSE,
  show_row_names = TRUE,
  heatmap_legend_param = list(at=c("0","1"), labels=c("Negative","Positive"))
)

# 拼接绘制
draw(ht1 + ht2, heatmap_legend_side="right", annotation_legend_side="right")
########################################################LASSO和RFE
########################################################
library(mixOmics)
library(glmnet)
library(caret)
library(pROC)
library(dplyr)
library(ggplot2)

# ======================================================
# 1️⃣ 数据准备 & PLS-DA VIP筛选
# ======================================================
peptides2 <- final_diff_peptides
zscore_sub <- zscore[rownames(zscore) %in% peptides2, ]
zscore_sub_t <- t(zscore_sub)

# 分组标签：L=1, N=0
labels <- ifelse(grepl("^L", rownames(zscore_sub_t)), 1,
                 ifelse(grepl("^N", rownames(zscore_sub_t)), 0, NA))
valid_idx <- !is.na(labels)
zscore_sub_t <- zscore_sub_t[valid_idx, , drop=FALSE]
labels <- labels[valid_idx]
table(labels)

# PLS-DA
X <- zscore_sub_t
Y <- factor(labels, levels=c(0,1))
plsda_res <- mixOmics::plsda(X, Y, ncomp = 2)
# VIP筛选
vip_scores <- vip(plsda_res)
vip_selected <- rownames(vip_scores)[vip_scores[,1] > 1]
cat("✅ VIP>1 特征数:", length(vip_selected), "\n")

zscore_vip <- zscore_sub_t[, vip_selected, drop=FALSE]

# 准备数据
vip_df <- data.frame(
  Feature = rownames(vip_scores),
  VIP = vip_scores[,1]
)

# 排序从大到小
vip_df <- vip_df[order(vip_df$VIP, decreasing=FALSE), ]  # 注意这里 decreasing=FALSE
vip_df$Feature <- factor(vip_df$Feature, levels = vip_df$Feature)

# 绘图
ggplot(vip_df, aes(x=Feature, y=VIP)) +
  geom_col(fill="#9fc7b4") +
  geom_hline(yintercept=1, color="red", linetype="dashed", size=1) +
  coord_flip() +
  labs(x="Feature", y="VIP Score", title="PLS-DA VIP Scores") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size=9,colour = "black",face = "bold"),
    plot.title = element_text(hjust=0.5)
  )
# ======================================================
# 2️⃣ LASSO Logistic 回归 + RFE
# ======================================================
features_current <- colnames(zscore_vip)
rfe_results <- data.frame(
  n_features = integer(),
  auc_mean = numeric(),
  auc_lower = numeric(),
  auc_upper = numeric()
)
rfe_features_all <- data.frame()  # 用于保存所有轮次特征及系数

step_large <- 3
step_small <- 1
n_repeat_rfe <- 20
nfold_rfe <- 5

while(length(features_current) > 0){
  auc_values <- c()
  
  # 交叉验证计算平均AUC
  for(seed in 1:n_repeat_rfe){
    set.seed(seed)
    folds <- createFolds(labels, k = nfold_rfe)
    
    for(fold_idx in 1:nfold_rfe){
      test_idx <- folds[[fold_idx]]
      train_idx <- setdiff(1:length(labels), test_idx)
      if(length(unique(labels[test_idx])) < 2) next
      
      x_train <- as.matrix(zscore_vip[train_idx, features_current, drop=FALSE])
      y_train <- labels[train_idx]
      x_test <- as.matrix(zscore_vip[test_idx, features_current, drop=FALSE])
      y_test <- labels[test_idx]
      
      cvfit <- cv.glmnet(x_train, y_train, family="binomial", alpha=1, nfolds=5)
      pred_prob <- predict(cvfit, x_test, s="lambda.min", type="response")
      auc_values <- c(auc_values, auc(roc(y_test, as.numeric(pred_prob))))
    }
  }
  
  auc_mean <- mean(auc_values)
  auc_ci <- quantile(auc_values, probs = c(0.025, 0.975))
  nf <- length(features_current)
  
  cat("✅", nf, "features:", sprintf("%.3f [%.3f - %.3f]\n", auc_mean, auc_ci[1], auc_ci[2]))
  
  rfe_results <- rbind(rfe_results, data.frame(
    n_features = nf,
    auc_mean = auc_mean,
    auc_lower = auc_ci[1],
    auc_upper = auc_ci[2]
  ))
  
  # 用全数据计算当前特征的LASSO系数
  fit_full <- cv.glmnet(as.matrix(zscore_vip[, features_current, drop=FALSE]),
                        labels, family="binomial", alpha=1, nfolds=5)
  coefs <- coef(fit_full, s="lambda.min")[-1,]  # 去掉截距
  coefs_abs <- abs(coefs)
  
  # 保存当前轮特征和系数到总表
  rfe_features_all <- rbind(rfe_features_all,
                            data.frame(
                              Round = paste0("n=", nf),
                              Feature = names(coefs),
                              Coefficient = as.numeric(coefs),
                              AbsCoefficient = coefs_abs
                            ))
  
  # 每轮删除特征数
  remove_n <- ifelse(nf > 15, step_large, step_small)
  
  # 按系数绝对值从小到大排序，删除贡献最小的特征
  to_remove <- names(sort(coefs_abs, decreasing = FALSE))[1:remove_n]
  
  features_current <- setdiff(features_current, to_remove)
  
  if(length(features_current) <= 1) break
}

# ============================================
# 导出所有轮次特征及系数到一个 CSV
write.csv(rfe_features_all, "RFE_all_rounds_features.csv", row.names = FALSE)
# ======================================================
# 3️⃣ 可视化 RFE 结果
# ======================================================
library(ggplot2)

# 自定义横坐标刻度
breaks_large <- seq(max(rfe_results$n_features), 16, by = -5)  # 大于15，每5个
breaks_small <- 15:1                                       # 小于等于15，每1个
x_breaks <- c(breaks_large, breaks_small)

ggplot(rfe_results, aes(x=n_features, y=auc_mean)) +
  geom_line(color="#1f77b4", size=1) +
  geom_point(color="#1f77b4", size=2) +
  geom_ribbon(aes(ymin=auc_lower, ymax=auc_upper), alpha=0.2, fill="#1f77b4") +
  geom_vline(xintercept = 6, color="red", linetype="dashed", size=1) +  # 红色竖线
  scale_x_reverse(breaks = x_breaks) +
  labs(x="Number of Features", y="Mean CV AUC",
       title="LASSO-RFE Performance (VIP>1 → RFE)") +
  theme_bw()
# ============================================
# 导出最终前9个最重要特征名称
# 找到剩余特征数为9的轮次
# 去掉 "n="，转换为数字
rfe_features_all$RoundNum <- as.numeric(sub("n=", "", rfe_features_all$Round))

# 假设你要导出 n=11 的特征
round_6 <- rfe_features_all[rfe_features_all$RoundNum == 6, ]

top6_features <- round_6[, c("Feature", "Coefficient", "AbsCoefficient")]
write.csv(top6_features, "RFE_top6_features.csv", row.names = FALSE)
df_top6 <- read.csv("Top6.csv", stringsAsFactors = FALSE)
#可视化
ggplot(df_top6, aes(x = AbsCoefficient, y = reorder(Feature, AbsCoefficient))) +
  geom_col(fill = "#ff7f0e") +
  labs(x = "Absolute LASSO Coefficient", y = "Feature",
       title = "Top 6 Features (LASSO)") +
  theme_bw()
#####################################ROC
library(glmnet)
library(pROC)
library(ggplot2)
library(dplyr)
library(caret)

# =============================
# 1️⃣ 数据准备
# =============================
top6_features <- read.csv("RFE_top6_features.csv", stringsAsFactors = FALSE)
features <- top6_features$Feature
x <- as.matrix(zscore_vip[, features])
y <- labels

# =============================
# 2️⃣ 重复交叉验证预测
# =============================
nfold <- 5
n_repeat <- 50
set.seed(123)

pred_all_matrix <- matrix(NA, nrow=length(y), ncol=n_repeat)

for(seed in 1:n_repeat){
  set.seed(seed)
  folds <- createFolds(y, k = nfold)
  pred_all <- rep(NA, length(y))
  
  for(f in 1:nfold){
    test_idx <- folds[[f]]
    train_idx <- setdiff(1:length(y), test_idx)
    
    x_train <- x[train_idx, ]
    y_train <- y[train_idx]
    x_test <- x[test_idx, ]
    
    # LASSO逻辑回归
    fit <- cv.glmnet(x_train, y_train, family="binomial", alpha=1, nfolds=5)
    pred_prob <- predict(fit, x_test, s="lambda.min", type="response")
    
    pred_all[test_idx] <- pred_prob
  }
  pred_all_matrix[, seed] <- pred_all
}

# 平均预测概率
pred_mean <- rowMeans(pred_all_matrix)

# =============================
# 3️⃣ ROC + 置信区间计算
# =============================
roc_obj <- roc(y, pred_mean)
auc_val <- auc(roc_obj)

# 每个FPR的TPR置信区间（兼容各种 pROC 版本）
ci_obj <- ci.se(roc_obj, specificities = seq(0, 1, l = 50), boot.n = 1000, conf.level = 0.95)

# 将 ci.se 返回值转换为数据框
ci_df <- as.data.frame(ci_obj)
colnames(ci_df) <- c("sens_low", "sens_mid", "sens_high")
ci_df$specificity <- seq(0, 1, l = nrow(ci_df))

# ROC曲线数据
roc_df <- data.frame(
  specificity = rev(roc_obj$specificities),
  sensitivity = rev(roc_obj$sensitivities)
)

# AUC文本
auc_val <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj, conf.level = 0.95)

auc_text <- paste0("AUC = ", round(auc_val, 3),
                   " [", round(auc_ci[1], 3), "-", round(auc_ci[3], 3), "]")

# =============================
# 4️⃣ 绘图
# =============================
ggplot() +
  # CI 阴影
  geom_ribbon(data = ci_df,
              aes(x = 1 - specificity, ymin = sens_low, ymax = sens_high),
              fill = "#FF1493", alpha = 0.2) +
  # ROC 曲线
  geom_line(data = roc_df,
            aes(x = 1 - specificity, y = sensitivity),
            color = "#FF1493", size = 1.2) +
  # 对角线
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "grey50") +
  labs(
    x = "False Positive Rate",
    y = "True Positive Rate",
    title = "Cross-validated ROC"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line    = element_line(colour = "black")
  ) +
  # 添加 AUC + CI 到右下角
  annotate("text", 
           x = 0.95, y = 0.05,
           label = auc_text, 
           hjust = 1, vjust = 0,
           size = 5, fontface = "bold", color = "black")
##############################################################3
############加载总体差异肽.RData
############################################按照阳性率+zscore分类一致性聚类
library(ConsensusClusterPlus)
library(cluster)  # daisy() 计算 Gower 距
# -------------------------------
# 1️⃣ 取 L/N 样本
LN_samples <- colnames(positive_matrix)[substr(colnames(positive_matrix),1,1) %in% c("L","N")]

# -------------------------------
# 2️⃣ 提取差异肽矩阵
LN_pos <- positive_matrix[rownames(positive_matrix) %in% candidate_peptides2, LN_samples]
LN_z   <- zscore[rownames(zscore) %in% candidate_peptides2, LN_samples]
# -------------------------------
# 3️⃣ 转置为样本 × 特征
pos_t <- t(LN_pos)
z_t   <- t(LN_z)
# -------------------------------
# 4️⃣ 合并为一个矩阵（连续 + 二分类）
# 注意：二分类变量需要转成 factor
pos_t_factor <- as.data.frame(pos_t)
pos_t_factor[] <- lapply(pos_t_factor, factor)
z_t_df <- as.data.frame(z_t)
dataMatrix <- cbind(z_t_df, pos_t_factor)  # 样本 × 特征
# -------------------------------
# 5️⃣ 计算 Gower 距离
gower_dist <- daisy(dataMatrix, metric = "gower")
# -------------------------------
# 6️⃣ 一致性聚类
outputDir <- "Consensus_final_diff_Gower"
if(!dir.exists(outputDir)) dir.create(outputDir)
maxK <- 6
reps <- 100
results <- ConsensusClusterPlus(
  d = as.matrix(gower_dist),
  maxK = maxK,
  reps = reps,
  pItem = 0.8,
  pFeature = 1,
  clusterAlg = "pam",
  distance = "euclidean",
  seed = 123,
  title = outputDir,
  plot = "png"
)
# 7️⃣ 计算 ICL
icl <- calcICL(results, title = outputDir, plot = "png")
# 8️⃣ 选择聚类数，比如 K=2
clusterNum <- 2
cluster <- results[[clusterNum]][["consensusClass"]]
# -------------------------------
# 9️⃣ 整理结果
cluster_df <- data.frame(Sample = names(cluster), Cluster = cluster)
letters_vec <- c("A","B","C","D","E","F")
uniqClu <- sort(unique(cluster_df$Cluster))
cluster_df$Cluster <- letters_vec[match(cluster_df$Cluster, uniqClu)]
# -------------------------------
# 10️⃣ 保存结果
write.csv(cluster_df, file = file.path(outputDir, "LN_candidate_cluster_assignments.csv"), row.names = FALSE)
cat("✅ 基于阳性率+Jaccard 距离的 L/N 一致性聚类完成，结果保存在:", outputDir, "\n")
#################################################################热图
library(ComplexHeatmap)
library(circlize)
# -------------------------------
# 1️⃣ 取 L/N 样本
LN_samples <- colnames(positive_matrix)[substr(colnames(positive_matrix),1,1) %in% c("L","N")]
LN_matrix <- positive_matrix[rownames(positive_matrix) %in% candidate_peptides2, LN_samples]
# -------------------------------
# 2️⃣ 去掉全零或全一的肽（行）
non_constant_rows <- apply(LN_matrix, 1, function(x) length(unique(x)) > 1)
LN_matrix_filtered <- LN_matrix[non_constant_rows, ]
# -------------------------------
# 3️⃣ 确保数值型
LN_matrix_filtered <- apply(LN_matrix_filtered, 2, as.numeric)
# -------------------------------
# 4️⃣ 样本分组信息
sample_clusters <- cluster_df$Cluster[match(colnames(LN_matrix_filtered), cluster_df$Sample)]
sample_groups <- ifelse(substr(colnames(LN_matrix_filtered),1,1)=="L", "LC", "NLC")
# -------------------------------
# 5️⃣ 按 Cluster(A/B) → Group(LC/NLC) 排序
order_idx <- order(sample_clusters, sample_groups)
LN_matrix_filtered <- LN_matrix_filtered[, order_idx]
sample_clusters <- sample_clusters[order_idx]
sample_groups <- sample_groups[order_idx]
# -------------------------------
# 6️⃣ 列注释
ha <- HeatmapAnnotation(
  Cluster = sample_clusters,
  Group = sample_groups,
  col = list(
    Cluster = c("A"="#FF8C42","B"="#377EB8"),
    Group = c("LC"="#4daf4a","NLC"="#984ea3")
  ),
  simple_anno_size = unit(2.5, "mm"),
  gap = unit(0.4, "mm"),
  show_annotation_name = TRUE
)
# -------------------------------
# 7️⃣ 绘制热图
ht <- Heatmap(
  LN_matrix_filtered,
  name = "Positive",
  col = c("0"="#e2f6ff", "1"="#57adde"),
  top_annotation = ha,
  cluster_rows = TRUE,       # 行聚类
  cluster_columns = FALSE,   # 列不聚类，按上面顺序
  show_row_names = FALSE,    # 不显示肽名称
  show_column_names = FALSE,
  heatmap_legend_param=list(at=c("0","1"), labels=c("Negative","Positive"))
)
draw(ht, heatmap_legend_side="right", annotation_legend_side="right")
######################################################
####################################################富集分析
#################################富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 2. 准备你的差异基因列表（Symbol）
gene_df <- read.csv("gene1.csv", header = FALSE, stringsAsFactors = FALSE)
genes <- na.omit(gene_df[[1]])
gene_list <- genes

# 3. ID 转换：将基因 Symbol 转换为 Entrez Gene ID
# org.Hs.eg.db 是人类基因注释数据库
# 从 'SYMBOL' 列映射到 'ENTREZID' 列
entrez_ids <- bitr(geneID = gene_list,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# 4. 进行 GO 和 KEGG 富集分析
# GO - 生物过程 (Biological Process, BP)
ego_bp <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)

# GO - 细胞组分 (Cellular Component, CC)
ego_cc <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05
)
# GO - 分子功能 (Molecular Function, MF)
ego_mf <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)
# KEGG 富集分析
options(timeout = 300)
ekegg <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                    organism     = 'hsa', # 人类 KEGG 库
                    pvalueCutoff = 0.05)
# 5. 导出所有 GO 和 KEGG 富集结果为 CSV 文件
# 导出 GO-BP 结果
if (!is.null(ego_bp)) {
  ego_bp_df <- as.data.frame(ego_bp@result)
  write.csv(ego_bp_df, "GO_BP_enrichment_results.csv", row.names = FALSE)
  print("GO生物过程富集结果已导出至 'GO_BP_enrichment_results.csv'")
} else {
  print("没有显著的GO生物过程富集通路，未导出文件。")
}
# 导出 GO-CC 结果
if (!is.null(ego_cc)) {
  ego_cc_df <- as.data.frame(ego_cc@result)
  write.csv(ego_cc_df, "GO_CC_enrichment_results.csv", row.names = FALSE)
  print("GO细胞组分富集结果已导出至 'GO_CC_enrichment_results.csv'")
} else {
  print("没有显著的GO细胞组分富集通路，未导出文件。")
}

# 导出 GO-MF 结果
if (!is.null(ego_mf)) {
  ego_mf_df <- as.data.frame(ego_mf@result)
  write.csv(ego_mf_df, "GO_MF_enrichment_results.csv", row.names = FALSE)
  print("GO分子功能富集结果已导出至 'GO_MF_enrichment_results.csv'")
} else {
  print("没有显著的GO分子功能富集通路，未导出文件。")
}

# 导出 KEGG 结果
if (!is.null(ekegg) && nrow(as.data.frame(ekegg@result)) > 0) {
  ekegg_df <- as.data.frame(ekegg@result)
  write.csv(ekegg_df, "KEGG_enrichment_results.csv", row.names = FALSE)
  print("KEGG富集结果已导出至 'KEGG_enrichment_results.csv'")
} else {
  print("KEGG富集结果为空，未导出文件。")
}
#############################################################KEGGda大分类
## 绘制 KEGG 功能分类的流星图
library(ggplot2)
library(ggforce)
library(readr)
library(dplyr)
# 读取 KEGG 数据
df <- read_csv("kegg.csv")
colors <- c(
  "Metabolism" = "#1F78B4",
  "Genetic Information Processing" = "#33A02C",
  "Environmental Information Processing" = "#E31A1C",
  "Cellular Processes" = "#FF7F00",
  "Organismal Systems" = "#6A3D9A",
  "Human Diseases" = "#B15928"
)
# 排序（按显著性）
df <- df %>%
  arrange(pvalue) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# 绘制流星图
p <- ggplot(df) +
  geom_link(aes(x = 0, y = Description,
                xend = -log10(pvalue), yend = Description,
                alpha = after_stat(index),
                color = Category,
                linewidth = after_stat(index)),
            n = 500, show.legend = c(color = TRUE, alpha = FALSE, linewidth = FALSE)) 
p1 <- p +
  geom_point(aes(x = -log10(pvalue), y = Description, fill = Category),
             color = "black",        # 黑边
             size = 6, shape = 21) + # 流星头
  geom_text(aes(x = -log10(pvalue), y = Description, label = Count),
            size = 3) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    strip.text = element_text(face = "bold.italic"),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    axis.text.x = element_text(face = "bold", size = 11, color = "black"),
    axis.title = element_text(size = 12),
    legend.key.size = unit(1, "lines"),
    legend.text = element_text(size = 10, face = "bold"),
    legend.title = element_text(size = 10, face = "bold"),
    legend.position = c(0.81, 0.7)  # legend 位置调整，靠近图，(x, y) 坐标，0-1
  ) +
  xlab("-Log10 Pvalue") + 
  ylab("") +
  # 图例方块和颜色一致
  scale_fill_manual(values = colors, 
                    guide = guide_legend(override.aes = list(shape = 22, size = 5))) +
  scale_color_manual(values = colors)  # 尾巴颜色和点填充颜色一致
# 保存结果
ggsave(filename = "KEGG_cometplot.tiff",
       width = 12, height = 7, plot = p1)
########################################################纵坐标黑色GO
library(ggplot2)
library(dplyr)
library(forcats)
# 读取数据
go_df <- read.csv("GO.csv", stringsAsFactors = FALSE)
# 按显著性排序 GO 条目
go_df <- go_df %>%
  group_by(Category) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  ungroup()
# 将 GO 条目按照 -log10(pvalue) 排序
go_df <- go_df %>%
  mutate(Description = factor(Description, levels = unique(Description)))
# 定义颜色
colors <- c(
  "Biological Process" = "#1f77b4",
  "Molecular Function" = "#ff7f0e",
  "Cellular Component" = "#2ca02c"
)
# 绘图
p <- ggplot(go_df, aes(x = -log10(pvalue), y = Description)) +
  geom_point(aes(size = Count, fill = Category), shape = 21, color = "black") +
  scale_fill_manual(values = colors) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw() +   # 四边框
  theme(
    axis.text.y = element_text(face = "bold", size = 13, color = "black"),  # 改成黑色
    axis.text.x = element_text(face = "bold", size = 13, color = "black"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold")
  ) +
  ylab("") +
  xlab("-Log10 P-value") +
  guides(
    fill = guide_legend(title = "Category", override.aes = list(size = 5))
  )
ggsave("GO_triplet_bubbleplot_colored_y1.tiff", plot = p, width = 11, height = 8)
#################################################自身抗体两簇间细胞因子比较
#####################################################
##########################自身抗体两个簇之间比较
library(gtsummary)
library(dplyr)
data <- read.csv("LC84AAb90-cluster.csv")
data <- data[3:68]
tbl_summary(data) #缺点是默认连续型变量不符合正态分布，统一采用非参数检验，分类变量统一采用皮尔逊卡方检验；优点是考虑缺失值并单独列出来
tbl_summary(data, by = autoantibody) # 根据LCC分层
res1 <- tbl_summary(data, by = autoantibody) %>% add_p()   #添加p值
#write.csv(res1, "autoantibody_AB.csv", row.names = FALSE)
####################################################自身抗体A簇
##########################自身抗体A簇两个组之间比较
library(gtsummary)
library(dplyr)
data <- read.csv("LC84AAb90-clusterA.csv")
data <- data[3:68]
tbl_summary(data) #缺点是默认连续型变量不符合正态分布，统一采用非参数检验，分类变量统一采用皮尔逊卡方检验；优点是考虑缺失值并单独列出来
tbl_summary(data, by = LC1) # 根据LCC分层
res1 <- tbl_summary(
  data,
  by = LC1,
  type = list(CSF3 ~ "continuous")
) %>%
  add_p()  #添加p值
#write.csv(res1, "autoantibodyA_LN.csv", row.names = FALSE)
########################################################3火山图
# 加载必要包
library(dplyr)
library(ggplot2)
library(ggrepel)

# 读取数据
df <- read.csv("LC84AAb90-clusterA.csv")
group <- df[[4]]  # 第二列为分组 LC（0/1）
expr <- df[, 5:47]  # 第3至第49列为细胞因子

# 差异分析 + log2FC 计算
res <- lapply(colnames(expr), function(gene) {
  group0 <- expr[group == 0, gene]
  group1 <- expr[group == 1, gene]
  test <- wilcox.test(group1, group0)
  log2fc <- log2(median(group1, na.rm = TRUE) / median(group0, na.rm = TRUE))
  data.frame(
    marker = gene,
    log2FC = log2fc,
    pvalue = test$p.value
  )
}) %>% bind_rows()

# 根据 p 值和 log2FC 标记显著性
res$threshold <- with(res, ifelse(pvalue < 0.05 & abs(log2FC) > 0.15,
                                  ifelse(log2FC > 0.15, "Up", "Down"), "NS"))

# 火山图绘制（基于 p 值）
ggplot(res, aes(x = log2FC, y = -log10(pvalue), color = threshold, label = marker)) +
  geom_point(size = 3, alpha = 0.8) +
  geom_text_repel(data = subset(res, threshold != "NS"), size = 3.5, max.overlaps = 15) +
  scale_color_manual(values = c("Up" = "#D73027", "Down" = "#4575B4", "NS" = "grey60")) +
  geom_vline(xintercept =0, linetype = "dashed") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed") +
  xlim(-2, 2) +  # 限制横坐标范围
  theme_bw(base_size = 14) +
  labs(title = NULL,
       x = "log2(Fold Change)",
       y = "-log10(p-value)",
       color = "")
#######################################################调整变量的单因素分析
library(logistf)
library(dplyr)
library(tibble)
library(ggplot2)

data <- read.csv("LC84AAb90-clusterA.csv")
data <- data[4:68]

data$LC1 <- as.factor(data$LC1)
data$treattype <- as.factor(data$treattype)

bad_vars <- names(data)[sapply(data, function(x) n_distinct(x) < 2)]
data <- data %>% select(-all_of(bad_vars))

vars <- setdiff(names(data), c("LC1", "treattype"))

# 初始化结果表
res <- tibble(
  Variable = character(),
  OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric(),
  Adjusted = logical()  # 是否调整treattype
)

# 1️⃣ 调整treattype
for (v in vars) {
  formula <- as.formula(paste("LC1 ~", v, "+ treattype"))
  
  fit <- tryCatch({
    logistf(formula, data = data)
  }, error = function(e) NULL)
  
  if(!is.null(fit)) {
    est <- fit$coefficients[v]
    ci <- confint(fit)[v, ]
    p <- fit$prob[v]
    
    res <- res %>% add_row(
      Variable = v,
      OR = exp(est),
      CI_lower = exp(ci[1]),
      CI_upper = exp(ci[2]),
      p_value = p,
      Adjusted = TRUE
    )
  }
}

# 2️⃣ 未调整treattype
for (v in vars) {
  formula <- as.formula(paste("LC1 ~", v))
  
  fit <- tryCatch({
    logistf(formula, data = data)
  }, error = function(e) NULL)
  
  if(!is.null(fit)) {
    est <- fit$coefficients[v]
    ci <- confint(fit)[v, ]
    p <- fit$prob[v]
    
    res <- res %>% add_row(
      Variable = v,
      OR = exp(est),
      CI_lower = exp(ci[1]),
      CI_upper = exp(ci[2]),
      p_value = p,
      Adjusted = FALSE
    )
  }
}

# 筛选显著因子（p<0.05）用于绘图
res_plot <- res %>%
  filter(!is.na(OR) & p_value < 0.05) %>%
  mutate(Variable = factor(Variable, levels = rev(unique(Variable))))

# 绘制森林图：调整treattype的红色，未调整的蓝色
ggplot(res_plot, aes(x = OR, y = Variable, color = Adjusted)) +
  geom_point(size = 3, position = position_dodge(width = 0.6)) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2,
                 position = position_dodge(width = 0.6)) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  scale_color_manual(values = c("TRUE" = "red", "FALSE" = "blue"),
                     labels = c("Adjusted,\np<0.05", "Unadjusted,\np<0.05")) +
  theme_bw() +
  labs(x = "Odds Ratio (log scale)", y = "", color = "") +
  theme(axis.text.y = element_text(size = 12,face = "bold",color = "black"),
        axis.text.x = element_text(size = 12,face = "bold",color = "black"),
        legend.key.height = unit(2.5, "lines"),   # 可以调大，1~1.5都好看
        legend.text = element_text(size = 12),
        axis.title.x = element_text(size = 12))
###################################################多因素分析
library(logistf)
library(dplyr)
library(tibble)

# 待分析变量
vars <- c("IL10","CXCL8","NGF","IL13","CCL7","IL3","TNF","IL17A")

# 创建空表
res2 <- tibble(
  Variable = character(),
  OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric()
)

# 构建公式（固定混杂变量）
formula <- as.formula(
  paste("LC1 ~", paste(vars, collapse = " + "),"+ treattype")
)

fit <- logistf(formula, data = data)

# 提取每个变量的 OR、CI 和 p
for(v in vars){
  est <- fit$coefficients[v]
  ci <- confint(fit)[v, ]
  p <- fit$prob[v]
  
  res2 <- res2 %>% add_row(
    Variable = v,
    OR = exp(est),
    CI_lower = exp(ci[1]),
    CI_upper = exp(ci[2]),
    p_value = p
  )
}

# 按 p 值排序
res2 <- res %>% arrange(p_value)
res2
###############################
##########################自身抗体B簇两个组之间比较
library(gtsummary)
library(dplyr)
data <- read.csv("LC84AAb90-clusterB.csv")
data <- data[3:68]
tbl_summary(data) #缺点是默认连续型变量不符合正态分布，统一采用非参数检验，分类变量统一采用皮尔逊卡方检验；优点是考虑缺失值并单独列出来
tbl_summary(data, by = LC1) # 根据LCC分层
res1 <- tbl_summary(data, by = LC1) %>% add_p()   #添加p值
#write.csv(res1, "autoantibodyB_LN.csv", row.names = FALSE)
##########################################
################B簇细胞因子之间没差异
library(logistf)
library(dplyr)
library(tibble)
library(ggplot2)

data <- read.csv("LC84AAb90-clusterB.csv")
data <- data[4:47]

data$LC1 <- as.factor(data$LC1)

bad_vars <- names(data)[sapply(data, function(x) n_distinct(x) < 2)]
data <- data %>% select(-all_of(bad_vars))

vars <- setdiff(names(data), c("LC1"))

# 初始化结果表
res <- tibble(
  Variable = character(),
  OR = numeric(),
  CI_lower = numeric(),
  CI_upper = numeric(),
  p_value = numeric()
)

# 逻辑回归循环（未调整treattype）
for (v in vars) {
  formula <- as.formula(paste("LC1 ~", v))
  
  fit <- tryCatch({
    logistf(formula, data = data)
  }, error = function(e) NULL)
  
  if(!is.null(fit)) {
    est <- fit$coefficients[v]
    ci <- confint(fit)[v, ]
    p <- fit$prob[v]
    
    res <- res %>% add_row(
      Variable = v,
      OR = exp(est),
      CI_lower = exp(ci[1]),
      CI_upper = exp(ci[2]),
      p_value = p
    )
  }
}

# 所有因子绘图，不筛选p
res_plot <- res %>%
  filter(!is.na(OR)) %>%
  mutate(
    Variable = factor(Variable, levels = rev(unique(Variable))),
    sig = ifelse(p_value < 0.05, "p<0.05", "ns")  # p>0.05用灰色
  )

# 按 OR 从小到大排序
#res_plot <- res_plot %>%
#  arrange(OR) %>%
#  mutate(Variable = factor(Variable, levels = Variable))  # 保持顺序

# 绘制森林图
ggplot(res_plot, aes(x = OR, y = Variable, color = sig)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = CI_lower, xmax = CI_upper), height = 0.2) +
  geom_vline(xintercept = 1, linetype = "dashed") +
  scale_x_log10() +
  scale_color_manual(values = c("p<0.05" = "red", "ns" = "grey50"),
                     labels = "Not\nsignificant") +
  theme_bw() +
  labs(x = "Odds Ratio (log scale)", y = "", color = "") +
  theme(
    axis.text.y = element_text(size = 11,face = "bold",color = "black"),
    axis.text.x = element_text(size = 11,face = "bold",color = "black"),
    legend.text = element_text(size = 11),
    axis.title.x = element_text(size = 11))
##################################ClusterA的ROC
#---------------------------------------
# 加载包
#---------------------------------------
library(glmnet)
library(pROC)
library(dplyr)

#---------------------------------------
# 读取数据
#---------------------------------------
df <- read.csv("LC84AAb90-clusterA.csv")
group <- df[[4]]         # LC1 分组 (0/1)
expr <- df[, 5:47]       # 表达矩阵

# 8个指定特征
features <- c("IL10","CXCL8","NGF","IL13","CCL7","IL3","TNF","IL17A")
x <- as.matrix(expr[, features])
y <- as.factor(group)

#---------------------------------------
# LASSO 逻辑回归 (对数几率模型)
#---------------------------------------
set.seed(1000)

cvfit <- cv.glmnet(
  x, y,
  family = "binomial",
  alpha = 1,         # LASSO
  nfolds = 5         # 小样本不宜折数太大
)

# 查看最优 lambda
cvfit$lambda.min
cvfit$lambda.1se

#---------------------------------------
# 提取最优模型下的非零变量（即最终选择的特征）
#---------------------------------------
coef_min <- coef(cvfit, s = "lambda.min")
selected <- rownames(coef_min)[coef_min[,1] != 0]
selected <- selected[selected != "(Intercept)"]
selected

cat("LASSO 选择的最佳特征：\n")
print(selected)
#############################
#---------------------------------------
# 系数柱状图
#---------------------------------------
library(ggplot2)

# 提取 lambda.min 下所有特征系数（去掉 Intercept）
coef_df <- data.frame(
  feature = rownames(coef_min)[-1], 
  coef = as.numeric(coef_min[-1,1])
)

# 排序（绝对值从大到小）
coef_df <- coef_df[order(abs(coef_df$coef), decreasing = TRUE), ]
coef_df$feature <- factor(coef_df$feature, levels=coef_df$feature)

# 绘制柱状图
ggplot(coef_df, aes(x=feature, y=coef, fill=coef>0)) +
  geom_bar(stat="identity") +
  coord_flip() +
  scale_fill_manual(values=c("TRUE"="#D73027","FALSE"="#4575B4"),
                    labels=c("Negative","Positive")) +
  theme_bw(base_size = 14) +
  labs(x="Feature", y="LASSO Coefficient", fill="Direction",
       title="LASSO Coefficients for 8 Features") +
  theme(legend.position = "right",
        axis.text.y = element_text(color="black") )

#---------------------------------------
# 计算最终模型的预测值（使用 lambda.min）
#---------------------------------------
library(ggplot2)
library(pROC)
library(dplyr)
library(glmnet)

# 选择特征
features_full <- c("NGF","CXCL8","IL3","IL10","CCL7")
features_single <- c("NGF")
y <- as.factor(group)

# -----------------------
# 多特征组合（5个）
# -----------------------
x_full <- as.matrix(expr[, features_full])
cvfit <- cv.glmnet(x_full, y, family="binomial", alpha=1, nfolds=5)
pred_full <- predict(cvfit, newx = x_full, s="lambda.min", type="response")

roc_full <- roc(y, as.numeric(pred_full))
roc_df_full <- data.frame(
  specificity = rev(roc_full$specificities),
  sensitivity = rev(roc_full$sensitivities)
)
auc_text_full <- paste0("5 factors AUC = ", round(auc(roc_full),3))

# -----------------------
# 单特征 NGF
# -----------------------
x_single <- as.matrix(expr[, features_single])
fit_single <- glm(y ~ x_single, family="binomial")
pred_single <- predict(fit_single, type="response")

roc_single <- roc(y, as.numeric(pred_single))
roc_df_single <- data.frame(
  specificity = rev(roc_single$specificities),
  sensitivity = rev(roc_single$sensitivities)
)
auc_text_single <- paste0("NGF only AUC = ", round(auc(roc_single),3))

# -----------------------
# 绘图
# -----------------------
ggplot() +
  # 多特征ROC曲线
  geom_line(data = roc_df_full,
            aes(x = 1 - specificity, y = sensitivity),
            color = "#BB0000", size = 1.2) +
  # 单特征ROC曲线
  geom_line(data = roc_df_single,
            aes(x = 1 - specificity, y = sensitivity),
            color = "#0072B2", size = 1.2) +
  # 对角线
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "grey50") +
  labs(
    x = "False-positive Rate",
    y = "True-positive Rate",
    title = "ROC Curves"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line    = element_line(colour = "black"),
    axis.text.y = element_text(colour = "black")
  ) +
  # 添加 AUC 文本到右下角
  annotate("text", x = 0.95, y = 0.05, label = auc_text_full,
           hjust = 1, vjust = 0, size = 5, color = "#BB0000") +
  annotate("text", x = 0.95, y = 0.15, label = auc_text_single,
           hjust = 1, vjust = 0, size = 5, color = "#0072B2")
###############################################################两簇自身抗体比较
##############################################################
#################################################################ALNC比较
library(dplyr)
library(tidyr)
library(reshape2)

# -------------------------------
# 读取簇信息
cluster_df <- read.csv("LN_candidate_cluster_assignments.csv", stringsAsFactors = FALSE)

# 读取阳性矩阵和mlxp
reads <- read.csv("reads.csv", row.names=1, check.names = FALSE)
mlxp  <- read.csv("mlxp.csv", row.names=1, check.names = FALSE)

# 构建阳性矩阵
positive_matrix <- (reads > 100 & mlxp > 1.3) * 1

# -------------------------------
# 样本分组
clu_A_samples <- cluster_df$Sample[cluster_df$Cluster == "A"]
LN_samples <- clu_A_samples[substr(clu_A_samples,1,1) %in% c("L","N")]

# C 样本：列名开头是 "C"
all_samples_reads <- colnames(reads)
C_samples <- all_samples_reads[substr(all_samples_reads,1,1) == "C"]

# 样本分组
all_samples <- c(LN_samples, C_samples)
sample_groups_all <- c(substr(LN_samples,1,1), rep("C", length(C_samples)))
names(sample_groups_all) <- all_samples
# -------------------------------
# 提取阳性矩阵对应样本
clu_pos_all <- positive_matrix[, all_samples]

# melt 矩阵
clu_pos_df_all <- as.data.frame(clu_pos_all)
clu_pos_df_all$Peptide <- rownames(clu_pos_df_all)
clu_long_all <- melt(clu_pos_df_all, id.vars="Peptide",
                     variable.name="Sample", value.name="Positive")
clu_long_all$Group <- sample_groups_all[clu_long_all$Sample]

# -------------------------------
# 读取 Apeptides.csv 中的肽
apeptides <- read.csv("Apeptides.csv", stringsAsFactors = FALSE)
cand_peptides <- apeptides$Peptide  # 假设第一列是 Peptide

# -------------------------------
# 只保留候选肽
clu_long_all_sel <- clu_long_all %>% filter(Peptide %in% cand_peptides)

# -------------------------------
# 计算每个肽在 L/N/C 的阳性率和总数
positive_rate_all <- clu_long_all_sel %>%
  group_by(Peptide, Group) %>%
  summarise(
    PositiveRate = mean(Positive),
    PositiveN = sum(Positive),
    TotalN = n(),
    .groups="drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(PositiveRate, PositiveN, TotalN),
    values_fill = 0
  )

# -------------------------------
# 多组阳性率比较（Chi-squared / Fisher）
fisher_res_all <- clu_long_all_sel %>%
  group_by(Peptide) %>%
  summarise(
    p_value = {
      tab <- table(Group, Positive)
      if(all(dim(tab)==2)) fisher.test(tab)$p.value
      else chisq.test(tab)$p.value
    },
    .groups="drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method="fdr"))

# -------------------------------
# 合并结果
final_table_C <- positive_rate_all %>%
  left_join(fisher_res_all, by="Peptide")

# -------------------------------
# 保存 CSV
write.csv(final_table_C, "clusterA_LNC_Apeptides_diff_peptides.csv", row.names=FALSE)
cat("✅ 簇 A L/N/C 对 Apeptides.csv 肽分析完成：clusterA_LNC_Apeptides_diff_peptides.csv\n")
###########################################################3
library(tidyverse)

df <- read.csv("A2.csv")

# 红色 FDR 显著性 (基于 p_value)
df <- df %>%
  mutate(sig_red = case_when(
    p_value < 0.001 ~ "***",
    p_value < 0.01  ~ "**",
    p_value < 0.05  ~ "*",
    TRUE ~ ""
  ))

# 蓝色原始 pvalue 显著性
df <- df %>%
  mutate(sig_blue = case_when(
    pvalue < 0.001 ~ "***",
    pvalue < 0.01  ~ "**",
    pvalue < 0.05  ~ "*",
    TRUE ~ ""
  ))

# 排序
order_vec <- df %>%
  mutate(delta = PositiveRate_L - PositiveRate_N) %>%
  arrange(delta) %>%
  pull(Peptide)

df_long <- df %>%
  pivot_longer(
    cols = c(PositiveRate_L, PositiveRate_N, PositiveRate_C),
    names_to = "Group",
    values_to = "PositiveRate"
  ) %>%
  mutate(
    Group = recode(Group,
                   "PositiveRate_L" = "LC",
                   "PositiveRate_N" = "NLC",
                   "PositiveRate_C" = "HC"),
    Peptide = factor(Peptide, levels = order_vec)
  )

p <- ggplot(df_long, aes(x = PositiveRate, y = Peptide)) +
  geom_line(aes(group = Peptide), color = "grey70", linewidth = 0.5) +
  geom_point(aes(color = Group), size = 3) +
  # ███ 让点变大
  geom_point(aes(color = Group), size = 4.5) +
  # 🔴 红色星号 —— 放在 NLC 点右边
  geom_text(
    data = df %>% mutate(Peptide = factor(Peptide, levels = order_vec)),
    aes(x = PositiveRate_N + 0.05, y = Peptide, label = sig_red),
    color = "red", size = 4
  ) +
  
  # 🔵 蓝色星号 —— 放在 LC 点左边
  geom_text(
    data = df %>% mutate(Peptide = factor(Peptide, levels = order_vec)),
    aes(x = PositiveRate_N - 0.05, y = Peptide, label = sig_blue),
    color = "blue", size = 4, hjust = 1
  ) +
  
  scale_color_manual(values = c(
    "LC" = "#E64B35",
    "NLC" = "#4DBBD5",
    "HC" = "#00A087"
  )) +
  theme_classic(base_size = 14) +
  labs(
    x = "Positive rate",
    y = "",
    color = ""
  ) +
  theme(
    legend.position = "top",
    axis.text = element_text(color = "black")
  )

p
#########################################A富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 2. 准备你的差异基因列表（Symbol）
gene_df <- read.csv("gene.csv", header = TRUE, stringsAsFactors = FALSE)
genes <- na.omit(gene_df[[1]])
gene_list <- genes

# 3. ID 转换：将基因 Symbol 转换为 Entrez Gene ID
# org.Hs.eg.db 是人类基因注释数据库
# 从 'SYMBOL' 列映射到 'ENTREZID' 列
entrez_ids <- bitr(geneID = gene_list,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# 4. 进行 GO 和 KEGG 富集分析

# GO - 生物过程 (Biological Process, BP)
ego_bp <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)

# GO - 细胞组分 (Cellular Component, CC)
ego_cc <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05
)

# GO - 分子功能 (Molecular Function, MF)
ego_mf <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)

# KEGG 富集分析
options(timeout = 300)
ekegg <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                    organism     = 'hsa', # 人类 KEGG 库
                    pvalueCutoff = 0.05)

# 5. 导出所有 GO 和 KEGG 富集结果为 CSV 文件

# 导出 GO-BP 结果
if (!is.null(ego_bp)) {
  ego_bp_df <- as.data.frame(ego_bp@result)
  write.csv(ego_bp_df, "GO_BP_enrichment_results.csv", row.names = FALSE)
  print("GO生物过程富集结果已导出至 'GO_BP_enrichment_results.csv'")
} else {
  print("没有显著的GO生物过程富集通路，未导出文件。")
}

# 导出 GO-CC 结果
if (!is.null(ego_cc)) {
  ego_cc_df <- as.data.frame(ego_cc@result)
  write.csv(ego_cc_df, "GO_CC_enrichment_results.csv", row.names = FALSE)
  print("GO细胞组分富集结果已导出至 'GO_CC_enrichment_results.csv'")
} else {
  print("没有显著的GO细胞组分富集通路，未导出文件。")
}

# 导出 GO-MF 结果
if (!is.null(ego_mf)) {
  ego_mf_df <- as.data.frame(ego_mf@result)
  write.csv(ego_mf_df, "GO_MF_enrichment_results.csv", row.names = FALSE)
  print("GO分子功能富集结果已导出至 'GO_MF_enrichment_results.csv'")
} else {
  print("没有显著的GO分子功能富集通路，未导出文件。")
}

# 导出 KEGG 结果
#if (!is.null(ekegg)) {
#  ekegg_df <- as.data.frame(ekegg)
#  write.csv(ekegg_df, "KEGG_enrichment_results.csv", row.names = FALSE)
#  print("KEGG富集结果已导出至 'KEGG_enrichment_results.csv'")
#} else {
#  print("没有显著的KEGG富集通路，未导出文件。")
#}
if (!is.null(ekegg) && nrow(as.data.frame(ekegg@result)) > 0) {
  ekegg_df <- as.data.frame(ekegg@result)
  write.csv(ekegg_df, "KEGG_enrichment_results.csv", row.names = FALSE)
  print("KEGG富集结果已导出至 'KEGG_enrichment_results.csv'")
} else {
  print("KEGG富集结果为空，未导出文件。")
}
#############################################小分类KEGG
library(ggplot2)
library(ggforce)
library(readr)
library(dplyr)

# 读取 KEGG 数据
df <- read_csv("kegg.csv")

# 假设 df 中有 Description, pvalue, Count, Category 列
# Category 替换为你自己的分类列表，例如：
# df$Category <- c("Infectious disease: viral", "Cell growth and death", ...)

# 定义分类颜色（示例，你可以换成自己喜欢的调色板）
category_colors <- c(
  "Infectious disease: viral" = "#1F78B4",
  "Cell growth and death" = "#33A02C",
  "Signal transduction" = "#E31A1C",
  "Aging" = "#FF7F00",
  "Cancer: specific types" = "#6A3D9A",
  "Infectious disease: bacterial" = "#B15928",
  "Drug resistance: antineoplastic" = "#FDBF6F",
  "Endocrine system" = "#CAB2D6",
  "Nervous system" = "#B2DF8A",
  "Neurodegenerative disease" = "#FB9A99",
  "Cellular community - eukaryotes" = "#A6CEE3"
)

# 排序（按 pvalue 从小到大）
df <- df %>%
  arrange(pvalue) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# 绘制流星图
p <- ggplot(df) +
  geom_link(aes(x = 0, y = Description,
                xend = -log10(pvalue), yend = Description,
                color = Category,
                linewidth = after_stat(index)),
            n = 500, show.legend = c(color = TRUE, alpha = FALSE, linewidth = FALSE)) +
  geom_point(aes(x = -log10(pvalue), y = Description, fill = Category),
             color = "black",        # 黑边
             size = 6, shape = 21) + # 流星头
  geom_text(aes(x = -log10(pvalue), y = Description, label = Count),
            size = 3) +
  theme_classic() +
  theme(
    panel.grid = element_blank(),
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    axis.text.x = element_text(face = "bold", size = 11, color = "black"),
    axis.title = element_text(size = 12),
    legend.key.size = unit(1, "lines"),
    legend.text = element_text(size = 9, face = "bold"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.position = c(0.8, 0.7)  # legend 位置调整
  ) +
  xlab("-Log10 Pvalue") + 
  ylab("") +
  # 图例方块和颜色一致
  scale_fill_manual(values = category_colors,
                    guide = guide_legend(override.aes = list(shape = 22, size = 6))) +
  scale_color_manual(values = category_colors)

# 保存结果
ggsave(filename = "KEGG_cometplot.tiff",
       width = 10, height = 7, plot = p)
###################################################B总体比较
###############################################BLNC组间比较
library(dplyr)
library(tidyr)
library(reshape2)

# -------------------------------
# 读取簇信息
cluster_df <- read.csv("LN_candidate_cluster_assignments.csv", stringsAsFactors = FALSE)

# 读取阳性矩阵和mlxp
reads <- read.csv("reads.csv", row.names=1, check.names = FALSE)
mlxp  <- read.csv("mlxp.csv", row.names=1, check.names = FALSE)

# 构建阳性矩阵
positive_matrix <- (reads > 100 & mlxp > 1.3) * 1

# -------------------------------
# 样本分组
clu_A_samples <- cluster_df$Sample[cluster_df$Cluster == "B"]
LN_samples <- clu_A_samples[substr(clu_A_samples,1,1) %in% c("L","N")]

# C 样本：列名开头是 "C"
all_samples_reads <- colnames(reads)
C_samples <- all_samples_reads[substr(all_samples_reads,1,1) == "C"]

# 样本分组
all_samples <- c(LN_samples, C_samples)
sample_groups_all <- c(substr(LN_samples,1,1), rep("C", length(C_samples)))
names(sample_groups_all) <- all_samples
# -------------------------------
# 提取阳性矩阵对应样本
clu_pos_all <- positive_matrix[, all_samples]

# melt 矩阵
clu_pos_df_all <- as.data.frame(clu_pos_all)
clu_pos_df_all$Peptide <- rownames(clu_pos_df_all)
clu_long_all <- melt(clu_pos_df_all, id.vars="Peptide",
                     variable.name="Sample", value.name="Positive")
clu_long_all$Group <- sample_groups_all[clu_long_all$Sample]

# -------------------------------
# 读取 Apeptides.csv 中的肽
apeptides <- read.csv("Bpeptides.csv", stringsAsFactors = FALSE)
cand_peptides <- apeptides$Peptide  # 假设第一列是 Peptide

# -------------------------------
# 只保留候选肽
clu_long_all_sel <- clu_long_all %>% filter(Peptide %in% cand_peptides)

# -------------------------------
# 计算每个肽在 L/N/C 的阳性率和总数
positive_rate_all <- clu_long_all_sel %>%
  group_by(Peptide, Group) %>%
  summarise(
    PositiveRate = mean(Positive),
    PositiveN = sum(Positive),
    TotalN = n(),
    .groups="drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(PositiveRate, PositiveN, TotalN),
    values_fill = 0
  )

# -------------------------------
# 多组阳性率比较（Chi-squared / Fisher）
fisher_res_all <- clu_long_all_sel %>%
  group_by(Peptide) %>%
  summarise(
    p_value = {
      tab <- table(Group, Positive)
      if(all(dim(tab)==2)) fisher.test(tab)$p.value
      else chisq.test(tab)$p.value
    },
    .groups="drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method="fdr"))

# -------------------------------
# 合并结果
final_table_C <- positive_rate_all %>%
  left_join(fisher_res_all, by="Peptide")

# -------------------------------
# 保存 CSV
write.csv(final_table_C, "clusterB_LNC_Bpeptides_diff_peptides.csv", row.names=FALSE)
cat("✅ 簇 A/N/C 对 Bpeptides.csv 肽分析完成：clusterB_LNC_Bpeptides_diff_peptides.csv\n")
###########################################################3
library(ComplexHeatmap)
library(circlize)
library(dplyr)
library(readr)

# -------------------------------------------------------
# 你已有的对象：
# positive_matrix
# sample_groups_all
# cand_peptides
# -------------------------------------------------------

# 重新组织分组（C = HC, L = LC, N = NLC）
group_df <- data.frame(
  Sample = names(sample_groups_all),
  Group  = sample_groups_all
)

# -------------------------------------------------------
# positive matrix 只保留候选肽
pos_mat <- positive_matrix[cand_peptides, names(sample_groups_all)]
pos_mat <- as.matrix(pos_mat)

# -------------------------------------------------------
# 加载 zscore
zscore <- read_csv("zscore.csv")
z_mat <- as.matrix(zscore[, -1])      # 如果第一列是peptide名
rownames(z_mat) <- zscore[[1]]

# 同样取候选肽与样本
z_mat <- z_mat[cand_peptides, names(sample_groups_all)]

# -------------------------------------------------------
# 排列顺序：C → L → N
col_order <- c(
  group_df$Sample[group_df$Group == "C"],
  group_df$Sample[group_df$Group == "L"],
  group_df$Sample[group_df$Group == "N"]
)

pos_mat <- pos_mat[, col_order]
z_mat   <- z_mat[, col_order]
group_df <- group_df[match(col_order, group_df$Sample), ]

# -------------------------------------------------------
# 调整 group_df 名称，使用 HC / LC / NLC
group_df$GroupLabel <- recode(group_df$Group,
                              "C" = "HC",
                              "L" = "LC",
                              "N" = "NLC")
# Positive 热图的列注释，不显示名称
col_ha_pos <- HeatmapAnnotation(
  Group = group_df$GroupLabel,
  col = list(Group = c(
    "HC"  = "#00A087",  # HC
    "LC"  = "#E64B35",  # LC
    "NLC" = "#4DBBD5"   # NLC
  )),
  show_annotation_name = FALSE   # 隐藏列注释名
)


# Z-score 热图的列注释，显示名称
col_ha_z <- HeatmapAnnotation(
  Group = group_df$GroupLabel,
  col = list(Group = c(
    "HC"  = "#00A087",  # HC
    "LC"  = "#E64B35",  # LC
    "NLC" = "#4DBBD5"   # NLC
  )),
  show_annotation_name = TRUE   # ✅ 显示注释名称
)

# Positive 热图
ht_pos <- Heatmap(
  pos_mat,
  name = "Positive",
  col = c("0"="#e2f6ff","1"="#57adde"),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  top_annotation = col_ha_pos,
  show_column_names = FALSE,
  show_row_names = FALSE
)

# Z-score 热图
ht_zscore <- Heatmap(
  z_mat,
  name = "Z-score",
  col = colorRamp2(c(-2, 0, 2), c("#61AACF", "white", "#D6604D")),
  cluster_rows = TRUE,
  cluster_columns = FALSE,
  top_annotation = col_ha_z,
  show_column_names = FALSE,
  show_row_names = FALSE
)

# 拼接两个热图
draw(ht_pos + ht_zscore)
################################################BLN组间比较+RFE筛选
library(dplyr)
library(reshape2)
library(tibble)
library(tidyr)

# -------------------------------
# 读取簇信息
cluster_df <- read.csv("LN_candidate_cluster_assignments.csv", stringsAsFactors = FALSE)

# 读取阳性矩阵和mlxp
reads <- read.csv("reads.csv", row.names=1, check.names = FALSE)
mlxp  <- read.csv("mlxp.csv", row.names=1, check.names = FALSE)

# 阳性矩阵
positive_matrix <- (reads > 100 & mlxp > 1.3) * 1

# 读取zscore
zscore <- read.csv("zscore.csv", row.names = 1, check.names = FALSE)
# -------------------------------
# 1️⃣ 取簇B的样本
clu_B_samples <- cluster_df$Sample[cluster_df$Cluster == "B"]
LN_samples <- clu_B_samples[substr(clu_B_samples,1,1) %in% c("L","N")]
sample_groups <- substr(LN_samples, 1, 1)
names(sample_groups) <- LN_samples

# -------------------------------
# 2️⃣ 提取阳性矩阵对应样本
clu_pos <- positive_matrix[, LN_samples]

# -------------------------------
# 3️⃣ melt 矩阵，计算阳性率
clu_pos_df <- as.data.frame(clu_pos)
clu_pos_df$Peptide <- rownames(clu_pos_df)

clu_long <- melt(clu_pos_df, id.vars = "Peptide",
                 variable.name = "Sample", value.name = "Positive")
clu_long$Group <- sample_groups[clu_long$Sample]

positive_rate <- clu_long %>%
  group_by(Peptide, Group) %>%
  summarise(
    PositiveRate = mean(Positive),
    PositiveN = sum(Positive),
    TotalN = n(),
    .groups = "drop"
  )

# -------------------------------
# 4️⃣ pivot_wider 整合成一行一个肽
positive_rate_wide <- positive_rate %>%
  pivot_wider(
    names_from = Group,
    values_from = c(PositiveRate, PositiveN, TotalN),
    values_fill = 0
  )

# -------------------------------
# 5️⃣ 筛选 L 或 N 阳性率 > 0.1
candidate_peptides <- positive_rate_wide %>%
  filter(PositiveRate_L > 0.1 | PositiveRate_N > 0.1)

# -------------------------------
# 6️⃣ 保存 CSV
#write.csv(candidate_peptides, "clusterB_LN_candidate_peptides_oneRow.csv", row.names = FALSE)
cat("✅ 簇 B候选肽已保存为一行一个肽：clusterB_LN_candidate_peptides_oneRow.csv\n")
#################################################
library(dplyr)
library(reshape2)
library(tibble)
library(tidyr)

# -------------------------------
# 取候选肽
cand_peptides <- candidate_peptides$Peptide

# -------------------------------
# 1️⃣ 阳性率比较 (Fisher)
clu_long_sel <- clu_long %>% filter(Peptide %in% cand_peptides, Group %in% c("L","N"))

positive_LN <- clu_long_sel %>%
  group_by(Peptide, Group) %>%
  summarise(
    PositiveN = sum(Positive),
    TotalN = n(),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Group, values_from = c(PositiveN, TotalN), values_fill = 0)

fisher_res <- positive_LN %>%
  rowwise() %>%
  do({
    pos_L <- .$PositiveN_L; pos_N <- .$PositiveN_N
    neg_L <- .$TotalN_L - pos_L; neg_N <- .$TotalN_N - pos_N
    if(all(c(pos_L,pos_N,neg_L,neg_N) == 0)) return(data.frame(Peptide=.$Peptide, p_value_fisher=NA))
    tab <- rbind(c(pos_L,pos_N), c(neg_L,neg_N))
    test <- fisher.test(tab)
    data.frame(Peptide=.$Peptide, p_value_fisher=test$p.value)
  }) %>%
  ungroup() %>%
  mutate(FDR_fisher = p.adjust(p_value_fisher, method="fdr"))

# -------------------------------
# 2️⃣ Zscore 长格式
clu_z_sel <- zscore[cand_peptides, LN_samples]
clu_z_df <- as.data.frame(clu_z_sel)
clu_z_df$Peptide <- rownames(clu_z_df)

clu_z_long <- pivot_longer(clu_z_df, -Peptide, names_to="Sample", values_to="Zscore") %>%
  mutate(Group = substr(Sample,1,1)) %>%
  filter(Group %in% c("L","N"))

# -------------------------------
# 2️⃣ Zscore 平均值 + Welch t-test + FC（新增）
zscore_res <- clu_z_long %>%
  group_by(Peptide) %>%
  summarise(
    Z_mean_L = mean(Zscore[Group=="L"], na.rm=TRUE),
    Z_mean_N = mean(Zscore[Group=="N"], na.rm=TRUE),
    p_value_z = t.test(Zscore ~ Group, var.equal = FALSE)$p.value,
    .groups="drop"
  ) %>%
  mutate(
    FDR_z = p.adjust(p_value_z, method="fdr"),
    FC = ifelse(is.na(Z_mean_N) | Z_mean_N==0, NA, abs(Z_mean_L / Z_mean_N)),
    log2FC = log2(FC),
    abs_log2FC = abs(log2FC)
  )

# -------------------------------
# 3️⃣ 合并结果
final_table <- fisher_res %>%
  left_join(zscore_res, by="Peptide") %>%
  left_join(positive_rate_wide %>% filter(Peptide %in% cand_peptides), by="Peptide") %>%
  mutate(Cluster = "B")


# -------------------------------
# 4️⃣ 保存 CSV
write.csv(final_table, "clusterB_LN_diff_peptides1.csv", row.names = FALSE)
# -------------------------------
# 读取包含 FDR 信息的差异肽文件
diff_external <- read.csv("clusterB_LNC_Bpeptides_diff_peptides.csv")

# 取 FDR < 0.05 的肽段
peptides_FDR <- diff_external$Peptide[diff_external$FDR < 0.05]

# 5️⃣ 筛选差异肽：Fisher、Zscore 都显著，并且外部 FDR < 0.05
diff_peptides <- final_table %>%
  filter(
    p_value_fisher < 0.05,
    p_value_z < 0.05,
    Peptide %in% peptides_FDR   # ★新增条件
  )
##############################################################
library(mixOmics)
library(glmnet)
library(caret)
library(pROC)
library(dplyr)
library(ggplot2)

# ======================================================
# 1️⃣ 取簇 B 的 L/N 样本
# ======================================================
clu_B_samples <- cluster_df$Sample[cluster_df$Cluster == "B"]
LN_samples <- clu_B_samples[substr(clu_B_samples,1,1) %in% c("L","N")]
sample_groups <- substr(LN_samples, 1, 1)
names(sample_groups) <- LN_samples

# ======================================================
# 2️⃣ 筛选差异肽
# ======================================================
peptides2 <- diff_peptides$Peptide

# 提取Z-score矩阵：只取B簇的L/N样本和差异肽
zscore_sub <- zscore[rownames(zscore) %in% peptides2, LN_samples]
zscore_sub_t <- t(zscore_sub)

# 分组标签：L=1, N=0
labels <- ifelse(grepl("^L", rownames(zscore_sub_t)), 1,
                 ifelse(grepl("^N", rownames(zscore_sub_t)), 0, NA))
valid_idx <- !is.na(labels)
zscore_sub_t <- zscore_sub_t[valid_idx, , drop=FALSE]
labels <- labels[valid_idx]
table(labels)

# ======================================================
# 3️⃣ PLS-DA VIP筛选
# ======================================================
X <- zscore_sub_t
Y <- factor(labels, levels=c(0,1))
plsda_res <- mixOmics::plsda(X, Y, ncomp = 2)

vip_scores <- vip(plsda_res)
vip_selected <- rownames(vip_scores)[vip_scores[,1] > 1]
cat("✅ VIP>1 特征数:", length(vip_selected), "\n")

zscore_vip <- zscore_sub_t[, vip_selected, drop=FALSE]

# VIP可视化
vip_df <- data.frame(
  Feature = rownames(vip_scores),
  VIP = vip_scores[,1]
)
vip_df <- vip_df[order(vip_df$VIP, decreasing=FALSE), ]
vip_df$Feature <- factor(vip_df$Feature, levels = vip_df$Feature)

ggplot(vip_df, aes(x=Feature, y=VIP)) +
  geom_col(fill="#9fc7b4") +
  geom_hline(yintercept=1, color="red", linetype="dashed", size=1) +
  coord_flip() +
  labs(x="Feature", y="VIP Score", title="PLS-DA VIP Scores (Bcluster)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size=9, colour="black", face="bold"),
    plot.title = element_text(hjust=0.5)
  )

# ======================================================
# 4️⃣ LASSO + RFE
# ======================================================
features_current <- colnames(zscore_vip)
rfe_results <- data.frame(
  n_features = integer(),
  auc_mean = numeric(),
  auc_lower = numeric(),
  auc_upper = numeric()
)
rfe_features_all <- data.frame()

step_large <- 3
step_small <- 1
n_repeat_rfe <- 50
nfold_rfe <- 5

while(length(features_current) > 0){
  auc_values <- c()
  
  for(seed in 1:n_repeat_rfe){
    set.seed(seed)
    folds <- createFolds(labels, k = nfold_rfe)
    
    for(fold_idx in 1:nfold_rfe){
      test_idx <- folds[[fold_idx]]
      train_idx <- setdiff(1:length(labels), test_idx)
      if(length(unique(labels[test_idx])) < 2) next
      
      x_train <- as.matrix(zscore_vip[train_idx, features_current, drop=FALSE])
      y_train <- labels[train_idx]
      x_test <- as.matrix(zscore_vip[test_idx, features_current, drop=FALSE])
      y_test <- labels[test_idx]
      
      cvfit <- cv.glmnet(x_train, y_train, family="binomial", alpha=1, nfolds=5)
      pred_prob <- predict(cvfit, x_test, s="lambda.min", type="response")
      auc_values <- c(auc_values, auc(roc(y_test, as.numeric(pred_prob))))
    }
  }
  
  auc_mean <- mean(auc_values)
  auc_ci <- quantile(auc_values, probs = c(0.025, 0.975))
  nf <- length(features_current)
  
  cat("✅", nf, "features:", sprintf("%.3f [%.3f - %.3f]\n", auc_mean, auc_ci[1], auc_ci[2]))
  
  rfe_results <- rbind(rfe_results, data.frame(
    n_features = nf,
    auc_mean = auc_mean,
    auc_lower = auc_ci[1],
    auc_upper = auc_ci[2]
  ))
  
  fit_full <- cv.glmnet(as.matrix(zscore_vip[, features_current, drop=FALSE]),
                        labels, family="binomial", alpha=1, nfolds=5)
  coefs <- coef(fit_full, s="lambda.min")[-1,]
  coefs_abs <- abs(coefs)
  
  rfe_features_all <- rbind(rfe_features_all,
                            data.frame(
                              Round = paste0("n=", nf),
                              Feature = names(coefs),
                              Coefficient = as.numeric(coefs),
                              AbsCoefficient = coefs_abs
                            ))
  
  remove_n <- ifelse(nf > 15, step_large, step_small)
  to_remove <- names(sort(coefs_abs, decreasing = FALSE))[1:remove_n]
  features_current <- setdiff(features_current, to_remove)
  
  if(length(features_current) <= 1) break
}

write.csv(rfe_features_all, "RFE_all_rounds_features_Bcluster.csv", row.names = FALSE)
###########################################################
# 3️⃣ 可视化 RFE 结果
# ======================================================
library(ggplot2)

# 自定义横坐标刻度
breaks_large <- seq(max(rfe_results$n_features), 16, by = -5)  # 大于15，每5个
breaks_small <- 15:1                                       # 小于等于15，每1个
x_breaks <- c(breaks_large, breaks_small)

ggplot(rfe_results, aes(x=n_features, y=auc_mean)) +
  geom_line(color="#1f77b4", size=1) +
  geom_point(color="#1f77b4", size=2) +
  geom_ribbon(aes(ymin=auc_lower, ymax=auc_upper), alpha=0.2, fill="#1f77b4") +
  geom_vline(xintercept = 5, color="red", linetype="dashed", size=1) +  # 红色竖线
  scale_x_reverse(breaks = x_breaks) +
  labs(x="Number of Features", y="Mean CV AUC",
       title="LASSO-RFE Performance (VIP>1 → RFE)") +
  theme_bw()
# ============================================
# 导出最终前9个最重要特征名称
# 找到剩余特征数为9的轮次
# 去掉 "n="，转换为数字
rfe_features_all$RoundNum <- as.numeric(sub("n=", "", rfe_features_all$Round))

# 假设你要导出 n=11 的特征
round_5 <- rfe_features_all[rfe_features_all$RoundNum == 5, ]

top5_features <- round_5[, c("Feature", "Coefficient", "AbsCoefficient")]
#write.csv(top5_features, "RFE_top5_features.csv", row.names = FALSE)
df_top5 <- read.csv("Top5.csv", stringsAsFactors = FALSE)
#可视化
ggplot(df_top5, aes(x = AbsCoefficient, y = reorder(Feature, AbsCoefficient))) +
  geom_col(fill = "#1f77b4") +  # 改为蓝色
  labs(x = "Absolute LASSO Coefficient", y = "Feature",
       title = "Top 5 Features (LASSO)") +
  theme_bw() +
  theme(
    axis.text.y = element_text(size = 12, colour = "black", face = "bold"),  # 纵坐标文字变大黑色
    axis.text.x = element_text(size = 12, colour = "black"),
    axis.title = element_text(size = 13, face = "bold"),
    plot.title = element_text(hjust = 0.5, face = "bold")
  )

#####################################ROC
library(glmnet)
library(pROC)
library(ggplot2)
library(dplyr)
library(caret)

# =============================
# 1️⃣ 数据准备
# =============================
top6_features <- read.csv("RFE_top5_features.csv", stringsAsFactors = FALSE)
features <- top6_features$Feature
x <- as.matrix(zscore_vip[, features])
y <- labels

# =============================
# 2️⃣ 重复交叉验证预测
# =============================
nfold <- 5
n_repeat <- 50
set.seed(123)

pred_all_matrix <- matrix(NA, nrow=length(y), ncol=n_repeat)

for(seed in 1:n_repeat){
  set.seed(seed)
  folds <- createFolds(y, k = nfold)
  pred_all <- rep(NA, length(y))
  
  for(f in 1:nfold){
    test_idx <- folds[[f]]
    train_idx <- setdiff(1:length(y), test_idx)
    
    x_train <- x[train_idx, ]
    y_train <- y[train_idx]
    x_test <- x[test_idx, ]
    
    # LASSO逻辑回归
    fit <- cv.glmnet(x_train, y_train, family="binomial", alpha=1, nfolds=5)
    pred_prob <- predict(fit, x_test, s="lambda.min", type="response")
    
    pred_all[test_idx] <- pred_prob
  }
  pred_all_matrix[, seed] <- pred_all
}

# 平均预测概率
pred_mean <- rowMeans(pred_all_matrix)

# =============================
# 3️⃣ ROC + 置信区间计算
# =============================
roc_obj <- roc(y, pred_mean)
auc_val <- auc(roc_obj)

# 每个FPR的TPR置信区间（兼容各种 pROC 版本）
ci_obj <- ci.se(roc_obj, specificities = seq(0, 1, l = 50), boot.n = 1000, conf.level = 0.95)

# 将 ci.se 返回值转换为数据框
ci_df <- as.data.frame(ci_obj)
colnames(ci_df) <- c("sens_low", "sens_mid", "sens_high")
ci_df$specificity <- seq(0, 1, l = nrow(ci_df))

# ROC曲线数据
roc_df <- data.frame(
  specificity = rev(roc_obj$specificities),
  sensitivity = rev(roc_obj$sensitivities)
)
# AUC文本
auc_val <- auc(roc_obj)
auc_ci <- ci.auc(roc_obj, conf.level = 0.95)

auc_text <- paste0("AUC = ", round(auc_val, 3),
                   " [", round(auc_ci[1], 3), "-", round(auc_ci[3], 3), "]")
# =============================
# 4️⃣ 绘图
# =============================
ggplot() +
  # CI 阴影
  geom_ribbon(data = ci_df,
              aes(x = 1 - specificity, ymin = sens_low, ymax = sens_high),
              fill = "#BB0000", alpha = 0.2) +
  # ROC 曲线
  geom_line(data = roc_df,
            aes(x = 1 - specificity, y = sensitivity),
            color = "#BB0000", size = 1.2) +
  # 对角线
  geom_abline(intercept = 0, slope = 1,
              linetype = "dashed", color = "grey50") +
  labs(
    x = "False-positive Rate",
    y = "True-positive Rate",
    title = "Cross-validated ROC"
  ) +
  theme_classic(base_size = 14) +
  theme(
    plot.title   = element_text(hjust = 0.5, face = "bold"),
    panel.border = element_rect(colour = "black", fill = NA, linewidth = 1),
    axis.line    = element_line(colour = "black")
  ) +
  # 添加 AUC + CI 到右下角
  annotate("text", 
           x = 0.95, y = 0.05,
           label = auc_text, 
           hjust = 1, vjust = 0,
           size = 5, fontface = "bold", color = "black")
############################################################最终B重要肽阳性率比较
library(dplyr)
library(tidyr)
library(reshape2)
library(openxlsx)
# -------------------------------
# 读取簇信息
cluster_df <- read.csv("LN_candidate_cluster_assignments.csv", stringsAsFactors = FALSE)
# 读取阳性矩阵和mlxp
reads <- read.csv("reads.csv", row.names=1, check.names = FALSE)
mlxp  <- read.csv("mlxp.csv", row.names=1, check.names = FALSE)
# 构建阳性矩阵
positive_matrix <- (reads > 100 & mlxp > 1.3) * 1
# 读取 zscore
zscore <- read.csv("zscore.csv", row.names=1, check.names=FALSE)
# -------------------------------
# 样本分组
clu_B_samples <- cluster_df$Sample[cluster_df$Cluster == "B"]
LN_samples <- clu_B_samples[substr(clu_B_samples,1,1) %in% c("L","N")]
all_samples_reads <- colnames(reads)
C_samples <- all_samples_reads[substr(all_samples_reads,1,1) == "C"]
all_samples <- c(LN_samples, C_samples)
sample_groups_all <- c(substr(LN_samples,1,1), rep("C", length(C_samples)))
names(sample_groups_all) <- all_samples
# -------------------------------
# 提取候选肽对应的原始矩阵
apeptides <- read.csv("Bpeptides5.csv", stringsAsFactors = FALSE)
cand_peptides <- apeptides$Peptide

positive_sel <- positive_matrix[cand_peptides, all_samples]
zscore_sel   <- zscore[cand_peptides, all_samples]
# -------------------------------
# melt 数据用于统计分析
clu_pos_df_all <- as.data.frame(positive_sel)
clu_pos_df_all$Peptide <- rownames(clu_pos_df_all)
clu_long_all <- melt(clu_pos_df_all, id.vars="Peptide",
                     variable.name="Sample", value.name="Positive")
clu_long_all$Group <- sample_groups_all[clu_long_all$Sample]

zscore_df <- as.data.frame(zscore_sel)
zscore_df$Peptide <- rownames(zscore_df)
zscore_long <- melt(zscore_df, id.vars="Peptide",
                    variable.name="Sample", value.name="Zscore")
zscore_long$Group <- sample_groups_all[zscore_long$Sample]

# -------------------------------
# 计算阳性率和总数
positive_rate_all <- clu_long_all %>%
  group_by(Peptide, Group) %>%
  summarise(
    PositiveRate = mean(Positive),
    PositiveN = sum(Positive),
    TotalN = n(),
    .groups="drop"
  ) %>%
  pivot_wider(
    names_from = Group,
    values_from = c(PositiveRate, PositiveN, TotalN),
    values_fill = 0
  )
# 多组 Zscore 比较（ANOVA）
zscore_res_all <- zscore_long %>%
  group_by(Peptide) %>%
  summarise(
    Z_mean_L = mean(Zscore[Group=="L"], na.rm=TRUE),
    Z_mean_N = mean(Zscore[Group=="N"], na.rm=TRUE),
    Z_mean_C = mean(Zscore[Group=="C"], na.rm=TRUE),
    p_value_z = if(length(unique(Group))>1) summary(aov(Zscore ~ Group))[[1]][["Pr(>F)"]][1] else NA,
    .groups="drop"
  ) %>%
  mutate(FDR_z = p.adjust(p_value_z, method="fdr"))
# 多组阳性率比较（Chi-squared / Fisher）
fisher_res_all <- clu_long_all %>%
  group_by(Peptide) %>%
  summarise(
    p_value = {
      tab <- table(Group, Positive)
      if(all(dim(tab)==2)) fisher.test(tab)$p.value
      else chisq.test(tab)$p.value
    },
    .groups="drop"
  ) %>%
  mutate(FDR = p.adjust(p_value, method="fdr"))
# 合并汇总结果
final_table_C <- positive_rate_all %>%
  left_join(fisher_res_all, by="Peptide") %>%
  left_join(zscore_res_all, by="Peptide")
# -------------------------------
# 保存到 Excel
wb <- createWorkbook()
# 汇总表
addWorksheet(wb, "Summary")
writeData(wb, "Summary", final_table_C)
# 原始阳性/阴性矩阵
addWorksheet(wb, "PositiveMatrix")
writeData(wb, "PositiveMatrix", positive_sel, rowNames=TRUE)
# 原始 Zscore 矩阵
addWorksheet(wb, "ZscoreMatrix")
writeData(wb, "ZscoreMatrix", zscore_sel, rowNames=TRUE)
# 保存 Excel 文件
saveWorkbook(wb, "clusterB_Bpeptides_full_data.xlsx", overwrite=TRUE)
##############################################富集分析B上调
#################################富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)
# 2. 准备你的差异基因列表（Symbol）
gene_df <- read.csv("gene1.csv", header = FALSE, stringsAsFactors = FALSE)
genes <- na.omit(gene_df[[1]])
gene_list <- genes

# 3. ID 转换：将基因 Symbol 转换为 Entrez Gene ID
# org.Hs.eg.db 是人类基因注释数据库
# 从 'SYMBOL' 列映射到 'ENTREZID' 列
entrez_ids <- bitr(geneID = gene_list,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# 4. 进行 GO 和 KEGG 富集分析

# GO - 生物过程 (Biological Process, BP)
ego_bp <- enrichGO(
  gene          = entrez_ids$ENTREZID,
  OrgDb         = org.Hs.eg.db,
  ont           = "BP",
  keyType       = "ENTREZID",
  pAdjustMethod = "BH",
  pvalueCutoff  = 0.3,   # 放宽到0.13
  qvalueCutoff  = 0.17     # qvalue<0.1   #这里的pvalue在内部是p.adjust，这里是为了后面做树状图，要不都卡掉了，不影响真正的pvalue
)

# GO - 细胞组分 (Cellular Component, CC)
ego_cc <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05
)

# GO - 分子功能 (Molecular Function, MF)
ego_mf <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)

# KEGG 富集分析
options(timeout = 300)
ekegg <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                    organism     = 'hsa', # 人类 KEGG 库
                    pvalueCutoff = 0.05)

# 5. 导出所有 GO 和 KEGG 富集结果为 CSV 文件

# 导出 GO-BP 结果
if (!is.null(ego_bp)) {
  ego_bp_df <- as.data.frame(ego_bp@result)
  write.csv(ego_bp_df, "GO_BP_enrichment_results.csv", row.names = FALSE)
  print("GO生物过程富集结果已导出至 'GO_BP_enrichment_results.csv'")
} else {
  print("没有显著的GO生物过程富集通路，未导出文件。")
}

# 导出 GO-CC 结果
if (!is.null(ego_cc)) {
  ego_cc_df <- as.data.frame(ego_cc@result)
  write.csv(ego_cc_df, "GO_CC_enrichment_results.csv", row.names = FALSE)
  print("GO细胞组分富集结果已导出至 'GO_CC_enrichment_results.csv'")
} else {
  print("没有显著的GO细胞组分富集通路，未导出文件。")
}

# 导出 GO-MF 结果
if (!is.null(ego_mf)) {
  ego_mf_df <- as.data.frame(ego_mf@result)
  write.csv(ego_mf_df, "GO_MF_enrichment_results.csv", row.names = FALSE)
  print("GO分子功能富集结果已导出至 'GO_MF_enrichment_results.csv'")
} else {
  print("没有显著的GO分子功能富集通路，未导出文件。")
}
if (!is.null(ekegg) && nrow(as.data.frame(ekegg@result)) > 0) {
  ekegg_df <- as.data.frame(ekegg@result)
  write.csv(ekegg_df, "KEGG_enrichment_results.csv", row.names = FALSE)
  print("KEGG富集结果已导出至 'KEGG_enrichment_results.csv'")
} else {
  print("KEGG富集结果为空，未导出文件。")
}
#############################################################KEGG柱状图
library(ggplot2)
library(readr)
library(dplyr)

# 读取 KEGG 数据
df <- read_csv("kegg.csv")

# 定义分类颜色
colors <- c(
  "Endocrine system" = "#3C5488",
  "Chromosome" = "#66B3E6",
  "Sensory system" = "#E68B81",
  "Immune system" = "#84C9A7",
  "Environmental adaptation" = "#EAAA60"
)

# 排序（按显著性）
df <- df %>%
  arrange(pvalue) %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# 绘制横向柱状图
p <- ggplot(df, aes(x = Description, y = -log10(pvalue), fill = Category)) +
  geom_col(color = "black", width = 0.7) +
  geom_text(aes(label = Count), 
            hjust = -0.2, size = 3.2, fontface = "bold") +
  coord_flip() +   # 横向显示
  theme_classic() +
  theme(
    axis.line = element_line(color = "black", linewidth = 0.6),
    axis.text = element_text(face = "bold"),
    axis.text.y = element_text(face = "bold", size = 11, color = "black"),
    axis.text.x = element_text(face = "bold", size = 11, color = "black"),
    axis.title = element_text(size = 12),
    legend.key.size = unit(1, "lines"),
    legend.text = element_text(size = 9, face = "bold"),
    legend.title = element_text(size = 9, face = "bold"),
    legend.position = "right",
    legend.direction = "vertical"
  ) +
  ylab("-Log10 Pvalue") +
  xlab("") +
  scale_fill_manual(values = colors)

# 保存
ggsave("KEGG_barplot_horizontal.tiff", p, width = 9, height = 4)

##################################################
library(clusterProfiler)
library(enrichplot)
library(ggplot2)

# 假设 ego_bp 是 enrichResult 对象
# -------------------------------
# 1️⃣ 指定要删除的 term 列表
remove_terms <- c(
  "negative regulation of intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator",
  "regulation of intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator",
  "negative regulation of signal transduction by p53 class mediator",
  "regulation of intrinsic apoptotic signaling pathway by p53 class mediator",
  "regulation of intrinsic apoptotic signaling pathway in response to DNA damage",
  "intrinsic apoptotic signaling pathway in response to DNA damage by p53 class mediator"
)

# -------------------------------
# 2️⃣ 过滤掉这些 term 并去重
ego_bp_filtered <- ego_bp
ego_bp_filtered@result <- ego_bp_filtered@result[!(ego_bp_filtered@result$Description %in% remove_terms), ]
ego_bp_filtered@result <- ego_bp_filtered@result[!duplicated(ego_bp_filtered@result$Description), ]

# -------------------------------
# 3️⃣ 计算 term 相似性
ego_bp_filtered <- pairwise_termsim(ego_bp_filtered)

# -------------------------------
# 4️⃣ 绘图
p <- treeplot(ego_bp_filtered,
              showCategory = 30,
              color = "pvalue",
              label_format = NULL,
              fontsize = 4,
              hilight.params = list(hilight = TRUE, align = "both"),
              offset.params = list(bar_tree = rel(3), tiplab = rel(4), extend = 0.2, hexpand = 0.2),
              cluster.params = list(method = "ward.D",
                                    n = 5,
                                    color =  c("#F0E442","#F4A6B9","#E69F00", "#56B4E9","#009E73"),##控制颜色框的颜色，c("#999999", "#E69F00", "#56B4E9", "#009E73", "#F0E442")
                                    label_words_n = 6,
                                    label_format = 30))#

# -------------------------------
# 5️⃣ cluster label 加粗
p$layers[[3]]$aes_params$fontface <- "bold"

p
#######################################################GOCC和MF气泡
library(ggplot2)
library(dplyr)
library(forcats)

# 读取数据
go_df <- read.csv("GO.csv", stringsAsFactors = FALSE)

# 按显著性排序 GO 条目
go_df <- go_df %>%
  group_by(Category) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  ungroup()

# 将 GO 条目按照 -log10(pvalue) 排序
go_df <- go_df %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# 定义颜色
colors <- c(
  "Molecular Function" = "#ff7f0e",
  "Cellular Component" = "#2ca02c"
)

# 绘图
p <- ggplot(go_df, aes(x = -log10(pvalue), y = Description)) +
  geom_point(aes(size = Count, fill = Category), shape = 21, color = "black") +
  scale_fill_manual(values = colors) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw() +   # 四边框
  theme(
    axis.text.y = element_text(face = "bold", size = 13, color = "black"),  # 改成黑色
    axis.text.x = element_text(face = "bold", size = 13, color = "black"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold")
  ) +
  ylab("") +
  xlab("-Log10 P-value") +
  guides(
    fill = guide_legend(title = "Category", override.aes = list(size = 5))
  )

ggsave("GO_triplet_bubbleplot_colored_y1.tiff", plot = p, width = 10, height = 8)
###############################################富集分析B下调
#################################富集分析
library(clusterProfiler)
library(org.Hs.eg.db)
library(enrichplot)
library(ggplot2)

# 2. 准备你的差异基因列表（Symbol）
gene_df <- read.csv("gene2.csv", header = FALSE, stringsAsFactors = FALSE)
genes <- na.omit(gene_df[[1]])
gene_list <- genes

# 3. ID 转换：将基因 Symbol 转换为 Entrez Gene ID
# org.Hs.eg.db 是人类基因注释数据库
# 从 'SYMBOL' 列映射到 'ENTREZID' 列
entrez_ids <- bitr(geneID = gene_list,
                   fromType = "SYMBOL",
                   toType = "ENTREZID",
                   OrgDb = org.Hs.eg.db)

# 4. 进行 GO 和 KEGG 富集分析

# GO - 生物过程 (Biological Process, BP)
ego_bp <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "BP",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)

# GO - 细胞组分 (Cellular Component, CC)
ego_cc <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "CC",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05
)

# GO - 分子功能 (Molecular Function, MF)
ego_mf <- enrichGO(gene          = entrez_ids$ENTREZID,
                   OrgDb         = org.Hs.eg.db,
                   ont           = "MF",
                   pAdjustMethod = "BH",
                   pvalueCutoff  = 0.05)

# KEGG 富集分析
options(timeout = 300)
ekegg <- enrichKEGG(gene         = entrez_ids$ENTREZID,
                    organism     = 'hsa', # 人类 KEGG 库
                    pvalueCutoff = 0.05)

# 5. 导出所有 GO 和 KEGG 富集结果为 CSV 文件

# 导出 GO-BP 结果
if (!is.null(ego_bp)) {
  ego_bp_df <- as.data.frame(ego_bp@result)
  write.csv(ego_bp_df, "GO_BP_enrichment_results.csv", row.names = FALSE)
  print("GO生物过程富集结果已导出至 'GO_BP_enrichment_results.csv'")
} else {
  print("没有显著的GO生物过程富集通路，未导出文件。")
}

# 导出 GO-CC 结果
if (!is.null(ego_cc)) {
  ego_cc_df <- as.data.frame(ego_cc@result)
  write.csv(ego_cc_df, "GO_CC_enrichment_results.csv", row.names = FALSE)
  print("GO细胞组分富集结果已导出至 'GO_CC_enrichment_results.csv'")
} else {
  print("没有显著的GO细胞组分富集通路，未导出文件。")
}

# 导出 GO-MF 结果
if (!is.null(ego_mf)) {
  ego_mf_df <- as.data.frame(ego_mf@result)
  write.csv(ego_mf_df, "GO_MF_enrichment_results.csv", row.names = FALSE)
  print("GO分子功能富集结果已导出至 'GO_MF_enrichment_results.csv'")
} else {
  print("没有显著的GO分子功能富集通路，未导出文件。")
}

if (!is.null(ekegg) && nrow(as.data.frame(ekegg@result)) > 0) {
  ekegg_df <- as.data.frame(ekegg@result)
  write.csv(ekegg_df, "KEGG_enrichment_results.csv", row.names = FALSE)
  print("KEGG富集结果已导出至 'KEGG_enrichment_results.csv'")
} else {
  print("KEGG富集结果为空，未导出文件。")
}
########################################################GOKEGG气泡
library(ggplot2)
library(dplyr)
library(forcats)

# 读取数据
go_df <- read.csv("GOKEGG.csv", stringsAsFactors = FALSE)

# 按显著性排序 GO 条目
go_df <- go_df %>%
  group_by(Category) %>%
  arrange(pvalue, .by_group = TRUE) %>%
  ungroup()

# 将 GO 条目按照 -log10(pvalue) 排序
go_df <- go_df %>%
  mutate(Description = factor(Description, levels = unique(Description)))

# 定义颜色
colors <- c(
  "KEGG" = "#BA0000",
  "GO: BP" = "#ff7f0e",#ff7f0e
  "GO: MF" = "#1f77b4",
  "GO: CC" = "#2ca02c"
)

# 绘图
p <- ggplot(go_df, aes(x = -log10(pvalue), y = Description)) +
  geom_point(aes(size = Count, fill = Category), shape = 21, color = "black") +
  scale_fill_manual(values = colors) +
  scale_size_continuous(range = c(2, 8)) +
  theme_bw() +   # 四边框
  theme(
    axis.text.y = element_text(face = "bold", size = 13, color = "black"),  # 改成黑色
    axis.text.x = element_text(face = "bold", size = 13, color = "black"),
    axis.title = element_text(size = 12),
    legend.title = element_text(size = 12, face = "bold"),
    legend.text = element_text(size = 12, face = "bold")
  ) +
  ylab("") +
  xlab("-Log10 P-value") +
  guides(
    fill = guide_legend(title = "Category", override.aes = list(size = 5))
  )

ggsave("GO_triplet_bubbleplot_colored_y1.tiff", plot = p, width = 10, height = 8)
############################################序列比对 分子模拟
library(Biostrings)
library(dplyr)
library(msa)

# -----------------------------
# 读取肽段名称
# -----------------------------
pep_names <- read.csv("peptides3.csv", header = TRUE, stringsAsFactors = FALSE)
pep_names <- pep_names[[1]]  # 第一列是肽段名称

# -----------------------------
# 读取FASTA格式肽段序列
# -----------------------------
all_seqs <- readAAStringSet("protein_tiles_AAG_P4625_20250714_toZX.txt", format="fasta")
# -----------------------------
# 挑出对应肽段的序列
# -----------------------------
selected_seqs <- all_seqs[names(all_seqs) %in% pep_names]

# 检查是否全部找到
cat("找到", length(selected_seqs), "条肽段序列\n")
# -----------------------------
# 保存成FASTA文件
# -----------------------------
writeXStringSet(selected_seqs,
                filepath = "diff_peptides3.fasta",
                format   = "fasta")

# -----------------------------
# 多序列比对
# -----------------------------
alignment <- msa(selected_seqs, method = "ClustalW")  # 或 "Muscle"
print(alignment)

# -----------------------------
# 导出比对结果
# -----------------------------
# 安装 ggseqlogo
#if(!requireNamespace("ggseqlogo", quietly = TRUE)) install.packages("ggseqlogo")
library(ggseqlogo)

# 转成字符向量
seqs_char <- as.character(selected_seqs)

# 绘制序列 logo 并保存为 png
png("peptides_logo.png", width=1200, height=400)
ggseqlogo(seqs_char, method="prob")  # method="prob" 用频率绘制
dev.off()

###########################################################
library(Biostrings)
# 病毒肽段
virus <- readAAStringSet("virus.fasta")
# 人蛋白肽段
human <- readAAStringSet("diff_peptides3.fasta")
library(Biostrings)  # 更高效的序列匹配

# human: AAStringSet 短肽
# virus: AAStringSet 全长病毒蛋白
# min_overlap: 最小匹配长度

matches <- data.frame(VirusID=character(),
                      HumanID=character(),
                      MatchSeq=character(),
                      stringsAsFactors=FALSE)

min_overlap <- 5

for(h_idx in seq_along(human)) {
  h_seq <- human[h_idx]
  h_id <- names(human)[h_idx]
  
  for(v_idx in seq_along(virus)) {
    v_seq <- virus[v_idx]
    v_id <- names(virus)[v_idx]
    
    # 短肽滑动窗口
    for(start in 1:(nchar(h_seq)-min_overlap+1)) {
      for(len in min_overlap:(nchar(h_seq)-start+1)) {
        sub_h <- substr(h_seq, start, start+len-1)
        
        if(grepl(sub_h, v_seq, fixed=TRUE)){
          matches <- rbind(matches,
                           data.frame(VirusID=v_id,
                                      HumanID=h_id,
                                      MatchSeq=sub_h,
                                      stringsAsFactors=FALSE))
        }
      }
    }
  }
}
write.csv(matches, "peptide_matches3.csv", row.names = FALSE)
##########################################################组织定位
library(dplyr)
library(readr)
library(stringr)

# 1. 读取文件
proteins <- read.csv("IDBsig.csv", stringsAsFactors = FALSE)
ihc <- read_tsv("normal_tissue.tsv", col_types = cols())
ihc_issue <- read_tsv("normal_ihc_tissues.tsv", col_types = cols())

# 2. 过滤出目标蛋白的 IHC 记录（注意列名是否恰好为 `Gene name` 和 proteins$name）
ihc_subset <- ihc %>%
  filter(`Gene name` %in% proteins$name)

# 3. 规范 Tissue 字段（去两端空白、统一大小写）以提高匹配成功率
ihc_subset <- ihc_subset %>%
  mutate(Tissue = str_trim(Tissue))

ihc_issue <- ihc_issue %>%
  mutate(Tissue = str_trim(Tissue))

# 4. 按 Tissue 把 Organ 合并到 ihc_subset（以 ihc_subset 为主表）
ihc_merged <- ihc_subset %>%
  left_join(ihc_issue %>% select(Tissue, Organ), by = "Tissue")

# 5. 检查哪些 Tissue 未匹配到 Organ（便于排查）
unmatched_tissues <- ihc_merged %>%
  filter(is.na(Organ)) %>%
  distinct(Tissue) %>%
  arrange(Tissue)

cat("未匹配到 Organ 的 Tissue 数量：", nrow(unmatched_tissues), "\n")
if (nrow(unmatched_tissues) > 0) {
  print(unmatched_tissues)
}

# 6. 写出结果
write.csv(ihc_merged, "protein_tissue_with_organ_sig.csv", row.names = FALSE)

# 查看前几行确认
head(ihc_merged)
#######################################upset图
library(dplyr)
library(tidyr)
library(UpSetR)

# 假设你的数据叫 df
df <- read.csv("定位sig.csv",header = TRUE)
head(df)
# 1. 对每个基因收集其所有 Organ
gene_organ <- df %>%
  select(Gene.name, Organ) %>%
  distinct() %>%
  group_by(Gene.name, Organ) %>%
  summarise(value = 1, .groups = "drop")

# 2. 转成 UpSetR 需要的 wide format（行=基因，列=organ）
mat <- gene_organ %>%
  pivot_wider(names_from = Organ, values_from = value, values_fill = 0) %>%
  as.data.frame()

# 3. 第一列 Gene.name 设为 rownames
rownames(mat) <- mat$Gene.name
mat$Gene.name <- NULL

# 4. 绘制 UpSet 图（上半部分变矮）
upset(mat,
      nsets = ncol(mat),
      nintersects = 20,
      order.by = "freq",
      sets.bar.color = "grey",
      main.bar.color = "grey40",
      mb.ratio = c(0.55, 0.45)   # ← 控制上下布局比例
)
#######################################################A簇mofa
########################################################
# 0️⃣ 加载依赖包
library(MOFA2)
library(dplyr)
library(tibble)
library(readr)
library(reticulate)

########################################################
# 1️⃣ 读取数据
cytokine <- read.csv("ccytokine.csv", header = TRUE, row.names = 1, check.names = FALSE)
positive <- read.csv("cpositive.csv", header = TRUE, row.names = 1, check.names = FALSE)
zscore <- read.csv("czscore.csv", header = TRUE, row.names = 1, check.names = FALSE)

# 1️⃣b 读取分簇信息
cluster_info <- read.csv("LN_candidate_cluster_assignments.csv", header = TRUE, stringsAsFactors = FALSE)
# 只保留 A 簇样本
samples_A <- cluster_info %>% filter(Cluster == "A") %>% pull(Sample)

# 1️⃣c 保留共同样本 & 只用 A 簇
common_samples <- Reduce(intersect, list(colnames(cytokine), colnames(positive), colnames(zscore), samples_A))
cytokine <- cytokine[, common_samples]
positive <- positive[, common_samples]
zscore <- zscore[, common_samples]

########################################################
# 2️⃣ 数据类型处理
cytokine_matrix <- as.matrix(cytokine)
storage.mode(cytokine_matrix) <- "double"

zscore_matrix <- as.matrix(zscore)
storage.mode(zscore_matrix) <- "double"

positive_matrix <- as.matrix(positive)
positive_matrix <- ifelse(positive_matrix != 0, 1, 0)
storage.mode(positive_matrix) <- "integer"

########################################################
# 3️⃣ 构建 MOFA 输入列表
data_list <- list(
  Cytokine = cytokine_matrix,
  Peptide_Zscore = zscore_matrix,
  Peptide_Positive = positive_matrix
)

########################################################
# 4️⃣ 创建 MOFA 对象
mofa_obj <- create_mofa(data_list)

# 5️⃣ 设置参数
data_opts <- get_default_data_options(mofa_obj)
model_opts <- get_default_model_options(mofa_obj)
train_opts <- get_default_training_options(mofa_obj)

model_opts$num_factors <- 7
train_opts$convergence_mode <- "medium"
model_opts$likelihoods <- c(
  Cytokine = "gaussian",
  Peptide_Zscore = "gaussian",
  Peptide_Positive = "bernoulli"
)

mofa_obj <- prepare_mofa(
  object = mofa_obj,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)

########################################################
# 6️⃣ 训练 MOFA
use_python("E:/Miniforge3/python.exe", required = TRUE)
reticulate::py_install("mofapy2", pip = TRUE)

mofa_trained <- run_mofa(mofa_obj, use_basilisk = FALSE)

########################################################
# 7️⃣ 提取 factor 矩阵（只包含 A 簇样本）
factor_df <- get_factors(mofa_trained, factors = "all", as.data.frame = TRUE)
########################################################
library(gtsummary)
library(dplyr)
data <- read.csv("clinical_LC.csv")
data <- data[,2:9]
tbl_summary(data) #缺点是默认连续型变量不符合正态分布，统一采用非参数检验，分类变量统一采用皮尔逊卡方检验；优点是考虑缺失值并单独列出来
tbl_summary(data, by = LC) # 根据LCC分层
tbl_summary(data, by = LC) %>% add_p()   #添加p值
########################################################
# 提取特征权重
weights_df <- get_weights(mofa_trained, views = "all", factors = "all", as.data.frame = TRUE)
# 每个 factor 下绝对值前 20 个特征，不管属于哪个 view
top_features <- weights_df %>%
  group_by(factor) %>%
  slice_max(order_by = abs(value), n = 30, with_ties = FALSE)  # with_ties=FALSE 保证正好 20 个
###########################################先改名字
library(dplyr)
library(stringr)
library(readr)
# 读取 name.csv
# 第一列 Index 对应 feature 中的 AAG编号，第二列 Name 是要替换的名称
name_map <- read.csv("name.csv", header = TRUE, stringsAsFactors = FALSE)
# 处理 top_features
top_features <- top_features %>%
  # 去掉 _Peptide_Zscore
  mutate(feature = str_remove(feature, "_Peptide_Zscore")) %>%
  # 提取编号部分和剩余部分
  mutate(
    AAG_id = str_extract(feature, "AAG\\d+"),
    rest = str_extract(feature, "\\|.*")  # |后面的部分，包括|
  ) %>%
  # 用 name.csv 替换 AAG编号
  left_join(name_map, by = c("AAG_id" = "Index")) %>%
  mutate(
    feature = ifelse(!is.na(Name), paste0(Name, rest), feature)
  ) %>%
  select(-AAG_id, -rest, -Name)
# 查看结果
head(top_features)
########################################################
# 9️⃣ 绘图
library(ggplot2)
ggplot(top_features, aes(x = reorder(feature, value), y = value, fill = view)) +
  geom_bar(stat = "identity") +
  facet_wrap(~factor, scales = "free") +
  coord_flip() +
  theme_bw() +
  labs(title = "Top features contributing to each MOFA factor (Cluster A)",
       x = "Feature",
       y = "Weight")
# 筛选 Factor4
top_features_f4 <- top_features %>%
  filter(factor == "Factor4")
# 绘图
ggplot(top_features_f4, aes(x = reorder(feature, value), y = value, fill = view)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # 横向显示 feature
  theme_bw() +
  labs(
    title = "Factor4 (Cluster B)",
    x = "Feature",
    y = "Weight"
  ) +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(colour = "black",face = "bold"),
    axis.text.y = element_text(size = 10,colour = "black",face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14)
  )
# 筛选 Factor7
top_features_f7 <- top_features %>%
  filter(factor == "Factor7")
# 绘图
ggplot(top_features_f7, aes(x = reorder(feature, value), y = value, fill = view)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # 横向显示 feature
  theme_bw() +
  labs(
    title = "Factor7 (Cluster B)",
    x = "Feature",
    y = "Weight"
  ) +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(colour = "black",face = "bold"),
    axis.text.y = element_text(size = 10,colour = "black",face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14)
  )
####################################################3
########################################################
# 0️⃣ 加载依赖包
library(MOFA2)
library(dplyr)
library(tibble)
library(readr)
library(reticulate)
########################################################
# 1️⃣ 读取数据
cytokine <- read.csv("ccytokine.csv", header = TRUE, row.names = 1, check.names = FALSE)
positive <- read.csv("cpositive.csv", header = TRUE, row.names = 1, check.names = FALSE)
zscore <- read.csv("czscore.csv", header = TRUE, row.names = 1, check.names = FALSE)
# 1️⃣b 读取分簇信息
cluster_info <- read.csv("LN_candidate_cluster_assignments.csv", header = TRUE, stringsAsFactors = FALSE)
# 只保留 A 簇样本
samples_A <- cluster_info %>% filter(Cluster == "B") %>% pull(Sample)
# 1️⃣c 保留共同样本 & 只用 A 簇
common_samples <- Reduce(intersect, list(colnames(cytokine), colnames(positive), colnames(zscore), samples_A))
cytokine <- cytokine[, common_samples]
positive <- positive[, common_samples]
zscore <- zscore[, common_samples]
########################################################
# 2️⃣ 数据类型处理
cytokine_matrix <- as.matrix(cytokine)
storage.mode(cytokine_matrix) <- "double"
zscore_matrix <- as.matrix(zscore)
storage.mode(zscore_matrix) <- "double"
positive_matrix <- as.matrix(positive)
positive_matrix <- ifelse(positive_matrix != 0, 1, 0)
storage.mode(positive_matrix) <- "integer"
########################################################
# 3️⃣ 构建 MOFA 输入列表
data_list <- list(
  Cytokine = cytokine_matrix,
  Peptide_Zscore = zscore_matrix,
  Peptide_Positive = positive_matrix
)
########################################################
# 4️⃣ 创建 MOFA 对象
mofa_obj <- create_mofa(data_list)
# 5️⃣ 设置参数
data_opts <- get_default_data_options(mofa_obj)
model_opts <- get_default_model_options(mofa_obj)
train_opts <- get_default_training_options(mofa_obj)
model_opts$num_factors <- 14
train_opts$convergence_mode <- "medium"
model_opts$likelihoods <- c(
  Cytokine = "gaussian",
  Peptide_Zscore = "gaussian",
  Peptide_Positive = "bernoulli"
)
mofa_obj <- prepare_mofa(
  object = mofa_obj,
  data_options = data_opts,
  model_options = model_opts,
  training_options = train_opts
)
########################################################
# 6️⃣ 训练 MOFA
use_python("E:/Miniforge3/python.exe", required = TRUE)
reticulate::py_install("mofapy2", pip = TRUE)
mofa_trained <- run_mofa(mofa_obj, use_basilisk = FALSE)
########################################################
# 7️⃣ 提取 factor 矩阵（只包含 A 簇样本）
factor_df <- get_factors(mofa_trained, factors = "all", as.data.frame = TRUE)
########################################################
library(gtsummary)
library(dplyr)
data <- read.csv("clinical_LC.csv")
data <- data[,2:16]
tbl_summary(data) #缺点是默认连续型变量不符合正态分布，统一采用非参数检验，分类变量统一采用皮尔逊卡方检验；优点是考虑缺失值并单独列出来
tbl_summary(data, by = LC) # 根据LCC分层
tbl_summary(data, by = LC) %>% add_p()   #添加p值
########################################################
# 提取特征权重
weights_df <- get_weights(mofa_trained, views = "all", factors = "all", as.data.frame = TRUE)
# 每个 factor 下绝对值前 20 个特征，不管属于哪个 view
top_features <- weights_df %>%
  group_by(factor) %>%
  slice_max(order_by = abs(value), n = 30, with_ties = FALSE)  # with_ties=FALSE 保证正好 20 个
head(top_features)
###########################################先改名字
library(dplyr)
library(stringr)
library(readr)
# 读取 name.csv
# 第一列 Index 对应 feature 中的 AAG编号，第二列 Name 是要替换的名称
name_map <- read.csv("name.csv", header = TRUE, stringsAsFactors = FALSE)
# 处理 top_features
top_features <- top_features %>%
  # 去掉 _Peptide_Zscore
  mutate(feature = str_remove(feature, "_Peptide_Zscore")) %>%
  # 提取编号部分和剩余部分
  mutate(
    AAG_id = str_extract(feature, "AAG\\d+"),
    rest = str_extract(feature, "\\|.*")  # |后面的部分，包括|
  ) %>%
  # 用 name.csv 替换 AAG编号
  left_join(name_map, by = c("AAG_id" = "Index")) %>%
  mutate(
    feature = ifelse(!is.na(Name), paste0(Name, rest), feature)
  ) %>%
  select(-AAG_id, -rest, -Name)
# 查看结果
head(top_features)
########################################################
# 9️⃣ 绘图
library(ggplot2)
ggplot(top_features, aes(x = reorder(feature, value), y = value, fill = view)) +
  geom_bar(stat = "identity") +
  facet_wrap(~factor, scales = "free") +
  coord_flip() +
  theme_bw() +
  labs(title = "Top features contributing to each MOFA factor (Cluster B)",
       x = "Feature",
       y = "Weight")
library(ggplot2)
library(dplyr)
# 筛选 Factor2
top_features_f2 <- top_features %>%
  filter(factor == "Factor2")
# 绘图
ggplot(top_features_f2, aes(x = reorder(feature, value), y = value, fill = view)) +
  geom_bar(stat = "identity") +
  coord_flip() +  # 横向显示 feature
  theme_bw() +
  labs(
    title = "Factor2 (Cluster B)",
    x = "Feature",
    y = "Weight"
  ) +
  theme(
    text = element_text(size = 12),
    axis.text.x = element_text(colour = "black",face = "bold"),
    axis.text.y = element_text(size = 10,colour = "black",face = "bold"),
    legend.text = element_text(size = 10),
    plot.title = element_text(size = 14)
  )
