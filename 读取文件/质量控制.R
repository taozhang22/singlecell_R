# 质量控制图--------------------------------------------------------------------
root <- sqrt(length(unique(seurat$Sample)))

# 绘制nCount_RNA与percent.mt的关系图
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
p_list <- list()
for (sample in unique(seurat$Sample)) {
  seu <- subset(seurat, subset = Sample == sample)
  p <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "percent.mt") + RotatedAxis() + NoLegend() + labs(subtitle = sample)
  p_list[[sample]] <- p
}
p <- wrap_plots(p_list, nrow = round(root))
p
ggsave(plot = p, filename = "nCount_mt.pdf", height = 5 * root, width = 4 * root)

# 绘制nCount_RNA与nFeature_RNA的关系图
p_list <- list()
for (sample in unique(seurat$Sample)) {
  seu <- subset(seurat, subset = Sample == sample)
  p <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + RotatedAxis() + NoLegend() + labs(subtitle = sample)
  p_list[[sample]] <- p
}
p <- wrap_plots(p_list, nrow = round(root))
p
ggsave(plot = p, filename = "nCount_nFeature.pdf", height = 4 * root, width = 4 * root)

# 绘制质量控制的小提琴图("nCount_RNA", "nFeature_RNA", "percent.mt")
p <- VlnPlot(object = seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), layer = "counts", ncol = 1); p
ggsave(plot = p, filename = "nCount_nFeature_mt.pdf", width = 0.6 * length(unique(seurat$Sample)), height = 20)

# 查看每个样本的细胞的数量
df <- table(seurat$Sample)
df <- data.frame(Sample = names(df), CellCount = as.integer(df), row.names = NULL)
ggplot(df, aes(x = reorder(Sample, CellCount), y = CellCount)) +
  geom_col()
