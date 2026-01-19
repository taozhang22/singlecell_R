# 质量控制图
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
ggsave(plot = p, filename = "nCount_mt.pdf", height = 5 * size, width = 4 * size)

# 绘制nCount_RNA与nFeature_RNA的关系图
p_list <- list()
for (sample in unique(seurat$Sample)) {
  seu <- subset(seurat, subset = Sample == sample)
  p <- FeatureScatter(seu, feature1 = "nCount_RNA", feature2 = "nFeature_RNA") + RotatedAxis() + NoLegend() + labs(subtitle = sample)
  p_list[[sample]] <- p
}
p <- wrap_plots(p_list, nrow = round(root))
p
ggsave(plot = p, filename = "nCount_nFeature.pdf", height = 4 * size, width = 4 * size)
