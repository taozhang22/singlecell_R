# 质量控制----------------------------------------------------------------------
seurat[["percent.mt"]] <- PercentageFeatureSet(seurat, pattern = "^MT-")
p1 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "percent.mt")
p2 <- FeatureScatter(seurat, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
p <- p1 + p2; p
ggsave(plot = p, filename = paste0(dir, "/qc_nCount_mt_nFeature.pdf"), width = 60, height = 20, units = "cm")

p <- VlnPlot(object = seurat, features = c("nCount_RNA", "nFeature_RNA", "percent.mt"), layer = "counts", ncol = 1); p
ggsave(plot = p, filename = paste0(dir, "/qc_nCount_nFeature_mt.pdf"), width = 0.6 * length(unique(seurat$Sample)), height = 20)
seurat <- subset(seurat, subset = percent.mt < 20) # 根据情况修改质量控制的条件
qs_save(seurat, paste0(dir, "qc.qs2"))
