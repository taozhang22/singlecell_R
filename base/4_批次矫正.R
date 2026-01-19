# 去除批次效应------------------------------------------------------------------
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
p <- ElbowPlot(object = seurat, ndims = 50); p
ggsave(plot = p, filename =paste0(dir, "/batch_correction_elbowplot.pdf"), width = 20, height = 20, units = "cm") # 保存肘状图
seurat <- RunUMAP(seurat, dims = 1:50, reduction = "pca")
p <- DimPlot(seurat, group.by = "Sample", dims = 1:50) 
ggsave(plot = p, filename =paste0(dir, "/batch_correction_before_umap.pdf"), width = 20, height = 20, units = "cm") # 保存批次矫正前的umap图

seurat <- IntegrateLayers(object = seurat, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = T) # 批次矫正

seurat[["RNA"]] <- JoinLayers(seurat[["RNA"]])
seurat <- RunUMAP(seurat, dims = 1:50, reduction = "integrated.cca")
p <- DimPlot(seurat, group.by = "Sample", dims = 1:50) 
ggsave(plot = p, filename =paste0(dir, "/batch_correction_after_umap.pdf"), width = 20, height = 20, units = "cm") # 保存批次矫正后的umap图

seurat <- FindNeighbors(object = seurat, reduction = "integrated.cca", dims = 1:50)
seurat <- FindClusters(seurat, resolution = seq(0.1, 1, by = 0.1))
p <- DimPlot(seurat, group.by = paste0("integrated.cca_snn_res.", seq(0.1, 1, by = 0.1)), dims = 1:50, ncol = 3)
ggsave(plot = p, filename =paste0(dir, "/batch_correction_after_umap.pdf"), width = 20, height = 20, units = "cm") # 保存批次矫正后的umap图
qs_save(seurat, paste0(dir, "batch_correction.qs2"))
