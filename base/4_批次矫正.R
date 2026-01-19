# 去除批次效应------------------------------------------------------------------
seurat <- NormalizeData(seurat)
seurat <- FindVariableFeatures(seurat, nfeatures = 3000)
seurat <- ScaleData(seurat)
seurat <- RunPCA(seurat)
p <- ElbowPlot(object = seurat, ndims = 50); p
ggsave(plot = p, filename =paste0(dir, "/batch_correction_elbowplot.pdf"), width = 20, height = 20, units = "cm")
DimPlot(seurat, group.by = "Sample")
ggsave(plot = p, filename =paste0(dir, "/batch_correction_before_elbowplot.pdf"), width = 20, height = 20, units = "cm")

seurat <- IntegrateLayers(object = seurat, method = CCAIntegration, orig.reduction = "pca", new.reduction = "integrated.cca", verbose = T)

seurat[["RNA"]] <- JoinLayers(seurat[["RNA"]])
seurat <- RunUMAP(seurat, group.by = paste0("integrated_snn_res.", seq(0.1, 1, by = 0.1)), dims = 1:50, reduction = "integrated.cca")
seurat <- FindNeighbors(object = seurat, reduction = "integrated.cca", dims = 1:50)
seurat <- FindClusters(seurat, resolution = seq(0.1, 1, by = 0.1))
qs_save(seurat, paste0(dir, "batch_correction.qs2"))
