# 去除双细胞--------------------------------------------------------------------
options(repr.plot.width = 5, repr.plot.height = 5, repr.plot.res = 300)
seurat <- SplitObject(seurat, split.by = "Sample")
seurat <- lapply((seurat), function(seu) {
  seu <- NormalizeData(seu)
  seu <- FindVariableFeatures(seu, nfeatures = 3000)
  seu <- ScaleData(seu)
  seu <- RunPCA(seu)
  seu <- RunUMAP(seu, dims = 1:50)
  seu <- FindNeighbors(seu, dims = 1:50)
  seu <- FindClusters(seu, resolution = 0.5)
  
  sweep <- paramSweep(seu, PCs = 1:50, sct = FALSE)
  sweep <- summarizeSweep(sweep, GT = FALSE)
  pk <- find.pK(sweep)
  pk <- as.numeric(as.character(pk$pK[which.max(pk$BCmetric)]))
  
  annotations <- seu@meta.data$seurat_clusters
  homotypic.prop <- modelHomotypic(annotations) # 估算同型双细胞的比例
  nExp_poi <- round(ncol(seu)/1000 * 0.008 * nrow(seu@meta.data)) # 算出双细胞的数量
  nExp_poi.adj <- round(nExp_poi * (1 - homotypic.prop)) # 去除同型双细胞的数量
  seu <- doubletFinder(seu, PCs = 1:50, pN = 0.25, pK = pk, nExp = nExp_poi.adj, reuse.pANN = NULL, sct = FALSE)
  
  seu$doublet.class <- seu[[paste0("DF.classifications_0.25_", pk, "_", nExp_poi.adj)]]
  seu[[paste0("DF.classifications_0.25_", pk, "_", nExp_poi.adj)]] <- NULL
  pann <- grep("pANN_", names(seu@meta.data), value = TRUE)
  seu$pANN <- seu[[pann]]; seu[[pann]] <- NULL
  
  return(seu)
})
seurat <- merge(seurat[[1]], seurat[-1])
seurat <- subset(seurat, doublet.class == "Singlet")
seurat <- DietSeurat(seurat, assays = "RNA", layers = "counts") 
qs_save(seurat, paste0(dir, "remove_doublets.qs2"))
