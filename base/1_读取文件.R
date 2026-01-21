# 读取文件----------------------------------------------------------------------
ids <- c("GSE132465", "GSE144735", "GSE144735")
seurat <- lapply(ids, function(id) {
  seu <- qs_read(paste0("../../database/", id, "/", id, ".qs2"))
  
  return(seu)
})
seurat <- merge(seurat[[1]], seurat[-1])
