rm(list = ls()); gc()
library(data.table)
library(tidyverse)
library(Seurat)
library(qs2)
library(patchwork)
setwd("D:/research/bioinformatics/singlecell/database/GSE144735") # 根据情况修改工作路径

# 参数，根据情况修改
filename_seurat <- "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz"
filename_metadata <- "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"
filename_out <- "GSE144735"

# 读取表达矩阵
seurat <- fread(filename_seurat) %>% 
  column_to_rownames(var = "Index") 

# 读取metadata数据
metadata <- fread(filename_metadata) %>% 
  mutate(Sample = word(Index, 1, sep = "_")) %>% 
  left_join(fread(paste0(filename_out, ".txt")), by = "Sample") %>% 
  column_to_rownames(var = "Index") %>% 
  transmute(Sample, MSI, Class)
metadata <- metadata[colnames(seurat), ]
table(rownames(metadata) == names(seurat))

# 创建seurat对象
seurat <- CreateSeuratObject(counts = seurat, project = filename_out, min.cells = 0, min.features = 0, meta.data = metadata)
seurat[["RNA"]] <- split(seurat[["RNA"]], f = seurat$Sample)
qs_save(seurat, paste0(filename_out, ".qs2"))
