rm(list = ls()); gc()
library(data.table)
library(tidyverse)
library(Seurat)
library(qs2)
setwd("D:/research/bioinformatics/singlecell/dataset/GSE144735")

# 读取表达矩阵
seurat <- fread("GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz") %>% 
  column_to_rownames(var = "Index") 

# 读取metadata数据
metadata <- fread("GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz") %>% 
  mutate(Sample = word(Index, 1, sep = "_")) %>% 
  left_join(fread("GSE144735.txt"), by = "Sample") %>% 
  column_to_rownames(var = "Index") %>% 
  transmute(Sample, MSI, Class)
metadata <- metadata[colnames(seurat), ] # 让metadata里面的行名与表达矩阵的列名完全的一致（包括顺序）

seurat <- CreateSeuratObject(counts = seurat, project = "GSE144735", min.cells = 3, min.features = 100, meta.data = metadata) # 创建seurat对象
qs_save(seurat, "GSE144735.qs2")
