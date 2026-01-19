rm(list = ls()); gc()
library(data.table)
library(tidyverse)
library(Seurat)
library(qs2)
library(patchwork)
setwd("D:/research/bioinformatics/singlecell/database/GSE144735")

filename_seurat <- "GSE144735_processed_KUL3_CRC_10X_raw_UMI_count_matrix.txt.gz"
filename_metadata <- "GSE144735_processed_KUL3_CRC_10X_annotation.txt.gz"
filename_out <- "GSE144735"

# 创建seurat对象
seurat <- fread(filename_seurat) %>% 
  column_to_rownames(var = "Index") 

metadata <- fread(filename_metadata) %>% 
  mutate(Sample = word(Index, 1, sep = "_")) %>% 
  left_join(fread(paste0(filename_out, ".txt")), by = "Sample") %>% 
  column_to_rownames(var = "Index") %>% 
  transmute(Sample, MSI, Class)
metadata <- metadata[colnames(seurat), ]
table(rownames(metadata) == names(seurat)) # 查看seurat的列名和metadata的行名是不是完全一致，包括顺序

seurat <- CreateSeuratObject(counts = seurat, project = "filename_out", min.cells = 3, min.features = 100, meta.data = metadata)
qs_save(seurat, paste0(filename_out, ".qs2"))
