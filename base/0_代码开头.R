rm(list = ls()); gc()
library(data.table)
library(tidyverse)
library(qs2)
library(Seurat)
library(patchwork)
setwd("D:/research/bioinformatics/singlecell/R/practice/") #根据情况修改工作路径

dir <- "result/base" # 根据情况修改文件夹的名字
dir.create(dir, showWarnings = FALSE, recursive = TRUE)
