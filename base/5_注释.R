# 注释--------------------------------------------------------------------------
# 根据疾病对marker基因进行修改
classical.markers <- c(
  "CD79A", "MS4A1", "BANK1",                     # B cell
  "PLVAP", "VWF",                                # Endothelial
  "PLP1", "CLU", "S100B",                        # Enteric glial
  "KRT18", "KRT8", "EPCAM",                      # Epithelial
  "COL1A1", "DCN", "COL3A1",                     # Fibroblast
  "MS4A2", "KIT", "TPSAB1",                      # MAST
  "LYZ", "CD68", "CD14",                         # McDC
  "G0S2", "FCGR3B", "CSF3R",                     # Neutrophil
  "LILRA4", "IL3RA", "CLEC4C",                   # pDC
  "MCAM", "PDGFRB", "CSPG4",                     # Pericyte
  "MZB1", "DERL3", "IGHG2",                      # Plasma B
  "CNN1", "TAGLN", "DES",                        # Smooth muscle cell
  "CD3D", "CD3E", "KLRB1"                        # T_NK
)

# 绘制气泡图
for(res in seq(0.1, 1.2, by = 0.1)) {
  p <- DotPlot(seurat, features = classical.markers, group.by = paste0("RNA_snn_res.", res)) + RotatedAxis()
  ggsave(plot = p, filename =paste0(dir, "/celltype_annotation_dotplot_", res, ".pdf"), width = 20, height = 0.5 * length(unique(seurat[[paste0("RNA_snn_res.", res)]])), units = "cm")
}

# 细胞注释阶段用到的全局变量，根据研究的实际情况进行修改
my_res = 0.3

