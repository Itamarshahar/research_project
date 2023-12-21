# filter object (with anael)
#load the markers
get_markers.of.ct.brain <- function(){
  # Merge of Gilad and https://www.pnas.org/doi/10.1073/pnas.2008762117
  list(
    Oligodendrocytes=c("MBP", "MOG", "MAG", "ST18", "PLP1", "CTNAA3", "MBP", "PIP4K2A") %>% unique(),
    Microglia=c("C3","PTPRC", "TREM2", "LRMDA", "DOCK8", "ARHGAP24", "PLXDC2") %>% unique(),
    Pericytes=c("PDGFRB","DCN") %>% unique(),
    OPCs=c("PCDH15", "MEGF11", "VCAN", "PDGFRA", "CSPG4") %>% unique(), # CSPG4 = NG2
    Astrocytes=c("SLC1A2","SLC1A3","APOE","GJA1","GFAP", "CD44", "ALDH1L1", "SLC1A2", "ADGRV1", "GPC5", "RYR3", "GFAP") %>% unique(),
    Endothelial=c("FLT1", "CLDN5","ABCB1","ATP10A", "CLDN5", "FLT1", "ABCB1", "EBF1") %>% unique(),
    GABAergic=c("MEG3", "PVALB","SST","VIP", "KIT","GAD2", "NXPH1", "LHFPL3", "GRIK1", "ADARB2") %>% unique(),
    Glutamatergic=c("SLC17A7", "RORB", "TOX","FOXP2", "CUX2", "RALYL", "KCNIP4", "CBLN2", "LDB2", "KCNQ5") %>% unique())
}

get_filtered_obj <- function(obj, path_to_plots, path_to_objs){
  obj <- NormalizeData(obj)
  DefaultAssay(obj) <- "RNA" 
  obj <- FindVariableFeatures(obj)
  obj <- ScaleData(obj)
  obj <- RunPCA(obj, npcs = 30)
  obj <- FindNeighbors(obj, dims = 1:15)
  obj <- RunUMAP(obj, dims = 1:15)
  obj <- FindClusters(obj, resolution = 0.2)
  
  #visualization the data
  ggsave("umap_microglia.pdf",
         plot = DimPlot(obj, label = T),
         width = 10,
         height = 7,
         limitsize = FALSE,
         path = path_to_plots)
  ggsave("dotplot_markers_microglia.pdf",
         plot = DotPlot(obj, cols = c("yellow", "purple"), features = get_markers.of.ct.brain()) + RotatedAxis(),
         width = 13,
         height = 7,
         limitsize = FALSE,
         path = path_to_plots)
  ggsave("VlnQuality_microglia.pdf",
         plot = VlnPlot(obj, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0),
         width = 9,
         height = 5,
         limitsize = FALSE,
         path = path_to_plots)
  #conclussion
  # 3 is doubletswith oligodendrocyte (dotplot + vlnplots)
  # 7 is ambient - dotplot?
  # 9 - doublets of astrocytes - mic
  
  des.8 <- FindMarkers(obj, ident.1 = 8,ident.2 = 0:5)
  # Macrophages
  des.4 <- FindMarkers(obj, ident.1 = 4,ident.2 = c(0, 1, 2, 3))
  des.6 <- FindMarkers(obj, ident.1 = 6,ident.2 = c(0, 1, 2, 3))
  des.6 %>% top_n(20)
  
  # Subset original object
  obj_subset <- subset(x = obj, subset = seurat_clusters %in% c(0, 1, 2, 4, 5, 6, 8))
  
  # run again after filtered
  obj_subset <- NormalizeData(obj_subset)
  DefaultAssay(obj_subset) <- "RNA" 
  obj_subset <- FindVariableFeatures(obj_subset)
  obj_subset <- ScaleData(obj_subset)
  obj_subset <- RunPCA(obj_subset, npcs = 30)
  obj_subset <- FindNeighbors(obj_subset, dims = 1:15)
  obj_subset <- RunUMAP(obj_subset, dims = 1:15)
  obj_subset <- FindClusters(obj_subset, resolution = 0.2)
  
  #visualization the data
  ggsave("umap_microglia_filtered.pdf",
         plot = DimPlot(obj_subset, label = T),
         width = 10,
         height = 7,
         limitsize = FALSE,
         path = path_to_plots)
  ggsave("dotplot_markers_microglia_filtered.pdf",
         plot = DotPlot(obj_subset, cols = c("yellow", "purple"), features = get_markers.of.ct.brain()) + RotatedAxis(),
         width = 13,
         height = 7,
         limitsize = FALSE,
         path = path_to_plots)
  ggsave("VlnQuality_microglia_filtered.pdf",
         plot = VlnPlot(obj_subset, features =  c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0),
         width = 9,
         height = 5,
         limitsize = FALSE,
         path = path_to_plots)
  
  #save the object after filtered as seurat file
  saveRDS(obj_subset,  paste0(path_to_objs,"filtered_microglia.rds"))
  SaveH5Seurat(
    obj_subset,
    filename = paste0(path_to_objs, "filtered_microglia.h5Seurat"),
    overwrite = FALSE,
    verbose = TRUE
  )
  
  return (obj_subset)
}

