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
    Glutamatergic=c("SLC17A7", "RORB", "TOX","FOXP2", "CUX2", "RALYL", "KCNIP4", "CBLN2", "LDB2", "KCNQ5") %>% unique(),
    Monocytes = c('CCR2', 'MSR1',  "S100A4", "CD63") %>% unique(), #"VCAN"
    Macrophages= c('CD86','CSF1R',"MRC1", "CD74", "CDK1") %>% unique()) #"CD63"
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
  # Create the first plot
  plot1 <- DimPlot(obj, label = TRUE)+theme(
    axis.line = element_blank(),   # Remove axis lines
    axis.text = element_blank(),   # Remove axis text
    axis.title = element_blank()   # Remove axis titles
  )
  # Create the second plot
  plot2 <- DimPlot(obj, group.by = "SampleID", label = TRUE, label.size = 3)+theme(
    axis.line = element_blank(),   # Remove axis lines
    axis.text = element_blank(),   # Remove axis text
    axis.title = element_blank()   # Remove axis titles
  )
  # Create the third plot
  obj@meta.data$Gender[obj@meta.data$Gender == ""] = "Unknown"
  plot3 <- DimPlot(obj, group.by = "Gender", label = TRUE, label.size = 3)+theme(
    axis.line = element_blank(),   # Remove axis lines
    axis.text = element_blank(),   # Remove axis text
    axis.title = element_blank()   # Remove axis titles
  )
  # Remove the group titles
  combined_plot <- grid.arrange(plot1, plot2, plot3, ncol = 3, top ="Two Dimensional Umap With Resolution 0.2 Colored by Cluster,Individual and Gender" )#+ ggtitle("Two Dimensional Umap With Resolution 0.2 Color by Cluster and by Individual")
  
  # Save the combined plot
  ggsave("umap_microglia.pdf",
         plot = combined_plot,
         width = 17, # Adjust as needed
         height = 5, # Adjust as needed
         limitsize = FALSE,
         path = path_to_plots)
  ggsave("dotplot_markers_microglia.pdf",
         plot = DotPlot(obj, cols = c("yellow", "purple"), features = get_markers.of.ct.brain()) + RotatedAxis()+theme(text =  element_text(size = 9), axis.text.x = element_text(size = 9)),
         width = 15,
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

  # Create the first plot
  plot1 <- DimPlot(obj_subset, label = TRUE)+theme(
    axis.line = element_blank(),   # Remove axis lines
    axis.text = element_blank(),   # Remove axis text
    axis.title = element_blank()   # Remove axis titles
  )
  # Create the second plot
  plot2 <- DimPlot(obj_subset, group.by = "SampleID_Diagnosis", label = TRUE, label.size = 2)+theme(
    axis.line = element_blank(),   # Remove axis lines
    axis.text = element_blank(),   # Remove axis text
    axis.title = element_blank()   # Remove axis titles
  )
  # Create the third plot
  obj_subset@meta.data$Gender[obj_subset@meta.data$Gender == ""] = "Unknown"
  plot3 <- DimPlot(obj_subset, group.by = "Gender", label = TRUE, label.size = 3)+theme(
    axis.line = element_blank(),   # Remove axis lines
    axis.text = element_blank(),   # Remove axis text
    axis.title = element_blank()   # Remove axis titles
  )
  # Remove the group titles
  combined_plot <- grid.arrange(plot1, plot2, plot3, ncol = 3, 
                                top ="Two Dimensional Umap With Resolution 0.2 Colored by Cluster,Individual and Gender" )#+ ggtitle("Two Dimensional Umap With Resolution 0.2 Color by Cluster and by Individual")

  # Save the combined plot
  ggsave("umap_microglia_filtered_comb.pdf",
         plot = combined_plot,
         width = 17, # Adjust as needed
         height = 5, # Adjust as needed
         limitsize = FALSE,
         path = path_to_plots)
  ggsave("dotplot_markers_microglia_filtered.pdf",
         plot = DotPlot(obj_subset, cols = c("yellow", "purple"), features = get_markers.of.ct.brain()) + RotatedAxis()+theme(text =  element_text(size = 9), axis.text.x = element_text(size = 9)),
         width = 15,
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


###added after the objects are saved

generate_all_samples_double_score_matrix <- function(path_to_double_score_matrix, relevant_cells) {
  merged_table <- data.frame(
    row.names = relevant_cells,
    doublet.score = rep(0,length(relevant_cells))
  )
  
  for (path in path_to_double_score_matrix) {
    temp <- h5read(path, "meta.data")
    
    # Update doublet scores in merged_table
    relevant_rows <- temp[["_index"]] %in% relevant_cells
    relevant_temp <- temp[relevant_rows, ]
    
    merged_table$doublet.score[rownames(relevant_temp)] <- relevant_temp$doublet.score
  }
  
  return(merged_table)
}

generate_all_samples_double_score_matrix1 <- function(obj, path_to_double_score_matrix) {
  rv_obj <- obj
  rv_obj@meta.data["doublet.score"] <- 0
  for (path in path_to_double_score_matrix) {
    
    temp <- h5read(path, "meta.data")
    #doublet_scores <- temp["_index","doublet.score"]
    doublet_scores_col1 <- temp[["_index"]]
    doublet_scores_col2 <- temp$doublet.score
    merged_table <- data.frame(
                               doublet.score = doublet_scores_col2,
                               row.names = doublet_scores_col1)
    
   
    matching_rownames <- intersect(rownames(rv_obj@meta.data), rownames(merged_table))#$cell_id)
    
    #matching_rows <- merged_table[merged_table$cell_id %in% matching_rownames, ]
    #rownames(matching_rows) <- matching_rows[["_index"]]
    #matching_rownames <- rownames(obj@meta.data)[rownames(obj@meta.data) %in% merged_table[["_index"]]]
    rv_obj@meta.data[matching_rownames, "doublet.score"] <- merged_table[matching_rownames, ]
  }
  
  return(obj)
}



path_to_double_score_matrix <- c(
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput6962.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput6991.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput6992.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput6998.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput7001.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput7162.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput7182.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput7258.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput7264-1.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput7320.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput7426.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput7436.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutput8126.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutputA-19-08.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutputA-19-79.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutputA20-22.h5seurat",
  "/Volumes/habib-lab/Shared/SuperAgers/objects/DoubletFinder/DoubletFinderOutputA21-135.h5seurat"
)

new_obj <- generate_all_samples_double_score_matrix(path_to_double_score_matrix, rownames(obj@meta.data))
plot <-FeaturePlot(object = obj, features = doublet.score, label=F)
