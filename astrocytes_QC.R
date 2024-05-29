# filter object (with anael)
#load the markers
get_markers.of.ct.brain <- function() {
  # Merge of Gilad and https://www.pnas.org/doi/10.1073/pnas.2008762117
  list(
    Oligodendrocytes = c("MBP", "MOG", "MAG", "ST18", "PLP1", "CTNAA3", "MBP", "PIP4K2A") %>% unique(),
    Microglia = c("C3", "PTPRC", "TREM2", "LRMDA", "DOCK8", "ARHGAP24", "PLXDC2") %>% unique(),
    Pericytes = c("PDGFRB", "DCN") %>% unique(),
    OPCs = c("PCDH15", "MEGF11", "VCAN", "PDGFRA", "CSPG4") %>% unique(), # CSPG4 = NG2
    Astrocytes = c("SLC1A2", "SLC1A3", "APOE", "GJA1", "GFAP", "CD44", "ALDH1L1", "SLC1A2", "ADGRV1", "GPC5", "RYR3", "GFAP") %>% unique(),
    Endothelial = c("FLT1", "CLDN5", "ABCB1", "ATP10A", "CLDN5", "FLT1", "ABCB1", "EBF1") %>% unique(),
    GABAergic = c("MEG3", "PVALB", "SST", "VIP", "KIT", "GAD2", "NXPH1", "LHFPL3", "GRIK1", "ADARB2") %>% unique(),
    Glutamatergic = c("SLC17A7", "RORB", "TOX", "FOXP2", "CUX2", "RALYL", "KCNIP4", "CBLN2", "LDB2", "KCNQ5") %>% unique(),
    Monocytes = c('CCR2', 'MSR1', "S100A4", "CD63") %>% unique(), #"VCAN"
    Macrophages = c('CD86', 'CSF1R', "MRC1", "CD74", "CDK1") %>% unique(),#"CD63"
    # Neurogenesis=c("DCX", "SOX2", "SOX4", "TOP2A", "EOMES", "NES", "VIM", "PAX6", "HES1", "ASCL1") %>% unique() # , "GFAP"
    RadialGlia=c("NES", "VIM", "PAX6", "HES1") %>% unique()
  )
}

run_QC_flow <- function(obj,
                        add_mt_precentage = FALSE) {
  if (add_mt_precentage) {
    obj [["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-", assay = "RNA")
  }
  return(obj)
}

generate_preprocessed_obj <- function(obj,
                                      resolution = 0.3,
                                      run_QC = FALSE,
                                      run_sctransform_data = FALSE,
                                      run_normalize_data = FALSE,
                                      run_scale_data = FALSE,
                                      run_find_variable_features = FALSE,
                                      run_pca = FALSE,
                                      run_find_neighbors = FALSE,
                                      run_umap = FALSE,
                                      run_find_clusters = FALSE) {
  if (run_QC) {
    log_info("Running QC")
    run_QC_flow(obj = obj, add_mt_precentage = TRUE) }
  if (run_sctransform_data) {
    log_info("Running SCTransform")
    obj <- SCTransform(obj = obj,
                       variable.features.n = 5000
    )
    DefaultAssay(obj) <- "SCT"
  }
  if (run_normalize_data && !run_sctransform_data) {
    log_info("Running NormalizeData")
    obj <- NormalizeData(obj)
  }
  if (run_scale_data && !run_sctransform_data) {
    log_info("Running ScaleData")
    obj <- ScaleData(obj)
  }
  if (run_find_variable_features) {
    log_info("Running FindVariableFeatures")
    obj <- FindVariableFeatures(obj)
  }
  if (run_pca) {
    log_info("Running RunPCA")
    obj <- RunPCA(obj, npcs = 30)
  }
  if (run_find_neighbors) {
    log_info("Running FindNeighbors")
    obj <- FindNeighbors(obj, dims = 1:15)
  }
  if (run_umap) {
    log_info("Running RunUMAP")
    obj <- RunUMAP(obj, dims = 1:15)
  }
  if (run_find_clusters) {
    log_info("Running FindClusters")
    obj <- FindClusters(object=obj, resolution = resolution)
  }
  # DefaultAssay(obj) <- "RNA"
  return(obj)
}

plot_QC_results <- function(obj, path_to_plots) {
  path <- "VlnQuality_astrocytes.pdf"
  ggsave(path,
         plot = VlnPlot(object = obj,
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        pt.size = 0),
         width = 9,
         height = 5,
         limitsize = FALSE,
         path = path_to_plots)
  plot1 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "percent.mt")
  plot2 <- FeatureScatter(obj, feature1 = "nCount_RNA", feature2 = "nFeature_RNA")
  ggsave("percent_mt_and_rna_counts.pdf",
         plot = plot1 + plot2,
         width = 9,
         height = 5,
         limitsize = FALSE,
         path = path_to_plots)
}

plot_umap_preproccess <- function(obj,
                                  path_to_plots="/Volumes/habib-lab/shmuel.cohen/astrocytes/plots/QC_with_sctransform/",
                                  resolution = 0.3) {
  plot1 <- DimPlot(obj, group.by = SEURAT_CLUSTERS, label = TRUE) + theme(
    axis.line = element_blank(),   # Remove axis lines
    axis.text = element_blank(),   # Remove axis text
    axis.title = element_blank()   # Remove axis titles
  )
  # Create the second plot
  plot2 <- DimPlot(obj, group.by = DIAGNOSIS_SAMPLE, label = TRUE, label.size = 2) + theme(
    axis.line = element_blank(),   # Remove axis lines
    axis.text = element_blank(),   # Remove axis text
    axis.title = element_blank()   # Remove axis titles
  )
  # Create the third plot
  obj@meta.data$Gender[obj@meta.data$Gender == ""] = "Unknown"
  plot3 <- DimPlot(obj, group.by = "Gender", label = TRUE, label.size = 3) + theme(
    axis.line = element_blank(),   # Remove axis lines
    axis.text = element_blank(),   # Remove axis text
    axis.title = element_blank()   # Remove axis titles
  )
  # Remove the group titles
  combined_plot <- grid.arrange(plot1, plot2, plot3, ncol = 3,
                                top = glue::glue("Two Dimensional Umap With Resolution ", resolution, "  Colored by Cluster,Individual and Gender")) #+ ggtitle("Two Dimensional Umap With Resolution 0.2 Color by Cluster and by Individual")

  ggsave(glue::glue("umap_astrocytes_resolution_", resolution, ".pdf"),
         plot = combined_plot,
         width = 17, # Adjust as needed
         height = 5, # Adjust as needed
         limitsize = FALSE,
         path = path_to_plots)
}

plot_preprocess_results <- function(obj,
                                    path_to_plots = "/Volumes/habib-lab/shmuel.cohen/astrocytes/plots",
                                    resolution = 0.3) {
  plot_umap_preproccess(obj = obj,
                        resolution = resolution)
  ggsave(glue::glue("dotplot_markers_astrocytes_resolution", "_", resolution, ".pdf"),
         plot = DotPlot(obj, cols = c("yellow", "purple"), features = get_markers.of.ct.brain()) +
           RotatedAxis() +
           theme(text = element_text(size = 9), axis.text.x = element_text(size = 9)),
         width = 15,
         height = 7,
         limitsize = FALSE,
         path = path_to_plots)
  ggsave(glue::glue("VlnQuality_astrocytes", "_", resolution, ".pdf"),
         plot = VlnPlot(obj, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), pt.size = 0),
         width = 9,
         height = 5,
         limitsize = FALSE,
         path = path_to_plots)
  ggsave(glue::glue("VlnQuality_astrocytes_group_by_diagnosis", "_", resolution, ".pdf"),
         plot = VlnPlot(obj = obj,
                        features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
                        pt.size = 0,
                        group.by = "Diagnosis"),
         width = 9,
         height = 5,
         limitsize = FALSE,
         path = path_to_plots)
}

###added after the objects are saved

generate_all_samples_double_score_matrix <- function(path_to_double_score_matrix, relevant_cells) {
  merged_table <- data.frame(
    row.names = relevant_cells,
    doublet.score = rep(0, length(relevant_cells))
  )

  for (path in path_to_double_score_matrix) {
    temp <- h5read(path, "meta.data")

    # Update doublet scores in merged_table
    relevant_rows <- temp[["_index"]] %in% relevant_cells
    relevant_temp <- temp[relevant_rows,]

    merged_table$doublet.score[rownames(relevant_temp)] <- relevant_temp$doublet.score
  }

  return(merged_table)
}

# TODO maybe remove

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


    matching_rownames <- intersect(rownames(rv_obj@meta.data), rownames(merged_table)) #$cell_id)

    #matching_rows <- merged_table[merged_table$cell_id %in% matching_rownames, ]
    #rownames(matching_rows) <- matching_rows[["_index"]]
    #matching_rownames <- rownames(obj@meta.data)[rownames(obj@meta.data) %in% merged_table[["_index"]]]
    rv_obj@meta.data[matching_rownames, "doublet.score"] <- merged_table[matching_rownames,]
  }

  return(obj)
}

#  TODO maybe remove
