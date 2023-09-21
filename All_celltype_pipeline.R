---
  title:"Basic Seurat Flow - SuperAgers"
output:html_notebook
---
  ---
    title:"5% DB"
output:html_notebook

---

load_libraries <- function() {
  library(Seurat)
  library(withr)
  library(Seurat)
  library(SeuratDisk)
  library(SeuratData)
  library(ggplot2)
  library(gridExtra)
  library(ggplot2)
  library(dplyr)
  library(plotly)
  library(hrbrthemes)
  source("utils.R")
}

markers_to_csv(cell_markers, "cell_marker_24_8_23.csv")

relevant_cell_clustering_draw_plots <- function() {
  Oligodendrocyte_precursor_cells_plot <- DotPlot(object = obj,
                                                  scale.min = 0,
                                                  features = cell_markers["Oligodendrocyte_precursor_cells"],
                                                  scale.max = 100,
  )
  Astrocytes_plot <- DotPlot(object = obj,
                             features = cell_markers["Astrocytes"],
                             scale.min = 0,
                             scale.max = 100
  )
  Microglia_plot <- DotPlot(object = obj,
                            features = cell_markers["Microglia"],
                            scale.min = 0,
                            scale.max = 100
  )
  Mature_neurons_plot <- DotPlot(object = obj,
                                 features = cell_markers["Mature_neurons"],
                                 scale.min = 0,
                                 scale.max = 100
  )
  Glutamatergic_neurons_plot <- DotPlot(object = obj,
                                        features = cell_markers["Glutamatergic_neurons"],
                                        scale.min = 0,
                                        scale.max = 100
  )
  GABAergic_neurons_plot <- DotPlot(object = obj,
                                    features = cell_markers["GABAergic_neurons"],
                                    scale.min = 0,
                                    scale.max = 100
  )
  vacular_plot <- DotPlot(object = obj,
                          features = cell_markers["vacular"],
                          scale.min = 0,
                          scale.max = 100
  )
}

draw_and_save_dot_plot <- function(obj)
{
  g <- DotPlot(obj,
               feature = cell_markers) &
    scale_color_viridis_c(option = "turbo") &
    RotatedAxis()

  ggsave("dot_plot.pdf",
         g,
         width = 35,
         height = 15, limitsize = FALSE)
}
plot_marker_expression <- function(cell_types_to_plot) {
  for (cell_name in cell_types_to_plot) {
    g <- FeaturePlot(obj,
                     features = cell_markers[[cell_name]],
                     label = TRUE) + NoAxes()

    #interactive = TRUE)                     + NoAxes()

    #+ NoAxes()
    ggsave(paste0(cell_name, ".pdf"),
           plot = g,
           width = 35,
           height = 15,
           limitsize = FALSE,
           path = "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/plots")
  }
}

plot_and_save_vln_plot <- function(features) {
  vln_plot <- VlnPlot(obj, features = features)
  ggsave(paste0(vln_plot, features, ".pdf"),
         plot = vln_plot,
         width = 35,
         height = 15,
         limitsize = FALSE,
         path = "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/plots")
}

run_label_clusters <- function () {
LabelClusters(plot = plot, id = obj@meta.data$Age)
}

assign_new_names_to_clusters <- function(obj, cluster_names, number_of_clusters) {
  # Check if the number of cluster names matches the number of clusters
  if (length(cluster_names) != number_of_clusters) {
    stop("The number of cluster names must match the number of clusters.")
  }

  # Create a mapping between seurat_clusters and cluster_names
  cluster_mapping <- setNames(cluster_names, seq_len(number_of_clusters))

  # Add a new column to the meta.data using the mapping
  obj@meta.data <- obj@meta.data %>%
    mutate(celltype = cluster_mapping[seurat_clusters])

  # Return the modified Seurat object
  return(obj)
}

plot_and_save_dim_plot_with_cells_names <- function() {
  g <- DimPlot(obj,
               reduction = "umap",
               label = TRUE,
               pt.size = 3,
               label.size = 20,
               group.by = "new_col") +
    NoAxes() +
    NoLegend() +
    ggtitle("Clusters") +
    theme(
      plot.title = element_text(size = 60),  # Adjust title size
      #axis.title = element_text(size = 30)   # Adjust axis label size
    )

  ggsave("umap_with_cells_name.pdf",
         plot = g,
         width = 30,
         height = 20,
         limitsize = FALSE,
         path = "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/3rd_year_project2/plots")
}


run_preprocess <- function(obj, path_to_plots) {
  # add the SampleID_Diagnosis column to the meta data

  obj@meta.data["SampleID_Diagnosis"] <- paste(obj$SampleID, obj$Diagnosis, sep = "_")
  obj [["SampleID_Diagnosis"]] <- paste(obj$SampleID, obj$Diagnosis, sep = "_")


  # QC
  # adding the MT precantage to the meta data
  obj [["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-", assay = "RNA")

  # For this stage there are multiple of ways to look for those outlayers:
  # We can use the `group.by` param and set it to: `Diagnosis` or `SampleID`.
  # In this case we choose to use `group.by = Diagnosis`

  #Visualize QC metrics as a violin plot
  VlnPlot(obj,
          features = c("nFeature_RNA", "nCount_RNA", "percent.mt"),
          ncol = 3,
          #group.by = "SampleID")
          group.by = "Diagnosis") #,
  #idents = 1:5)
  ggsave(filename = "Vln_QC_group_by_Diagnosis.pdf",
         path = path_to_plots,
  )
  sort(unique(obj@active.ident))

  # FeatureScatter is typically used to visualize feature-feature relationships, but can be used for anything calculated by the object, i.e. columns in object metadata, PC scores etc.

  plot1 <- FeatureScatter(obj,
                          feature1 = "nCount_RNA",
                          feature2 = "percent.mt",
                          #group.by = "Diagnosis",
                          group.by = "SampleID")

  plot2 <- FeatureScatter(obj,
                          feature1 = "nCount_RNA",
                          feature2 = "nFeature_RNA",
                          #group.by = "Diagnosis",
                          group.by = "SampleID")

  plot1 + plot2
  ggsave(filename = "FeatureScatter_QC_group_by_Diagnosis.pdf",
         path = path_to_plots,
  )
  summary(obj@meta.data$percent.mt)


  # generate SubSet
  # this is very specific to the data set, make sure you change it accordingly
  # here we will remove cells that have more than 10% of their reads mapped to mitochondrial genes
  obj <- subset(obj, subset = nFeature_RNA > 200 &
    nFeature_RNA < 13000 &
    percent.mt < 10)
  return(obj)
}

run_transformation <- function(obj) {
  obj <- SCTransform(obj,
                     #conserve.memory = TRUE,
                     vars.to.regress = "percent.mt",
  )
  return(obj)
}

run_and_evaluate_pca <- function(obj) {
  #Next we perform PCA on the scaled data.
  # By default, only the previously determined variable features are used as
  # input, but can be defined using `features` argument if you wish to choose
  # a different subset.
  obj <- RunPCA(obj, features =
    VariableFeatures(object = obj))
  return(obj)

  # Examine and visualize PCA results a few different ways
  print(obj[['pca']], dims = 1:5, nfeatures = 5)
  VizDimLoadings(obj, dims = 1:2, reduction = 'pca')
  # Determine the 'dimensionality' of the dataset
  ElbowPlot(obj)
  ggsave(filename = "ElbowPlot.pdf",
         path = path_to_plots,
  )
  # alternative heuristic method
  # NOTE: This process can take a long time for big datasets, comment out for expediency. More approximate techniques such as those implemented in ElbowPlot() can be used to reduce computation time
  #obj <- JackStraw(obj, num.replicate = 100)
  #obj <- ScoreJackStraw(obj, dims = 1:20)

}

run_clustering <- function (obj, number_of_clusters){
  obj <- RunUMAP(obj, dims = 1:number_of_clusters)

# note that you can set `label = TRUE` or use the LabelClusters function to help label individual clusters
DimPlot(obj,
        reduction = 'umap',
        label = T
        #group.by ="Diagnosis",
        #group.by ="SampleID",
) &
  NoAxes() &
  NoLegend()
  ## Cluster the cells
  obj <- FindNeighbors(obj, dims = 1:number_of_clusters)
  obj <- FindClusters(obj, resolution = 0.2)

  Command(obj, "FindClusters")

}

run_find_markers <- function(){
  cluster2.markers <- FindMarkers(obj,
                                ident.1 = 0,
                                min.pct = 0.25)
head(cluster2.markers, n = 5)


cluster2.markers["pct1_div_pct2"] <- abs(cluster2.markers$pct.1 - cluster2.markers$pct.2)
head(cluster2.markers, n = 5)

sorted_cluster2.markers <- cluster2.markers[order(cluster2.markers$pct1_div_pct2, decreasing = TRUE),]
head(sorted_cluster2.markers, n = 40)
}

main <- function(path_to_obj) {
  load_libraries()
  obj <- LoadH5Seurat(path_to_obj)
  obj <- run_preprocess(obj)
  obj <- run_transformation(obj)
  obj <- run_and_evaluate_pca(obj)
  number_of_clusters <- 10
  obj <- run_clustering(obj)
  obj <- assign_new_names_to_clusters(obj, cluster_names, number_of_clusters)
  obj < plot_and_save_dim_plot_with_cells_names(obj, path_to_plots)

}

main()



