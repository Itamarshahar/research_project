check_dependencies <- function() {
  # Check if BiocManager package is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
}

load_libraries <- function() {
  library(Seurat)
  library(withr)
  library(SeuratDisk)
  #library(SeuratData)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(plotly)
  library(hrbrthemes)
  source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/utils/utils.R")
  library(Matrix)
  library(fastTopics)
  library(cowplot)
  library(MASS)
  library(rhdf5)
}

run_transformation <- function(obj) {
  obj <- SCTransform(obj,
                     #conserve.memory = TRUE,
                     vars.to.regress = "percent.mt",
  )
  return(obj)
}

run_generate_box_plot <- function(obj, fit) {
  source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/boxplot.R")
  generate_box_plot(obj, fit)
}

run_de_flow <- function(hippocampus_exist = FALSE) {

  generate_de <- function(obj, fit, assay = "RNA", organism = "hsa", z_score_log_val = 1, percent_cells_val = 0.01) {
    reweighted_f <- topic_reweight_f(fit$F)
    df <- df_pos_des(obj, reweighted_f, fit, assay = assay, organism = organism)
    des_df <- df %>% dplyr::filter(z_score_log > z_score_log_val & percent_cells > percent_cells_val)
    return(des_df)
  }

  if (!hippocampus_exist) {
    de_hippocampus_X15 <- generate_de(obj = readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/filtered_microglia.rds"),
                                      fit = readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds"))
    saveRDS(de_hippocampus_X15, "/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_hippocampus_X15_.rds")
  }
  else {
    print("de object exists")
  }

}

run_correlation_with_cortex_flow <- function() {
  path <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/microglia_topic_correlation_with_cortex.R"
  source(path)
  run_corralation_with_cortex(fit_hippo = "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds",
                              fit_cortex = "/Volumes/habib-lab/shmuel.cohen/microglia/objects/cortex_all_microglia_topic_fit.15.RDS",
                              de_hippocampus = "/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_hippocampus_X15.rds",
                              de_cortex = "/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_cortex_X15.csv")
}

run_pathways_flow <- function() {
  source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/microglia_pathways.R")
  run_pathways()
}

main <- function(path_to_obj, run_filler = "NA", run_subset = "NA", generate_plots = TRUE) {
  path_to_obj <- "/Users/shmuel/SuperAgerRemoveSample7264-2Microglia.h5seurat"
  path_to_plots <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/plots/QC"
  path_to_objs <- "/Volumes/habib-lab/shmuel.cohen/microglia/objects/"

  #filder microglia
  #if run_subset{}
  # obj <- generate_filtered_obj # TODO combine make generate_filtered_obj() to be the  run_QC_flow() and plot_QC_results()
  # obj <- run_QC_flow(path_to_obj)
  # plot_QC_results(obj, path_to_plots)

  check_dependencies()
  load_libraries()
  # for (resolution in c(0.05, 0.1, 0.15, 0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5, 0.55, 0.6, 0.65, 0.7, 0.75, 0.8, 0.85, 0.9, 0.95)) {
  for (resolution in c(0.1, 0.3, 0.6, 0.9)) {
    obj <- generate_preprocessed_obj(obj = obj,
                                     resolution = resolution,
                                     normalize_data = FALSE,
                                     find_variable_features = FALSE,
                                     scale_data = FALSE,
                                     run_pca = FALSE,
                                     find_neighbors = FALSE,
                                     run_umap = FALSE,
                                     find_clusters = TRUE)
    plot_preprocess_results(obj = obj,
                            resolution = resolution
    )
  }

  # source("~/Desktop/project/research_project/microglia_filter_redundence.R")
  obj_subset <- generate_preprocessed_obj(obj, path_to_plots, path_to_objs)
  obj_subset <- run_preprocess(obj_subset) # todo maybe run before make subset
  # obj_subset <- run_transformation(obj_subset)

  #### here we need to run the fit_topic_model but it takes time so we skip it ####

  #list of the links to the fit files


  #run visualization topic model
  # source("~/Desktop/project/research_project/microglia_visualization_topic_model.R")
  # for (fit_file in fit_files_paths) {
  #   fit <- readRDS(fit_file)
  #   k <- as.integer(dim(fit$L)[2])
  #   path_to_plots_for_k <- glue(path_to_plots, "Plots_for_k={k}/")
  #   run_main_flow(obj_subset, fit, path_to_plots_for_k) #from ~/Desktop/project/research_project/microglia_visualization_topic_model.R"
  # }

  #run correlation between the topics
  # source("~/Desktop/project/research_project/microglia_topic_correlation.R")
  # run_topic_evaluation(fit_files_paths, path_to_plots) #from "~/Desktop/project/research_project/microglia_topic_correlation.R"

  #check correlation vs cortex topics
  #todo - remove or enter to condition

  # run_de_flow()
  # run_correlation_with_cortex_flow()
  # generate boxplot
  # run_generate_box_plot()

  #run pathways
  # run_pathways_flow
  #...

}

# path_to_obj <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/SuperAgerRemoveSample7264-2Astrocytes.h5seurat"
path_to_obj <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/super_agers_astrocytes.h5Seurat"
load_libraries()
obj <- LoadH5Seurat(path_to_obj)
# The [[ operator can add columns to object metadata. This is a great place to stash QC stats
backup_obj <- obj


path_to_rds_super_agers_astrocytes <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/super_agers_astrocytes.rds"
# path_to_rds_super_agers_astrocytes <- "/Users/itamar_shahar/Downloads/super_agers_astrocytes.rds"
# read_rds <-
#"super_agers_astrocytes"
#"superAgers_astrocytes_filtered"
# saveRDS(obj,
#         path_to_rds_super_agers_astrocytes)
# path_to_rds_obj <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/super_agers_astrocytes.rds"
# # path_to_rds_obj <- "/Users/itamar_shahar/Downloads/super_agers_astrocytes.rds"
# obj <- readRDS(path_to_rds_obj)

obj <- NormalizeData(obj)
DefaultAssay(obj) <- "RNA"
obj <- FindVariableFeatures(obj)
obj <- ScaleData(obj)
obj <- RunPCA(obj, npcs = 30)

#obj <- FindNeighbors(obj, dims = 1:15) # TODO look
#obj <- RunUMAP(obj, dims = 1:15)

saveRDS(obj,
        "/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/super_agers_astrocytes.rds")

path_to_h5 <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/super_agers_astrocytes.h5Seurat"
# h5createFile(path_to_h5)
# h5write(obj, path_to_h5, "obj")


SaveH5Seurat(object = obj,
             filename = path_to_h5,
             overwrite = TRUE)

# dims_before <- get_number_of_cells("/Volumes/habib-lab/shmuel.cohen/microglia/objects/SuperAgerRemoveSample7264-2Microglia.h5seurat")
# dims_after <-  get_number_of_cells("/Volumes/habib-lab/shmuel.cohen/microglia/objects/archive/filtered_microglia.h5Seurat")