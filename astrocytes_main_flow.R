check_dependencies <- function() {
  # Check if BiocManager package is installed
  if (!requireNamespace("BiocManager", quietly = TRUE)) {
    install.packages("BiocManager")
  }
  if (!requireNamespace("glmGamPoi", quietly = TRUE)) {
    BiocManager::install("glmGamPoi")
  }
}

load_libraries <- function() {
  library(logger)
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

generate_astrocytes_filltered_obj <- function(obj, path_to_save=Null, to_save=TRUE) {
  obj_subset <- subset(x = obj, subset = seurat_clusters %in% c(0, 1, 2, 3, 5, 6, 9, 10, 11, 12))
  obj_subset <- subset(obj_subset, subset = nFeature_RNA > 200 & percent.mt < 10)
  if (to_save & path_to_save) {
    saveRDS(obj_subset, path_to_save)
  }
  return (obj_subset)
}
load_rds_or_seurat <- function(path) {
  if (grepl(".rds", path)) {
    log_info(paste("Loading RDS Object", resolution))
    obj <- readRDS(path)
    return (obj)
  }
  if (grepl(".h5seurat", path)) {
    log_info(paste("Loading Serat Object", resolution))
    obj <- LoadH5Seurat(path)
        return (obj)
  }


main <- function(path_to_obj, run_filler = "NA", run_subset = "NA", generate_plots = TRUE) {
  resolution <- 0.3 # Hyper parametter
  path_to_obj <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/super_agers_astrocytes_sctransform.rds"
  path_to_plots <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/plots/QC_with_sctransform/"
  check_dependencies()
  load_libraries()
  obj <- load_rds_or_seurat(path_to_obj)
  obj <- generate_preprocessed_obj(obj = obj,
                                   resolution = resolution,
                                   run_sctransform_data = TRUE,
                                   run_find_variable_features = TRUE,
                                   run_pca = TRUE,
                                   run_find_neighbors = TRUE,
                                   run_umap = TRUE,
                                   run_find_clusters = FALSE)
    log_info(paste("Running resolution", resolution))
    obj <- generate_preprocessed_obj(obj = obj,
                                     resolution = resolution,
                                     run_find_clusters = TRUE)
    log_info(paste("Generating Plots", resolution))
    plot_preprocess_results(obj = obj, path_to_plots = path_to_plots, resolution = resolution)
    log_info(paste("Done Generating Plots", resolution))
    log_info(paste("Done with running astrocytes_main_flow.R", resolution))
  }
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



# TODO remove this
helper_for_checks <- function() {
  check_dependencies()
  load_libraries()
  path_to_obj <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/super_agers_astrocytes_sctransform.rds"
  path_to_plots <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/plots/QC_with_sctransform/"
  check_dependencies()
  load_libraries()
  obj <- readRDS(path_to_obj)

  obj <- generate_preprocessed_obj(obj = obj,
                                   resolution = resolution,
                                   run_find_variable_features = TRUE,
                                   run_pca = TRUE,
                                   run_find_neighbors = TRUE,
                                   run_umap = TRUE,
                                   run_find_clusters = FALSE)
  for (resolution in c(0.1, 0.3, 0.6, 0.9)) {
    log_info(paste("Loading Object", resolution))
    log_info(paste("Running resolution", resolution))
    obj <- generate_preprocessed_obj(obj = obj,
                                     resolution = resolution,
                                     run_find_clusters = TRUE)
    log_info(paste("Generating Plots", resolution))
    plot_preprocess_results(obj = obj, path_to_plots = path_to_plots, resolution = resolution)
    log_info(paste("Done Generating Plots", resolution))
    log_info(paste("Done with the run", resolution))
  }
}

