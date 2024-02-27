

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
  source("utils.R")
  library(Matrix)
  library(fastTopics)
  library(cowplot)
  library(MASS)
}

run_preprocess <- function(obj, path_to_plots) {
  # add the SampleID_Diagnosis column to the meta data
  obj@meta.data["SampleID_Diagnosis"] <- paste(obj$SampleID, obj$Diagnosis, sep = "_")
  obj [["SampleID_Diagnosis"]] <- paste(obj$SampleID, obj$Diagnosis, sep = "_")
  # adding the MT precantage to the meta data
  obj [["percent.mt"]] <- PercentageFeatureSet(obj, pattern = "^MT-", assay = "RNA")

  return(obj)
}

run_transformation <- function(obj) {
  obj <- SCTransform(obj,
                     #conserve.memory = TRUE,
                     vars.to.regress = "percent.mt",
  )
  return(obj)
}


main <- function(path_to_obj, run_filler= "NA", run_subset ="NA") {
  path_to_obj <- "/Users/shmuel/SuperAgerRemoveSample7264-2Microglia.h5seurat"
  path_to_plots <- "/Volumes/habib-lab/shmuel.cohen/microglia/plots/"
  path_to_objs <- "/Volumes/habib-lab/shmuel.cohen/microglia/objects/"
  
  load_libraries()
  obj <- LoadH5Seurat(path_to_obj)
  obj <-set_sample_and_diagnosis_order(obj)
  
  #filder microglia
  #if run_subset{}
  source("~/Desktop/project/research_project/microglia_filter_redundence.R")
  obj_subset <- get_filtered_obj(obj, path_to_plots, path_to_objs)
  obj_subset <- run_preprocess(obj_subset) # todo maybe run before make subset 
  # obj_subset <- run_transformation(obj_subset)
  
  #### here we need to run the fit_topic_model but it takes time so we skip it ####

  #list of the links to the fit files
  fit_files_paths <- c(
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_8.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_9.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_10.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_11.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_12.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_13.rds",
    #"/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_14.rds",
    "/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds"
    
  )
  
  #run visualization topic model
  source("~/Desktop/project/research_project/microglia_visualization_topic_model.R")
  for (fit_file in fit_files_paths) {
    fit <- readRDS(fit_file)
    k <- as.integer(dim(fit$L)[2])
    path_to_plots_for_k <- glue(path_to_plots, "Plots_for_k={k}/")
    run_main_flow(obj_subset, fit, path_to_plots_for_k) #from ~/Desktop/project/research_project/microglia_visualization_topic_model.R"
  }
  
  #run correlation between the topics
  source("~/Desktop/project/research_project/microglia_topic_correlation.R")
  run_topic_evaluation(fit_files_paths, path_to_plots) #from "~/Desktop/project/research_project/microglia_topic_correlation.R"
  
  #check correlation vs cortex topics
  #todo - remove or enter to condition
  source("~/Desktop/project/research_project/correlation_by_genes.R")
  cortex_fit_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/all_microglia_topic_fit.15.RDS")
  correlation_with_cortex(obj = obj,
                          cortex_fit_15 = cortex_fit_15,
                          fits_list = fit_files_paths,
                          path_to_plots = "/Users/shmuel/microglia/plots/gene_correlation/correlation_with_cortex_1000t/",
                          path_to_predicted = "/Volumes/habib-lab/shmuel.cohen/microglia/objects/predict_1000t.RDS")
  
  #run de, use Adi function
  source("~/Desktop/project/research_project/compare_by_de.R")
  hippocampus_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")
  cortex_15 <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/cortex_all_microglia_topic_fit.15.RDS")
  for (name in c("column", "row", "NA") ){
    generate_intersection_gene_heatmap(obj= obj,
              hippocampus =  hippocampus_15,
              cortex = cortex_15,
              #path_to_plots = "/Users/shmuel/microglia/plots/gene_correlation/correlation_de_adi_15/",
              path_to_plots = "/Volumes/habib-lab/shmuel.cohen/microglia/plots/correlation/de/",
              scale = name)  
    
  }
  
  #run pathways
  #...
  
}
main()

###########3

###################draft
obj <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/filtered_microglia.rds")

