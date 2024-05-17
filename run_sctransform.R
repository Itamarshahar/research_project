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
  library(Seurat)
  library(withr)
  library(SeuratDisk)
  #library(SeuratData)
  library(ggplot2)
  library(gridExtra)
  library(dplyr)
  library(plotly)
  library(hrbrthemes)
  # source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/utils/utils.R")
  library(Matrix)
  library(fastTopics)
  library(cowplot)
  library(MASS)
  library(rhdf5)
}

check_dependencies()
load_libraries()
path_to_obj <- "/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/SuperAgerRemoveSample7264-2Astrocytes.h5seurat"
obj <- readH5Seurat(path_to_obj)
obj <- SCTransform(obj = obj,
                       variable.features.n = 5000)
SaveH5Seurat(obj, "/Volumes/habib-lab/shmuel.cohen/astrocytes/objects/super_agers_sctransform.h5seurat")
