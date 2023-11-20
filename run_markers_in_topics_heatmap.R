
rm(list = ls())

path_to_markers <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/cell_marker_24_8_23.csv"
path_to_fit <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/objects/fitted_k_10.rds"
color_gradient <- colorRampPalette(c("white","#800000"))


load_libraries <- function() {



  library(ggplot2)
  library(ggrepel)
  library(RColorBrewer)
  library(ggExtra)
  library(tidyverse)
  library(pheatmap)
  library(RColorBrewer)
  library(pheatmap)
  source("utils.R")
}

modify_markers <- function (markers) {
  # markers is a table. for each row we have a cell type and the markers for this cell type.
  # we will create new table that each row will be a marker and the cell type will be the column
  return_val <- data.frame()
  for (i in 1:length(markers)) {
    tmp <- markers[[i]]
    tmp <- data.frame(tmp)
    tmp$celltype <- names(markers)[i]
    return_val <- rbind(return_val, tmp)
  }
  return (return_val)
}

main <- function() {
  #load_libraries()
    fit <- readRDS(path_to_fit)
    markers <- csv_to_markers(path_to_markers)
    marker_anotation <-modify_markers(markers)
    m  <- modify_markers(markers)
mat_f <- data.frame(fit$F)
mat_f$gene_name <- as.factor(rownames(mat_f))
desired_row_names <- marker_anotation$tmp
mat_f <- mat_f[rownames(mat_f) %in% desired_row_names, ]
mat_f$celltype <- marker_anotation$celltype
mat_f <- mat_f %>% group_by(celltype)
gene_names <- mat_f$gene_name
mat_f <- as.matrix(mat_f)
row.names(mat_f) <- gene_names
rownames(marker_anotation) <- marker_anotation$tmp


h <- mat_f[, -c(ncol(mat_f) - 1, ncol(mat_f))]
new_h <- matrix(as.numeric(h), nrow = nrow(h), ncol = ncol(h))
row.names(new_h) <- gene_names
colnames(new_h) <- c("k1", "k2", "k3", "k4", "k5", "k6", "k7", "k8", "k9", "k10")


# generate the anottation which is mat_f with only the last column
final_anotation <- data.frame(mat_f[, ncol(mat_f)])
colnames(final_anotation) <- c("celltype")
p <- pheatmap(new_h,
         cluster_cols = FALSE,
         cluster_rows = FALSE,
         annotation_row = final_anotation,
         color = color_gradient(5000),
         angle_col = 0,
        #annotation_colors =
              #show_rownames = T,
              # cellwidth = 10,
  #cellheight = 10,
  fontsize_row = 15,
  fontsize_col = 15,

)


  }


main()
