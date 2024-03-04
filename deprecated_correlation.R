load_libraries <- function()
{
  library(RColorBrewer)
  library(dplyr)
  library(glue)
  library(pheatmap)
  library(circlize)
  library(grid)

}


extract_des <- function(de_table, fit_hippo, fit_cortex) {
  intersection_gene <- intersect(rownames(fit_hippo$F), rownames(fit_cortex$F))

  des <- lapply(1:15, function(i) {
    de_topic <- de_table[de_table$topic == paste0("k", i),]
    genes <- de_topic$gene[de_topic$gene %in% intersection_gene]
    return(genes)
  })

  names(des) <- paste0("k", 1:15)

  return(des) }

generate_jaccrd_matrix <- function(de_hippocampus, de_cortex) {
  des_hippo <- extract_des(de_table = de_hippocampus,
                           fit_hippo = fit_hippo,
                           fit_cortex = fit_cortex)
  des_cortex <- extract_des(de_table = de_cortex,
                            fit_hippo = fit_hippo,
                            fit_cortex = fit_cortex)

  # Compute the Jaccard similarity matrix
  similarity_matrix <- matrix(0, nrow = length(des_hippo), ncol = length(des_cortex), dimnames = list(names(des_hippo), names(des_cortex)))
  for (i in seq_along(des_hippo)) {
    for (j in seq_along(des_cortex)) {
      sets = list(
        des_hippo[[i]],
        des_cortex[[j]]
      )
      similarity_matrix[i, j] <- mlr3measures::jaccard(sets)
    }
  }
  return(similarity_matrix)
}

generate_hypergeomtric_matrix <- function(de_hippocampus, de_cortex, fit_hippo, fit_cortex) {

  extract_des <- function(de_table) {
    des <- lapply(1:15, function(i) {
      de_table[de_table$topic == glue("k{i}"),] %>%
        pull(gene)
    }) %>% setNames(., paste("k", 1:15, sep = ""))
    return(des)
  }

  des_hippo <- extract_des(de_hippocampus)
  des_cortex <- extract_des(de_cortex)

  #calculate the bachground genes
  #intersection_gene <- (intersect(rownames(fit_hippo$F), rownames(fit_cortex$F)))
  des_gene <- list()
  for (i in seq_along(des_hippo)) {
    des_gene <- c(des_gene, des_hippo[[i]])
  }
  for (i in seq_along(des_cortex)) {
    des_gene <- c(des_gene, des_cortex[[i]])
  }
  des_gene <- unique(des_gene)
  #calculate the matrix
  sum_gene <- length(unlist(des_gene))
  hypergeomtric_matrix <- matrix(0, nrow = length(des_hippo), ncol = length(des_cortex), dimnames = list(names(des_hippo), names(des_cortex)))

  for (i in seq_along(des_hippo)) {
    for (j in seq_along(des_cortex)) {
      x <- length(intersect(des_hippo[[i]], des_cortex[[j]]))
      m <- length(des_cortex[[j]])
      n <- sum_gene - m
      k <- length(des_hippo[[i]])
      hypergeomtric_matrix[i, j] <- phyper(q = x, m = m, n = n, k = k, lower.tail = FALSE, log.p = FALSE) * 225
    }
  }
  return(hypergeomtric_matrix)
}


save_to_pdf <- function(hypergeomtric_matrix, jaccard_matrix) {
  file_name <- glue("/Volumes/habib-lab/shmuel.cohen/microglia/plots/Plots_for_k=15/p_value.pdf")
  source("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/research_project/utils/utils.R")
  pdf(file_name)
  pheatmap::pheatmap(t(hypergeomtric_matrix),
                     cluster_rows = F,
                     cluster_cols = F,
                     #cell_fun = HeatmapHelper_add_values_to_display(correlation = t(hypergeomtric_matrix),
                     # OnlyPositive = TRUE))
  )
  dev.off()

  file_name <- glue("/Volumes/habib-lab/shmuel.cohen/microglia/plots/Plots_for_k=15/jaccard.pdf")
  pdf(file_name)
  pheatmap(t(jaccard_matrix), cluster_rows = F, cluster_cols = F, cell_fun = HeatmapHelper_add_values_to_display(correlation = t(jaccard_matrix), OnlyPositive = TRUE))
  dev.off()
}

generate_color_palentte <- function() {
    #col_fun <- colorRamp2(c(0, 0.8),  c("lightblue", "blue"))

  col_fun <- colorRampPalette(rev(brewer.pal(n = 7, name = "RdYlBu")))(100)
  return(col_fun)
}

run_corralation_with_cortex <- function(fit_hippo="/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds",
                                        fit_cortex="/Volumes/habib-lab/shmuel.cohen/microglia/objects/cortex_all_microglia_topic_fit.15.RDS",
                                        de_hippocampus="/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_hippocampus_X15.rds",
                                        de_cortex="/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_cortex_X15.csv") {
  jaccard_matrix <- generate_jaccrd_matrix(de_hippocampus, de_cortex)
  hypergeomtric_matrix <- generate_hypergeomtric_matrix(de_hippocampus, de_cortex, fit_hippo, fit_cortex)
  save_to_pdf(hypergeomtric_matrix, jaccard_matrix)


  load_libraries()
  fit_hippo <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/microglia_fitted_topic_model_k_15.rds")

  de_hippocampus <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_hippocampus_X15.rds")
  de_cortex <- read.csv("/Volumes/habib-lab/shmuel.cohen/microglia/objects/DE_cortex_X15.csv")
  fit_cortex <- readRDS("/Volumes/habib-lab/shmuel.cohen/microglia/objects/cortex_all_microglia_topic_fit.15.RDS")
  jaccard_matrix <- generate_jaccrd_matrix(de_hippocampus, de_cortex)
  hypergeomtric_matrix <- generate_hypergeomtric_matrix(de_hippocampus, de_cortex, fit_hippo, fit_cortex)
  col_fun<- generate_color_palentte()

  # lgd_list <- list(
  #   ComplexHeatmap::Legend(
  #     labels = c("<0.01"), title = "pvalue",
  #     graphics = list(
  #       function(x, y, w, h) grid.text("*", x = x, y = y,
  #                                      gp = gpar(fill = "black")))
  #   ))


  cell_fun <- function(j, i, x, y, w, h, fill) {
    fontsize <- 6
    if (hypergeomtric_matrix[i, j] < 0.01) {
      grid.text(paste0(sprintf("%.2f", jaccard_matrix[i, j]), "*"), x, y, gp = gpar(fontsize = fontsize))
    }

  }

  row_title_side <- c("right")

  pdf("/Volumes/habib-lab/shmuel.cohen/microglia/plots/Plots_for_k=15/Correlation Between Hippocampus and Cortex.pdf")
  ComplexHeatmap::draw(ComplexHeatmap::Heatmap(matrix = as.matrix(jaccard_matrix),
                                name = "jaccard_score",
                                column_title = "Cortex",
                                cluster_rows = FALSE,
                                cluster_columns = FALSE,
                                col = col_fun,
                                cell_fun = cell_fun,
                                #annotation_legend_list = lgd_list,
                                rect_gp = gpar(col = "white", lwd = 1),
                                row_title = "Hippocampus",
                                column_names_rot = 45,
                                row_title_side = c("right"),
                                row_title_rot = switch(row_title_side[1], "left" = 90, "right" = 270),
                                right_annotation = T
                                #row_title_rot = 0
                                #column_names_gp = gpar(fontsize = 10, fontface = "bold", col = "black", fontangle = 315)
                                # cell_fun = HeatmapHelper_add_values_to_display(correlation = matrix, OnlyPositive = TRUE),
  )
  )
  # modify the legend loction
  lgd <- ComplexHeatmap::Legend(
    labels = c("<0.01"),
    title = "pvalue",
    graphics = list(
      function(x, y, w, h) grid.text(label = "*",
                                     x = x,
                                     y = y,
                                     gp = gpar(fill = "black"))))
  ComplexHeatmap::draw(lgd,
                       x = unit(0.93, "npc"),
                       y = unit(0.66, "npc"),
                       just = c("right", "top"))
  dev.off()
}




# run example:
# run_corralation_with_cortex()