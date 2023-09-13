rm(list = ls())

load_libraries <- function() {
  library(ComplexHeatmap)
  library(circlize)
  source("utils.R")
}
run_main_flow <- function() {
  print("Running main flow")
  print("Loading libraries")
  load_libraries()
  path_to_fit_k_15 <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/objects/fitted_topic_model_k_15.rds"
  path_to_fit_k_10 <- "/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/objects/fitted_k_10.rds"
  fit_k_15 <- readRDS(path_to_fit_k_15)
  fit_k_10 <- readRDS(path_to_fit_k_10)
  col_fun = colorRamp2(c(-1, 0, 1), c("blue", "white", "red"))
  col_fun(seq(-20, 20))
  correlation <- claculate_correlation(fit_k_15$L, fit_k_10$L)

  pdf("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/k15/plots/corrlation_k=10_with_K=15.pdf") #,width=3.25,height=3.25)#,units="in",res=1200)
  draw(Heatmap(correlation,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          col = col_fun,
          column_title = paste0("The Corralation Between the K = ", 10, " and k = ", 15, " Topics"),

          column_title_gp = gpar(fontsize = 20, fontface = "bold"),
          name = "Correlation",
          rect_gp = gpar(col = "white", lwd = 2),
          column_names_rot = 45,
          cell_fun = HeatmapHelper_add_values_to_display(correlation = correlation,
                                                         OnlyPositive = TRUE)
  )
  )
  dev.off()
  load_libraries()

  correlation <- claculate_correlation(fit_k_10$L, fit_k_10$L)

  pdf("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/k15/plots/corrlation_k=10_with_k=10.pdf") #,width=3.25,height=3.25)#,units="in",res=1200)
  #k_10_vs_k_10 <-
  draw(Heatmap(correlation,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          col = col_fun,
          column_title = paste0("The Corralation Between the K = ", 10, " and K = ", 10, " Topics"),

          column_title_gp = gpar(fontsize = 20, fontface = "bold"),
          name = "Correlation",
          rect_gp = gpar(col = "white", lwd = 2),
          column_names_rot = 45,
          cell_fun = HeatmapHelper_add_values_to_display(correlation = correlation,
                                                         OnlyPositive = TRUE)
  ))
  dev.off()
  load_libraries()
  correlation <- claculate_correlation(fit_k_15$L, fit_k_15$L)

  pdf("/Users/itamar_shahar/Library/CloudStorage/GoogleDrive-itamar.shahar2@mail.huji.ac.il/My Drive/University/General/3rd_year_project/k15/plots/corrlation_k=15_with_k=15.pdf") #,width=3.25,height=3.25)#,units="in",res=1200)
  #k_15_vs_k_15 <-
  draw(Heatmap(correlation,
          cluster_rows = FALSE,
          cluster_columns = FALSE,
          col = col_fun,
          column_title = paste0("The Corralation Between the K = ", 15, " and K = ", 15, " Topics"),
          column_title_gp = gpar(fontsize = 20, fontface = "bold"),
          name = "Correlataion",
          rect_gp = gpar(col = "white", lwd = 2),
          column_names_rot = 45,
          cell_fun = HeatmapHelper_add_values_to_display(correlation = correlation,
                                                         OnlyPositive = TRUE)
  ))
  dev.off()
  load_libraries()
  print("Done with the heatmaps")

}
run_main_flow()




correlation <- claculate_correlation(fit_k_15$L, fit_k_10$L)



