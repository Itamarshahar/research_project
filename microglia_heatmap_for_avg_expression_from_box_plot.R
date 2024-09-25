library(pheatmap)
# install.packages('wordspace')
library(wordspace)
matrix <- matrix(c(
  0.06, 0.04, 0.05, 0.06, 0.03,
  0.03, 0.03, 0.04, 0.05, 0.03,
  0.07, 0.16, 0.01, 0.01, 0.01,
  0.02, 0.02, 0.01, 0.02, 0.01,
  0.01, 0.02, 0.02, 0.03, 0.03,
  0.17, 0.12, 0.05, 0.05, 0.01,
  0.05, 0.06, 0.17, 0.18, 0.13,
  0.04, 0.03, 0.02, 0.04, 0.01,
  0.08, 0.12, 0.11, 0.07, 0.41,
  0.03, 0.04, 0.05, 0.03, 0.02,
  0.01, 0.02, 0.09, 0.05, 0.01,
  0.09, 0.06, 0.14, 0.05, 0.06,
  0.24, 0.18, 0.08, 0.05, 0.03,
  0.04, 0.02, 0.04, 0.2, 0.17,
  0.08, 0.09, 0.12, 0.13, 0.02
  ), nrow = 15, byrow = TRUE)

column_sums <- rowSums(matrix)
matrix <- matrix / column_sums
rownames(matrix) <- paste0("Program.", 1:15)
colnames(matrix) <- c("Alzheimer", "Mild cognitive impairment", "Healthy aging", "Superagers", "Young")
pheatmap(matrix,
         cluster_rows = T,
         cluster_cols = F,
         display_numbers = FALSE,
         color = colorRampPalette(c("white", "red"))(50),
         main = "Heatmap of Mean Values",
         fontsize = 30,
         angle_col= 45,
         edge_color = "white",
          )
# ComplexHeatmap::Heatmap(matrix, name = "Mean Values", col = circlize::colorRamp2(c(-1, 0, 1), c("blue", "white", "red")), show_row_names = TRUE, show_column_names = TRUE, row_names_gp = gpar(fontsize = 10), column_names_gp = gpar(fontsize = 10), column_title = "Clusters", row_title = "Cell Types", row_title_side = "left", column_title_side = "top", column_title_gp = gpar(fontsize = 15), row_title_gp = gpar(fontsize = 15), width = unit(10, "cm"), height = unit(10, "cm"), top_annotation = HeatmapAnnotation(foo = anno_text("Mean Values", angle = 90)), bottom_annotation = HeatmapAnnotation(foo = anno_text("Mean Values", angle = 90)))
# pdf("/Users/itamar_shahar/Downloads/Microglia_Mean Values per diagnosis.pdf",
#     width = 12,
#     height = 8)

ComplexHeatmap::Heatmap(matrix,
                        name = "Mean Values per diagnosis",
                        cluster_columns=F,
                        cluster_rows = T,
                        column_title_side = "top",
                        row_title_side = "right",
                        row_names_gp = grid::gpar(fontsize = 20, font=2),
                        column_names_rot=45,
                        column_names_gp = grid::gpar(fontsize = 18, fontface="bold"),
                        column_title_gp = grid::gpar(fontsize = 25, fontface="bold"),
                        show_row_dend = F,
                      # column_title_gp = grid::gpar(fontsize = 35),
                      #       row_title_gp = grid::gpar(fontsize = 35),
                        # right_annotation=is_shared.k1,
                        column_title ="Mean Values per diagnosis",
                        # col = circlize::colorRamp2(c(0, 0.3076923, 0.6153846), c("white","lightblue", "blue"))) #"firebrick4")))
                        col = circlize::colorRamp2(c(0,0.5), c("white","firebrick4"))
)

  # dev.off()

# "k1" <-
# Mean= 0.06
# Mean= 0.04
# Mean= 0.05
# Mean= 0.06
# Mean= 0.03
# k2 <-
#   Mean= 0.03
# Mean= 0.03
# Mean= 0.04
# Mean= 0.05
# Mean= 0.03
#
#
# k3
# Mean= 0.07
# Mean= 0.16
# Mean= 0.01
# Mean= 0.01
# Mean= 0.01
#
#
# k4
# Mean= 0.02
# Mean= 0.02
# Mean= 0.01
# Mean= 0.02
# Mean= 0.01
#
#
# k5
# Mean= 0.01
# Mean= 0.02
# Mean= 0.02
# Mean= 0.03
# Mean= 0.03
#
# k6
# Mean = 0.17
# Mean= 0.12
# Mean= 0.05
# Mean= 0.05
# Mean= 0.01
#
# k7
# Mean= 0.05
# Mean= 0.06
# Mean= 0.17
# Mean= 0.18
# Mean= 0.13
#
# k8
# Mean= 0.04
# Mean= 0.03
# Mean= 0.02
# Mean= 0.04
# Mean= 0.01
#
#
# k9
#
#
# Mean= 0.08
# Mean= 0.12
# Mean= 0.11
# Mean= 0.07
# Mean= 0.41
#
# k10
# Mean= 0.03
# Mean= 0.04
# Mean= 0.05
# Mean= 0.03
# Mean= 0.02
#
# k11
# Mean= 0.01
# Mean= 0.02
# Mean= 0.09
# Mean= 0.05
# Mean= 0.01
#
# k12
# Mean= 0.09
# Mean= 0.06
# Mean= 0.14
# Mean= 0.05
# Mean= 0.06
#
#
# k13
# Mean= 0.24
# Mean= 0.18
# Mean= 0.08
# Mean= 0.05
# Mean= 0.03
#
# k14
# Mean= 0.04
# Mean= 0.02
# Mean= 0.04
# Mean= 0.2
# Mean= 0.17
#
# k15
# Mean= 0.08
# Mean= 0.09
# Mean= 0.12
# Mean= 0.13
# Mean= 0.02