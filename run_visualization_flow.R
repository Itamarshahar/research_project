main <- function () {
  rm(list = ls())

  source("run_topic_evaluation.R")
  run_topic_evaluation()

  source("run_visualization_topic_model.R")
  run_visualization_topic_model(15)
  run_visualization_topic_model(10)

  #source("run_de_visualzation.R")
  #run_de_visualzation()
}

main()






