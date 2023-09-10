# Research Project
This package is the 3rd year research project as part of the BSc in Computer Science and Computational Biology degree at the Hebrew University in Jerusalem,Israel.

The project use supervised by Dr. Naomi Habib. 

## Pre-requisites
### Packages


## Project Description
### 1. Clustering flow 
### 2. run_fit_topics_model.R:
This script runs the fit_topics_model flow of the SuperAgers database. 
The script is using the [fastTopics](https://github.com/stephenslab/fastTopics).

The script gets 5 arguments:
 1. ```input_path``` - the path to the input file. Should be [Seurat](https://satijalab.org/seurat/) object (```.h5seurat``` postfix).
 2. ```output_path``` - the path to the output file.
 3. ```K``` - the number of topics.
 4. ```NC``` - the number of cores to use for the process. 
 5. ```TRANSPOSE``` - 1 if the _counts matrix_ need to be transposed, 0 otherwise.


### 3. run_de.R:
This script runs the differential expression flow of the SuperAgers database.
  1. ```path_to_obj``` - the path to the input file. Should be [Seurat](https://satijalab.org/seurat/) object (```.h5seurat``` postfix).
  2. ```path_to_fit``` - the path to the fitted model (Can be the ```run_fit_topics_model.R``` output file.)
  3. ```path_to_output``` - the path to the output file.
  4. ```LFC_type``` - the type of the calculation of the LogFoldChange (vsnull, k, le) more info in the [fastTopics](https://stephenslab.github.io/fastTopics/reference/de_analysis.html))
  5. ```NC``` - the number of cores to use for the process.



## Authors
