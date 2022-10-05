list.of.packages <- c("dplyr", "plyr", "stringr", "ggplot2", "ggpubr")
new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
install.packages(new.packages)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

list.of.bioconductor.packages <- c("TCGAbiolinks", "SummarizedExperiment", 
                                   "DESeq2", "DEGreport", "apeglm",
                                   "EnhancedVolcano", "enrichplot")
new.bioconductor.packages <- list.of.bioconductor.packages[!(list.of.bioconductor.packages %in% installed.packages()[,"Package"])]
if(length(new.bioconductor.packages)) BiocManager::install(new.bioconductor.packages)
