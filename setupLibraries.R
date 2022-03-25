list.of.packages <- c("dplyr", "plyr", "stringr", "TCGAbiolinks", "SummarizedExperiment", 
                      "ggplot2", "ggpubr", "DESeq2", "DEGreport", "apeglm",
                      "EnhancedVolcano", "enrichplot")

new.packages <- list.of.packages[!(list.of.packages %in% installed.packages()[,"Package"])]
if(length(new.packages)) install.packages(new.packages)
