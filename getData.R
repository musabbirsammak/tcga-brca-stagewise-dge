library(dplyr)
library(TCGAbiolinks)
library(SummarizedExperiment)

if (!dir.exists('Data')) {
  dir.create('Data')
}

# Downloads TCGA-BRCA expression data
query.exp <- GDCquery(project = "TCGA-BRCA",
                      data.category = "Transcriptome Profiling",
                      data.type = "Gene Expression Quantification", 
                      workflow.type = "HTSeq - Counts")
GDCdownload(query.exp)
dataExpression <- GDCprepare(query.exp)

# Separates clinical, expression, and annotation data
clinicalData <- as.data.frame(colData(dataExpression))
expData <- as.data.frame(assay(dataExpression))
annotData <- as.data.frame(rowRanges(dataExpression))

clinicalData <- select(clinicalData, -c(treatments, disease_type, primary_site))
expData <- select(expData, clinicalData$barcode)

# Saves intermediate data
write.csv(annotData, './Data/annotations.csv')
write.csv(clinicalData, "./Data/clinical.csv")
write.csv(expData, "./Data/expression.csv")

# Clears memory
rm(list = ls())
