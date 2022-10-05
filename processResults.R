library(tidyverse)
library(plyr)
library(DESeq2)
library(ggplot2)
library(EnhancedVolcano)

annotation <- read.csv('./Data/annotations.csv')
annotation <- annotation[!duplicated(annotation[, 'gene_name']), ]

# Reduced 
result <- readRDS('./Results/dds_lrt.Rds')
result_lrt <- results(result)

result_lrt_tb <- result_lrt %>%
  data.frame() %>%
  mutate(gene = rownames(result)) %>% 
  as_tibble()

sig_genes <- result_lrt_tb %>% 
  filter(padj < 0.05)

sig_genes$gene_name <- mapvalues(
  sig_genes$gene,
  from = annotation$gene_id,
  to = annotation$gene_name
)
result_lrt_tb$gene_name <- mapvalues(
  result_lrt_tb$gene,
  from = annotation$gene_id,
  to = annotation$gene_name
)

p1 <- result_lrt_tb %>%
  ggplot() +
  geom_point(aes(x = log2FoldChange, y = padj))
p1

write.csv(sig_genes, './Results/sig_genes.csv')
write.csv(result_lrt_tb, './Results/result.csv')

clusters <- readRDS('./Results/clusters.Rds')
clusters_df <- clusters$df
clusters_df$gene_name <- mapvalues(
  clusters_df$genes,
  from = annotation$gene_id,
  to = annotation$gene_name
)

write.csv(clusters_df, './Results/gene_clusters.csv')
clusters_df %>%
  group_by(cluster) %>%
  View()

