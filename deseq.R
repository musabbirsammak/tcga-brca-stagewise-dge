library(ggplot2)
library(ggpubr)
library(DEGreport)
library(DESeq2)
library(apeglm)
library(dplyr)
library(stringr)
library(patchwork)

# Process clinical data
clin <- read.csv('./Data/clinical_selected.csv')
clin$barcode <- str_replace_all(clin$barcode, '-', '.')
rownames(clin) <- clin$barcode
clin <- filter(clin, prior_treatment == 'No')
clin <- filter(clin, prior_malignancy == 'no')
clin <- select(clin, c('stage', 'shortLetterCode', 'paper_BRCA_Subtype_PAM50'))

# Final Stage-Distribution
table(clin$stage)

# Process expression data
expr <- read.csv('./Data/expression.csv')
row.names(expr) <- expr$X
expr <- select(expr, intersect(colnames(expr), rownames(clin)))

# Correct order
expr <- expr[, order(colnames(expr))]
clin <- clin[order(rownames(clin)), ]

# DESeq2
dds <- DESeqDataSetFromMatrix(countData = as.matrix(expr),
                              colData = clin,
                              design = ~ paper_BRCA_Subtype_PAM50 + shortLetterCode + stage + shortLetterCode:stage)
rm(clin)
rm(expr)

dds_lrt <- DESeq(dds, test="LRT", reduced = ~ paper_BRCA_Subtype_PAM50 + shortLetterCode + stage)

saveRDS(dds_lrt, './Results/dds_lrt.Rds')
dds_lrt <- readRDS('./Results/dds_lrt.Rds')

# Significant Genes
res <- results(dds_lrt)
res_tb <- res %>%
  data.frame() %>%
  mutate(gene = rownames(res)) %>%
  as_tibble()
sig_genes <- res_tb %>% 
  filter(padj < 0.05)


sig_dds <- dds[sig_genes$gene, ]
rlog_mat <- vst(sig_dds, nsub = 100, fitType = 'local')

clusters <- degPatterns(as.matrix(assay(rlog_mat)), metadata = colData(dds),
                        time="stage", col="shortLetterCode", minc = 10)

clusters$df$cluster <- factor(clusters$df$cluster,
                              levels = c(1, 2, 6),
                              labels = c('1', '2', '3'))


clusters$df %>%
  group_by(cluster) %>%
  summarise(freq = length(genes))

clusters$plot$data$title <- factor(clusters$plot$data$title,
                                   labels = c('Group 1 - 22 Genes',
                                              'Group 2 - 40 Genes',
                                              'Group 3 - 21 Genes'))

plot1 <- ggplot(clusters[["normalized"]], aes(shortLetterCode, value, color = shortLetterCode, fill = stage)) +
  geom_boxplot() +
  geom_point(position = position_jitterdodge(dodge.width = 0.9), alpha = 0.5) +
  labs(color = 'Sample Type', fill = 'Stages') +
  xlab('Sample Type') + ylab('Z-Score of Gene Abundance') +
  scale_color_manual(values = c("blue", "black"), aesthetics = "color") +
  theme_pubclean()
plot1
ggsave('./Figures/zscore-stagewise.png', height = 6, width = 8, units = 'in')

plot2 <- clusters$plot + theme_pubclean() +
  ylab('Z-Score of Gene Abundance') +
  labs(colour = 'Sample Type', fill = 'Sample Type') +
  scale_color_manual(values = c("blue", "black"), aesthetics = "colour") +
  scale_color_manual(values = c("blue", "black"), aesthetics = "fill")
plot2
ggsave('./Figures/cluster_groups.tiff', height = 4, width = 8, units = 'in')


combined <- plot1 + plot2 + plot_annotation(tag_levels = 'A')
combined
ggsave('./Figures/combined_groups_zscore.png', height = 6, width = 16, units = 'in')

saveRDS(clusters, './Results/clusters.Rds')
