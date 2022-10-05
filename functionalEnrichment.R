library(tidyverse)
library(ggpubr)
library(ggplot2)
library(stringr)
library(patchwork)

c1_reactome <- read.csv('./Results/cluster_1_reactome.txt', sep = '\t')
c1_reactome$Adjusted.P.value <- -log10(c1_reactome$Adjusted.P.value)
c1_reactome <- arrange(c1_reactome, desc(Adjusted.P.value), desc(Combined.Score))
c1_reactome <- c1_reactome[1:10, ]
c1_reactome$Term <- str_wrap(c1_reactome$Term, 45)
c1_reactome$Term <- factor(c1_reactome$Term,
                           levels = c1_reactome$Term[order(c1_reactome$Adjusted.P.value, c1_reactome$Combined.Score, decreasing = FALSE)])
c1_reactome_plot <- c1_reactome %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c1_reactome_plot
ggsave('./Figures/cluster_1_reactome.tiff', height = 1200, width = 2200, units = 'px', dpi = 300)



c2_reactome <- read.csv('./Results/cluster_2_reactome.txt', sep = '\t')
c2_reactome$Adjusted.P.value <- -log10(c2_reactome$Adjusted.P.value)
c2_reactome <- arrange(c2_reactome, desc(Adjusted.P.value), desc(Combined.Score))
c2_reactome <- c2_reactome[1:10, ]
c2_reactome$Term <- str_wrap(c2_reactome$Term, 45)
c2_reactome$Term <- factor(c2_reactome$Term,
                           levels = c2_reactome$Term[order(c2_reactome$Adjusted.P.value, c2_reactome$Combined.Score, decreasing = FALSE)])

c2_reactome_plot <- c2_reactome %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(Adjusted P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c2_reactome_plot
ggsave('./Figures/cluster_2_reactome.tiff', height = 1400, width = 2200, units = 'px', dpi = 300)



c3_reactome <- read.csv('./Results/cluster_3_reactome.txt', sep = '\t')
c3_reactome$Adjusted.P.value <- -log10(c3_reactome$Adjusted.P.value)
c3_reactome <- arrange(c3_reactome, desc(Adjusted.P.value), desc(Combined.Score))
c3_reactome <- c3_reactome[1:10, ]
c3_reactome$Term <- str_wrap(c3_reactome$Term, 45)
c3_reactome$Term <- factor(c3_reactome$Term,
                           levels = c3_reactome$Term[order(c3_reactome$Adjusted.P.value, c3_reactome$Combined.Score, decreasing = FALSE)])

c3_reactome_plot <- c3_reactome %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(Adjusted P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c3_reactome_plot
ggsave('./Figures/cluster_3_reactome.tiff', height = 1400, width = 2200, units = 'px', dpi = 300)


combined_reactome <- c1_reactome_plot + c2_reactome_plot + c3_reactome_plot
combined_reactome + plot_annotation(title = 'Reactome Pathway Analysis of Identified Gene Clusters')
combined_reactome + plot_annotation(tag_levels = 'A')
ggsave('./Figures/cluster_reactome.tiff', height = 1400, width = 5000, units = 'px', dpi = 300)










c1_bp <- read.csv('./Results/cluster_1_go_bp.txt', sep = '\t')
c1_bp$Adjusted.P.value <- -log10(c1_bp$Adjusted.P.value)
c1_bp <- arrange(c1_bp, desc(Adjusted.P.value), desc(Combined.Score))
c1_bp <- c1_bp[1:10, ]
c1_bp$Term <- str_wrap(c1_bp$Term, 45)
c1_bp$Term <- str_to_title(c1_bp$Term, 45)
c1_bp$Term <- factor(c1_bp$Term,
                     levels = c1_bp$Term[order(c1_bp$Adjusted.P.value, c1_bp$Combined.Score, decreasing = FALSE)])
c1_bp_plot <- c1_bp %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c1_bp_plot
ggsave('./Figures/cluster_1_bp.tiff', height = 1200, width = 2200, units = 'px', dpi = 300)



c2_bp <- read.csv('./Results/cluster_2_go_bp.txt', sep = '\t')
c2_bp$Adjusted.P.value <- -log10(c2_bp$Adjusted.P.value)
c2_bp <- arrange(c2_bp, desc(Adjusted.P.value), desc(Combined.Score))
c2_bp <- c2_bp[1:10, ]
c2_bp$Term <- str_wrap(c2_bp$Term, 45)
c2_bp$Term <- str_to_title(c2_bp$Term, 45)
c2_bp$Term <- factor(c2_bp$Term,
                     levels = c2_bp$Term[order(c2_bp$Adjusted.P.value, c2_bp$Combined.Score, decreasing = FALSE)])
c2_bp_plot <- c2_bp %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c2_bp_plot
ggsave('./Figures/cluster_2_bp.tiff', height = 1200, width = 2200, units = 'px', dpi = 300)



c3_bp <- read.csv('./Results/cluster_3_go_bp.txt', sep = '\t')
c3_bp$Adjusted.P.value <- -log10(c3_bp$Adjusted.P.value)
c3_bp <- arrange(c3_bp, desc(Adjusted.P.value), desc(Combined.Score))
c3_bp <- c3_bp[1:10, ]
c3_bp$Term <- str_wrap(c3_bp$Term, 45)
c3_bp$Term <- str_to_title(c3_bp$Term, 45)
c3_bp$Term <- factor(c3_bp$Term,
                     levels = c3_bp$Term[order(c3_bp$Adjusted.P.value, c3_bp$Combined.Score, decreasing = FALSE)])
c3_bp_plot <- c3_bp %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c3_bp_plot
ggsave('./Figures/cluster_3_bp.tiff', height = 1200, width = 2200, units = 'px', dpi = 300)


combined_bp <- c1_bp_plot + c2_bp_plot + c3_bp_plot
combined_bp + plot_annotation(title = 'GO Biological Process Analysis of Identified Gene Clusters')
combined_bp + plot_annotation(tag_levels = 'A')
ggsave('./Figures/cluster_bp.tiff', height = 1400, width = 5000, units = 'px', dpi = 300)






c1_mf <- read.csv('./Results/cluster_1_go_mf.txt', sep = '\t')
c1_mf$Adjusted.P.value <- -log10(c1_mf$Adjusted.P.value)
c1_mf <- arrange(c1_mf, desc(Adjusted.P.value), desc(Combined.Score))
c1_mf <- c1_mf[1:10, ]
c1_mf$Term <- str_wrap(c1_mf$Term, 45)
c1_mf$Term <- str_to_title(c1_mf$Term, 45)
c1_mf$Term <- factor(c1_mf$Term,
                     levels = c1_mf$Term[order(c1_mf$Adjusted.P.value, c1_mf$Combined.Score, decreasing = FALSE)])
c1_mf_plot <- c1_mf %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c1_mf_plot
ggsave('./Figures/cluster_1_mf.tiff', height = 1200, width = 2200, units = 'px', dpi = 300)



c2_mf <- read.csv('./Results/cluster_2_go_mf.txt', sep = '\t')
c2_mf$Adjusted.P.value <- -log10(c2_mf$Adjusted.P.value)
c2_mf <- arrange(c2_mf, desc(Adjusted.P.value), desc(Combined.Score))
c2_mf <- c2_mf[1:10, ]
c2_mf$Term <- str_wrap(c2_mf$Term, 45)
c2_mf$Term <- str_to_title(c2_mf$Term, 45)
c2_mf$Term <- factor(c2_mf$Term,
                     levels = c2_mf$Term[order(c2_mf$Adjusted.P.value, c2_mf$Combined.Score, decreasing = FALSE)])
c2_mf_plot <- c2_mf %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c2_mf_plot
ggsave('./Figures/cluster_2_mf.tiff', height = 1200, width = 2200, units = 'px', dpi = 300)




c3_mf <- read.csv('./Results/cluster_3_go_mf.txt', sep = '\t')
c3_mf$Adjusted.P.value <- -log10(c3_mf$Adjusted.P.value)
c3_mf <- arrange(c3_mf, desc(Adjusted.P.value), desc(Combined.Score))
c3_mf <- c3_mf[1:10, ]
c3_mf$Term <- str_wrap(c3_mf$Term, 45)
c3_mf$Term <- str_to_title(c3_mf$Term, 45)
c3_mf$Term <- factor(c3_mf$Term,
                     levels = c3_mf$Term[order(c3_mf$Adjusted.P.value, c3_mf$Combined.Score, decreasing = FALSE)])
c3_mf_plot <- c3_mf %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c3_mf_plot
ggsave('./Figures/cluster_3_mf.tiff', height = 1200, width = 2200, units = 'px', dpi = 300)


combined_mf <- c1_mf_plot + c2_mf_plot + c3_mf_plot
combined_mf + plot_annotation(title = 'GO Molecular Function Analysis of Identified Gene Clusters')
combined_mf + plot_annotation(tag_levels = 'A')
ggsave('./Figures/cluster_mf.tiff', height = 1400, width = 5000, units = 'px', dpi = 300)







c1_cc <- read.csv('./Results/cluster_1_go_cc.txt', sep = '\t')
c1_cc$Adjusted.P.value <- -log10(c1_cc$Adjusted.P.value)
c1_cc <- arrange(c1_cc, desc(Adjusted.P.value), desc(Combined.Score))
c1_cc <- c1_cc[1:10, ]
c1_cc$Term <- str_wrap(c1_cc$Term, 45)
c1_cc$Term <- str_to_title(c1_cc$Term, 45)
c1_cc$Term <- factor(c1_cc$Term,
                     levels = c1_cc$Term[order(c1_cc$Adjusted.P.value, c1_cc$Combined.Score, decreasing = FALSE)])
c1_cc_plot <- c1_cc %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c1_cc_plot
ggsave('./Figures/cluster_1_cc.tiff', height = 1200, width = 2200, units = 'px', dpi = 300)



c2_cc <- read.csv('./Results/cluster_2_go_cc.txt', sep = '\t')
c2_cc$Adjusted.P.value <- -log10(c2_cc$Adjusted.P.value)
c2_cc <- arrange(c2_cc, desc(Adjusted.P.value), desc(Combined.Score))
c2_cc <- c2_cc[1:10, ]
c2_cc$Term <- str_wrap(c2_cc$Term, 45)
c2_cc$Term <- str_to_title(c2_cc$Term, 45)
c2_cc$Term <- factor(c2_cc$Term,
                     levels = c2_cc$Term[order(c2_cc$Adjusted.P.value, c2_cc$Combined.Score, decreasing = FALSE)])
c2_cc_plot <- c2_cc %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c2_cc_plot
ggsave('./Figures/cluster_2_cc.tiff', height = 1200, width = 2200, units = 'px', dpi = 300)



c3_cc <- read.csv('./Results/cluster_3_go_cc.txt', sep = '\t')
c3_cc$Adjusted.P.value <- -log10(c3_cc$Adjusted.P.value)
c3_cc <- arrange(c3_cc, desc(Adjusted.P.value), desc(Combined.Score))
c3_cc <- c3_cc[1:10, ]
c3_cc$Term <- str_wrap(c3_cc$Term, 45)
c3_cc$Term <- str_to_title(c3_cc$Term, 45)
c3_cc$Term <- factor(c3_cc$Term,
                     levels = c3_cc$Term[order(c3_cc$Adjusted.P.value, c3_cc$Combined.Score, decreasing = FALSE)])
c3_cc_plot <- c3_cc %>%
  ggplot(aes(y = Adjusted.P.value, x = Term)) +
  geom_bar(stat = 'identity') +
  geom_col(aes(fill = Adjusted.P.value)) +
  scale_fill_gradientn(colours = c("firebrick1", "firebrick3", "firebrick4")) +
  labs(fill='-log10(P Value)') +
  coord_flip() +
  ylab("-log10(Adjusted P Value)") + xlab("Term") +
  theme_pubclean() + theme(legend.position = 'none')
c3_cc_plot
ggsave('./Figures/cluster_3_cc.tiff', height = 1200, width = 2200, units = 'px', dpi = 300)



combined_cc <- c1_cc_plot + c2_cc_plot + c3_cc_plot
combined_cc + plot_annotation(title = 'GO Cellular Components Analysis of Identified Gene Clusters')
combined_cc + plot_annotation(tag_levels = 'A')
ggsave('./Figures/cluster_cc.tiff', height = 1400, width = 5000, units = 'px', dpi = 300)
