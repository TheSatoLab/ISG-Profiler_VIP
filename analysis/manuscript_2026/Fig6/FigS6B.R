library(ComplexHeatmap)
library(circlize)
library(tidyverse)　## 
grid::unit

gs.interest.v <- c("Interferon Alpha Response",
                   "IL_6/JAK/STAT3 Signaling",
                   "Xenobiotic Metabolism",
                   "Adipogenesis")

project_Id.v <- c("PRJNA774885","PRJEB63475")


# change working directory
setwd("/Users/kyokokurihara/iLab/itolab_backup/backup-latest/Lab/projects/2507blastx/output/250909_4474_samples/final_test/hepatovirus_260113_MSigDB_Hallmark_2020_stat/heatmap_hepatovirus/")

# read tables
mat_df  <- read.csv("gene_project_log2FC_matrix.csv", row.names = 1, check.names = FALSE)
term_bin_df <- read.csv("term_anno.csv", row.names = 1, check.names = FALSE)

mat_df.long <- mat_df %>% mutate(symbol = rownames(mat_df)) %>% gather(key = project_Id, value = log2FC, -symbol)
term_bin_df.long <- term_bin_df %>% mutate(symbol = rownames(mat_df)) %>% gather(key = gs, value = membership, -symbol) %>% filter(membership == 1)

term_bin_df.long <- term_bin_df.long %>% filter(gs %in% gs.interest.v) %>% mutate(gs = factor(gs, levels = gs.interest.v)) %>% arrange(gs, symbol)

symbol_order.v <- term_bin_df.long %>% pull(symbol) %>% unique()

term_bin_df.long <- term_bin_df.long %>% mutate(symbol = factor(symbol, levels = rev(symbol_order.v)))



mat_df.long <- mat_df.long %>% filter(symbol %in% symbol_order.v) %>%
  mutate(symbol = factor(symbol, levels = rev(symbol_order.v)),
         project_Id = factor(project_Id, levels = project_Id.v))


col_fun <- colorRamp2(c(-2, 0, 2), c("#1B9E9E", "white", "#C51B7D"))

g1 <- ggplot(mat_df.long, aes(x = project_Id, y = symbol, fill = log2FC)) +
  geom_tile() +
  scale_fill_gradientn(
    colours = col_fun(seq(-2, 2, length.out = 100)),
    limits = c(-2, 2),
    oob = scales::squish
  ) +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    )
  )

g2 <- ggplot(term_bin_df.long, aes(x = gs, y = symbol, fill = membership)) +
  geom_tile() +
  coord_fixed(ratio = 1) +
  theme_classic() +
  theme(
    axis.text.x = element_text(
      angle = 90,
      hjust = 1,
      vjust = 0.5
    ),
    panel.grid.major = element_line(color = "grey90", linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

pdf("stat_log2FC_heatmap.pdf", width = 6, height = 6)
g1 + g2
dev.off()