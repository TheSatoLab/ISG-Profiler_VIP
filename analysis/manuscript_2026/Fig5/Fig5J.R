library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(scales)

# change working directory
setwd("/Users/kyokokurihara/iLab/itolab_backup/backup-latest/Lab/projects/2507blastx/output/250909_4474_samples/final_test/chaphama_260206_MSigDB_Hallmark_2020_stat_blastx_filtered/dotplot_chaphama/")

# read table
df_all <- read_csv("NES_termshared_top10_FDR05.csv") |>
  mutate(
    Term = gsub("-", "_", Term),
    sig  = if_else(`FDR q-val` < 0.05, "FDR < 0.05", "ns"),
    sig  = factor(sig, levels = c("FDR < 0.05", "ns"))
  )

# term order
term_levels <- df_all |>
  distinct(Term) |>
  pull(Term)

# project order
project_levels <- df_all |>
  distinct(Project) |>
  pull(Project)

df_all <- df_all |>
  mutate(
    Term    = factor(Term, levels = term_levels),
    Project = factor(Project, levels = project_levels)
  )

global_size_max <- max(abs(df_all$mlog10_fdr), na.rm = TRUE)

p <- ggplot(df_all, aes(
  x = Project, y = Term,
  size = abs(mlog10_fdr),
  fill = NES,
  color = sig
)) +
  geom_point(shape = 21, stroke = 0.6) +
  scale_fill_gradient2(
    low = "#2166AC", mid = "#F7F7F7", high = "#B2182B",
    midpoint = 0, limits = c(-3, 3), oob = squish, name = "NES"
  ) +
  scale_size_continuous(
    limits = c(0, global_size_max),
    range = c(1, 6),
    breaks = pretty(c(0, global_size_max), n = 3),
    name = "-log10(FDR)"
  ) +
  scale_color_manual(
    values = c("FDR < 0.05" = "black", "ns" = "white"),
    breaks = c("FDR < 0.05"),
    name = "FDR"
  ) +
  labs(x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(color = "black", angle = 90, hjust = 1, vjust = 1),
    axis.text.y = element_text(color = "black"),
    legend.title = element_text(size = 9),
    legend.text  = element_text(size = 8),
    plot.margin = margin(5, 5, 5, 5),
    panel.grid.major = element_line(linewidth = 0.3),  # TODO
    panel.grid.minor = element_line(linewidth = 0.3)  # TODO
  ) +
  guides(
    fill = guide_colorbar(barheight = unit(20, "mm"), barwidth = unit(4, "mm"))
  )

ggsave("stat_dotplot_1panel.pdf", p, width = 4.2, height = 6.5, units = "in")