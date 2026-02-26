rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ggplot2")
library("dplyr")
library(pryr)
library("caret")
library(ggtree)
library(ape)
library(patchwork)
library(ggtext)
library(ggpubr)
library(multcompView)
library(tidytext)
library(lme4)
library(broom.mixed)


#import data
ISG_meta <- read.table("Aves_Mam_mbio_250613_ISGscore_logan.txt", header = T, sep = "\t")

####species level plot#######
species_inf_summary <- ISG_meta %>%
  ungroup() %>%
  filter(Infection %in% c("Positive", "Negative")) %>% 
  group_by(species, Infection) %>%
  summarise(n = n(), .groups = "drop") %>%
  tidyr::pivot_wider(
    names_from = Infection,
    values_from = n,
    values_fill = 0
  ) %>%
  mutate(
    total = Positive + Negative,
    positive_ratio = Positive / total*100
  ) %>%
  arrange(desc(total))

species_inf_summary

sp_class_info <- ISG_meta %>%
  ungroup() %>%
  select(species, class) %>%
  distinct()

# 2. 結合して class を追加
species_inf_summary <- species_inf_summary %>%
  left_join(sp_class_info, by = "species")

topN <- 15
top_sp <- species_inf_summary %>%
  arrange(desc(total)) %>%
  slice_head(n = topN) %>%
  pull(species)

top_part_sp <- species_inf_summary %>%
  filter(species %in% top_sp) %>%
  arrange(desc(total))

top_part_sp <- top_part_sp %>%
  mutate(species = factor(species, levels = rev(top_sp)))

p_sp_num <- ggplot(top_part_sp, aes(x = species, y = total, fill = class)) +
  geom_col() +
  coord_flip() +
  labs(x = "Order", y = "Number of samples") +
  theme_classic(base_size = 10)
p_sp_num
p_sp_ratio <- ggplot(top_part_sp, aes(x = "positive_ratio", y = species, fill = positive_ratio)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  scale_fill_gradient(
    low = "white",
    high = "darkgreen",
    limits = c(0, 50),
    name = "Positive ratio (%)"
  ) +
  theme_minimal(base_size = 11) +
  labs(
    x = NULL,
    y = "Order",
    title = "Positive ratio per host order"
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

combined_sp_plot <- p_sp_ratio + p_sp_num + plot_layout(ncol = 2, widths = c(1, 3))
combined_sp_plot

output_dir <- "."
sp_bar <- paste0(output_dir, "/", "sp15z_num_barplot_260115.pdf")
#ggsave(sp_bar, plot = combined_sp_plot, height=7, width = 7)

#positivity rate top 15
p_topN <- 15
p_top_sp_10 <- species_inf_summary %>%
  filter(total >= 30) %>%
  arrange(desc(positive_ratio)) %>%
  slice_head(n = p_topN) %>%
  pull(species)

p_top_part_sp_10 <- species_inf_summary %>%
  mutate(log = log10(total)) %>%
  filter(species %in% p_top_sp_10) %>%
  arrange(desc(total))

p_top_part_sp_10 <- p_top_part_sp_10 %>%
  mutate(species = fct_reorder(species, positive_ratio, .desc = FALSE))

p_p_sp_rate_10 <- ggplot(p_top_part_sp_10, aes(x = positive_ratio, y = species, fill = class)) +
  geom_col() +
  labs(x = "Infection ratio", y = NULL, title = "Positive ratio per host order") +
  theme_classic(base_size = 10) +
  theme(axis.text.x = element_text())

p_p_sp_10_num <- ggplot(p_top_part_sp_10, aes(x = "total", y = species, fill = total)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  scale_fill_gradient(
    low = "white",
    high = "red",
    oob = scales::squish,
    limits = c(0, 300),
    name = "sample size"
  ) +
  theme_minimal(base_size = 11) +
  labs(
    x = NULL,
    y = NULL,
    title = NULL
  ) +
  theme(
    axis.text.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid = element_blank()
  )

p_combined_sp_10_plot <- p_p_sp_rate_10 + p_p_sp_10_num + plot_layout(widths  = c(5, 1))
p_combined_sp_10_plot

output_dir <- "."
p_combined_sp_10 <- paste0(output_dir, "/", "sp_inf10_infratio_260116.pdf")
#ggsave(p_combined_sp_10, plot = p_combined_sp_10_plot, height=7, width = 7)

#BioProject_IDの偏り
bp_counts <- ISG_meta %>%
  group_by(BioProject_ID) %>%
  summarise(n_samples = n())

bp_cumu <- bp_counts %>%
  arrange(desc(n_samples)) %>%
  mutate(
    rank = row_number(),
    cum_sum = cumsum(n_samples),
    cum_frac = cum_sum / sum(n_samples) * 100
  )

bp_bins <- bp_counts %>%
  mutate(bin = cut(n_samples,
                   breaks = c(0, 10, 50, 100, 500, 1000, Inf),
                   labels = c("1–10", "11–50", "51–100", "101–500", "501–1000", "1000+"),
                   right = TRUE)) %>%
  count(bin)

#plot
bp_bins_p <- ggplot(bp_bins, aes(x = bin, y = n)) +
  geom_col(fill = "gray50") +
  labs(
    title = "Number of BioProject",
    x = "Number of Samples",
    y = "Number of BioProject"
  ) +
  theme_classic(base_size = 12) +
  theme(
      axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
      panel.grid = element_blank()
  )
plot(bp_bins_p)
output_dir <- "."
bp_bins_p_file <- paste0(output_dir, "/", "BioProject_bin.pdf")
#ggsave(bp_bins_p_file, plot = bp_bins_p, height=5, width = 5)
