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

ISGcntl_norm.list <- read.table("sp10_260216_forin.txt", header = T, sep = "\t")
ISG_meanscore_sp <- ISGcntl_norm.list %>%
  filter(type == "ISG") %>%
  ungroup() %>%
  group_by(type, species, Induction, hum_symbol) %>%
  summarise(ISG_score_mean = mean(ISG_score, na.rm = TRUE)) %>%
  ungroup() %>%
  select(species, hum_symbol, Induction, ISG_score_mean) %>%
  pivot_wider(
    names_from = Induction,
    values_from = ISG_score_mean
  ) %>%
  mutate(ISG_score_diff = Positive - Negative)

ISGcntl_log2FC <- read.table("isgcntl_241003_log2FC.list", header = T, sep = "\t") %>%
  mutate(species = dplyr::recode(species,
                                 "Hsapiens" = "Homo sapiens",
                                 "Ggallus" = "Gallus gallus",
                                 "Clupus" = "Canis lupus",
                                 "Btaurus" = "Bos taurus",
                                 "Ecaballus" = "Equus caballus",
                                 "Mlucifugus" = "Myotis lucifugus",
                                 "Oaries" = "Ovis aries",
                                 "Pvampyrus" = "Pteropus vampyrus",
                                 "Rnorvegicus" = "Rattus norvegicus",
                                 "Sscrofa" = "Sus scrofa"
  ))

ISG_log2FC_mean <- ISGcntl_log2FC %>%
  filter(type == "ISG") %>%
  group_by(type, species, hum_symbol) %>%
  summarise(log2FC_mean = mean(log2FoldChange)) %>%
  ungroup() %>%
  select(-type)  %>%
  mutate(
    log2FC_mean = if_else(species == "Pteropus vampyrus", -log2FC_mean, log2FC_mean)
  )

ISG_zscore_log2FC <- ISG_meanscore_sp %>%
  full_join(ISG_log2FC_mean, by = c("species", "hum_symbol"))

manual_order <- c("Hsapiens", "Rnorvegicus", "Btaurus", "Oaries", "Sscrofa", "Ecaballus", "Clupus", "Pvampyrus", "Mlucifugus", "Ggallus")
species_order <- c("Gallus gallus", "Homo sapiens", "Rattus norvegicus", "Myotis lucifugus",
                   "Pteropus vampyrus","Sus scrofa", "Ovis aries", "Bos taurus",
                   "Equus caballus", "Canis lupus")
# ✅ factor() ISG_score_10sp_n
ISG_zscore_log2FC <- ISG_zscore_log2FC %>%
  mutate(species = factor(species, levels = species_order))

library(ggplot2)
library(ggpubr)

log2_zscore_p <- ggplot(ISG_zscore_log2FC, aes(x = log2FC_mean, y = ISG_score_diff)) +
  geom_point(alpha = 0.8, size = 2) +
  geom_smooth(method = "lm", se = FALSE, color = "gray40", linetype = "solid") +  
  geom_vline(xintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +  
  geom_hline(yintercept = 0, linetype = "dotted", color = "black", linewidth = 0.5) +  
  stat_cor(method = "pearson", label.x = 0, label.y = 1.5, size = 3) +
  coord_cartesian(xlim = c(0, 16), ylim = c(-0.4, 4.2)) +
  facet_wrap(~ species, nrow = 1, scales = "free") +
  theme_classic(base_size = 12) +
  labs(
    x = "Mean log2 Fold Change, STAR",
    y = "Z-score based ISG score difference",
    title = "Comparison of ISG Expression: log2FC vs. Z-score ISG Score"
  )

plot(log2_zscore_p)
output_dir <- "."
log2_zscore_p_name <- paste0(output_dir, "/", "sp10_ISG_log2FC_Zscoredif.pdf")
#ggsave(log2_zscore_p_name, plot = log2_zscore_p, width = 12, height = 2)

ISGscore_ISGcntl_mean.list <- ISGcntl_norm.list %>%
  group_by(ID, type) %>%
  mutate(ISG_score_mean = mean(ISG_score, na.rm = TRUE)) %>%
  dplyr::select(ID, ISG_score_mean, species, Induction, type, cntl_sum, taxid) %>%
  unique()

ISGscore_ISG_mean.list <- ISGscore_ISGcntl_mean.list %>%
  filter(type == "ISG") %>%
  mutate(BioProject_ID = "PRJEB21332") %>%
  select(ID, BioProject_ID, taxid,species,Induction,cntl_sum, ISG_score_mean)

#write.table(ISGscore_ISG_mean.list, "sp10_ISG.txt", row.names = F, sep = "\t", quote = F)

wilcox_results <- ISGscore_ISG_mean.list %>%
  group_by(species) %>%
  summarise(
    p_value = wilcox.test(ISG_score_mean ~ Induction)$p.value,
    .groups = "drop"
  )

library(ggpubr)

# 表示したい順序をベクトルで定義
species_order <- c("Gallus gallus", "Homo sapiens", "Rattus norvegicus", "Myotis lucifugus",
                   "Pteropus vampyrus","Sus scrofa", "Ovis aries", "Bos taurus",
                   "Equus caballus", "Canis lupus")

# species 列を順序付き factor に変換
ISGscore_ISG_mean.list$species <- factor(ISGscore_ISG_mean.list$species, levels = species_order)

sp10_z_plot_ISG <- ggplot(ISGscore_ISG_mean.list, aes(x = Induction, y = ISG_score_mean, fill = Induction)) +
  # geom_jitter(width = 0.2, alpha = 0.5) +
  # geom_boxplot(outlier.shape = NA, alpha = 1.0) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
  geom_jitter(aes(color = Induction), shape = 21, color = "black", width = 0.2, alpha = 1.0, size = 1.5) + 
  stat_compare_means(aes(group = Induction), method = "t.test", label = "p.signif") +
  theme_classic() +
  facet_grid(. ~ species) +
  theme(axis.text.x = element_text()) +
  labs(title = "ISG_score by Induction status", x = NULL, y = "ISG Score") +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) +
  scale_x_discrete(
    limits = c("Negative", "Positive"),
    labels = c("Negative" = "-", "Positive" = "+") 
  )
plot(sp10_z_plot_ISG)

output_dir <- "."
file_name <- paste0(output_dir, "/", "sp10_salmon_meanISG_onlyISG_260115.pdf")
#ggsave(file_name, plot = sp10_z_plot_ISG, width = 8, height = 5)
