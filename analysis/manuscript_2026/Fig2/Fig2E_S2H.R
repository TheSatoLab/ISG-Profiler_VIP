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

#import data
ISG_meta <- read.table("Aves_Mam_mbio_250613_ISGscore_logan.txt", sep = "\t", header = T)
Aves_Mam_mbio_meta_filt <- read.table("Aves_Mam_mbio_metadata_ISG_logan_ML_250616_2.txt", sep = "\t", header = T)
Aves_Mam_mbio_genomad_filt <- read.table("Aves_Mam_mbio_genomad_filt_ML_250630.txt", sep = "\t", header = T)

ISG_virus <- ISG_meta %>%
  filter(Infection == "Positive") %>%
  select(-ISG_level) %>%
  inner_join(Aves_Mam_mbio_genomad_filt, by = c("ID"), relationship = "many-to-many") %>%
  select(-n_contig)

ISG_neg <- ISG_meta %>%
  filter(Infection == "Negative") %>%
  mutate(taxonomy = "Negative", nucleotide = "Negative")

ISG_virus_neg <- rbind(ISG_virus, ISG_neg)

#taxonomyごと
med_levels <- ISG_virus_neg %>%
  group_by(taxonomy) %>%
  summarize(med_ISG_med = median(ISG_mean, na.rm = TRUE)) %>%
  arrange(desc(med_ISG_med)) %>%
  pull(taxonomy)

ISG_virus_neg <- ISG_virus_neg %>%
  mutate(taxonomy = factor(taxonomy, levels = med_levels))

ISG_virus_neg_mod <- ISG_virus_neg %>%
  filter(taxonomy != "Negative") %>%
  separate(col = taxonomy, into=c("Virus", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";") %>%
  filter(Family != "")

med_levels <- ISG_virus_neg_mod %>%
  group_by(Family) %>%
  summarize(med_ISG_med = median(ISG_mean, na.rm = TRUE)) %>%
  arrange(desc(med_ISG_med)) %>%
  pull(Family) 

ISG_virus_neg_mod <- ISG_virus_neg_mod %>%
  mutate(Family = factor(Family, levels = med_levels))

################全てのサンプル########################
# 各ペアの比較＋p値計算
kruskal.test(ISG_mean ~ nucleotide, data = ISG_virus_neg)
anova_result <- aov(ISG_mean ~ nucleotide, data = ISG_virus_neg)
summary(anova_result)

# 並び順を指定
group_levels <- c("Negative", "DNA", "RNA")

# 並び順をfactorにして固定
ISG_virus_neg <- ISG_virus_neg %>%
  mutate(nucleotide = factor(nucleotide, levels = group_levels))

# ANOVAとTukeyHSD
anova_result <- aov(ISG_mean ~ nucleotide, data = ISG_virus_neg)
tukey_result <- TukeyHSD(anova_result)

# グループ記号の抽出
tukey_letters <- multcompLetters4(anova_result, tukey_result)
group_labels <- as.data.frame.list(tukey_letters$nucleotide)
group_labels$nucleotide <- rownames(group_labels)
colnames(group_labels)[1] <- "label"

# ラベルを配置するためのy位置
label_df <- ISG_virus_neg %>%
  group_by(nucleotide) %>%
  summarise(y_pos = max(ISG_mean, na.rm = TRUE) + 0.3) %>%
  left_join(group_labels, by = "nucleotide")

# TukeyHSDからペアごとのp値抽出 → シンボルへ変換
pval_lines <- as.data.frame(tukey_result$nucleotide) %>%
  rownames_to_column(var = "comparison") %>%
  separate(comparison, into = c("group1", "group2"), sep = "-") %>%
  mutate(
    signif_label = case_when(
      `p adj` < 0.001 ~ "***",
      `p adj` < 0.01  ~ "**",
      `p adj` < 0.05  ~ "*",
      TRUE            ~ "ns"
    ),
    x1 = match(group1, group_levels),
    x2 = match(group2, group_levels),
    x_mid = (x1 + x2) / 2,
    y.position = max(ISG_virus_neg$ISG_mean, na.rm = TRUE) + seq(0.6, 1.2, length.out = n())
  )

# 最終プロット
DNARNAneg_plot<- ggplot(ISG_virus_neg, aes(x = nucleotide, y = ISG_mean, fill = nucleotide)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  geom_segment(data = pval_lines,
               aes(x = x1, xend = x2,
                   y = y.position/2, yend = y.position/2),
               inherit.aes = FALSE, linewidth = 0.6) +
  geom_text(data = pval_lines,
            aes(x = x_mid, y = y.position/2*1.005, label = signif_label),
            inherit.aes = FALSE, size = 5) +
  coord_cartesian(ylim = c(-2, 2.5)) +
  scale_fill_manual(values = c(
    "DNA" = "#619CFF",
    "RNA" = "#F8766D",
    "Negative" = "#BBBBBB"
  )) +
  theme_classic() +
  labs(title = "ISG score comparison with TukeyHSD (groups & significance), Mammal&Aves",
       x = "Infection condition/Nucleotide type", y = "ISG score")

plot(DNARNAneg_plot)
output_dir <- "."
file_name <- paste0(output_dir, "/", "DNARNAneg_zscore_box_260110.pdf")
#ggsave(file_name, plot = DNARNAneg_plot, width = 4, height = 4)


#top10でheatmap作成
species_freq <- ISG_virus_neg %>%
  group_by(species, nucleotide) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = nucleotide, values_from = n, values_fill = 0) %>%
  filter(DNA >= 15, RNA >= 15, Negative >= 15)

species_freq_top10 <- ISG_virus_neg %>%
  group_by(species, nucleotide) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = nucleotide, values_from = n, values_fill = 0) %>%
  mutate(total = DNA + RNA + Negative) %>%        
  slice_max(total, n = 10)                        

# 対象種だけ抽出
ISG_virus_neg_multi_10 <- ISG_virus_neg %>%
  filter(species %in% species_freq_top10$species)

ISG_virus_neg_multi_10_mean <- ISG_virus_neg_multi_10 %>%
  group_by(species, nucleotide) %>%
  summarise(ISG_mean_mean = mean(ISG_mean, na.rm = TRUE), .groups = "drop") %>%
  ungroup()

# lineplot作成
p_line_tax <- ggplot(ISG_virus_neg_multi_10_mean,
                     aes(x = nucleotide, y = ISG_mean_mean, group = species, color = species)) +
  geom_line(linewidth = 1) +
  geom_point(size = 2) +
  theme_classic(base_size = 12) +
  labs(
    x = "Infection status",
    y = "Mean ISG score",
    title = "Mean ISG score across infection status (Top species)"
  ) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    panel.grid.minor = element_blank()
  )
p_line_tax
output_dir <- "."
p_line_tax_name <- paste0(output_dir, "/", "DNARNAneg_line_260116.pdf")
#ggsave(p_line_tax_name, plot = p_line_tax, width = 5, height = 4)

