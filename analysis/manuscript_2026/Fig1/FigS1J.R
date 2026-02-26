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
library(viridis)
library(ggpubr)

ISG_cleaned_isgscore <- read.table("tool5_260217.txt", header = T, sep = "\t")

#ロング形式に変換
ISG_long_score <- ISG_cleaned_isgscore %>%
  group_by(Induction, species, hum_symbol) %>%
  summarise(RPM_mean = mean(RPM),
            kallisto_mean = mean(ISG_score_kallisto),
            kma_mean = mean(ISG_score_kma),
            salmon_mean = mean(ISG_score_salmon),
            star_mean = mean(ISG_score_star),
            bowtie_mean = mean(ISG_score_bowtie)
  ) %>%
  ungroup() %>%
  select(hum_symbol, species, Induction, kallisto_mean, kma_mean, salmon_mean, star_mean, bowtie_mean, RPM_mean) %>%
  pivot_longer(cols = c(kallisto_mean, kma_mean, salmon_mean, star_mean, bowtie_mean), names_to = "method", values_to = "normalized_value") %>%
  mutate(method = str_remove(method, "^norm_"))


ISG_score_diff_df <- ISG_long_score %>%
  group_by(hum_symbol, species, method, Induction) %>%
  summarise(normalized_value = mean(normalized_value, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(
    names_from = Induction,
    values_from = normalized_value
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



ISG_zscore_log2FC <- ISG_score_diff_df %>%
  full_join(ISG_log2FC_mean, by = c("species", "hum_symbol"))

manual_order <- c("Hsapiens", "Rnorvegicus", "Btaurus", "Oaries", "Sscrofa", "Ecaballus", "Clupus", "Pvampyrus", "Mlucifugus", "Ggallus")
species_order <- c("Gallus gallus", "Homo sapiens", "Rattus norvegicus", "Myotis lucifugus",
                   "Pteropus vampyrus","Sus scrofa", "Ovis aries", "Bos taurus",
                   "Equus caballus", "Canis lupus")
# ✅ factor() ISG_score_10sp_n
ISG_zscore_log2FC <- ISG_zscore_log2FC %>%
  mutate(species = factor(species, levels = species_order)) %>%
  filter(!is.na(ISG_score_diff))

cor_df <- ISG_zscore_log2FC %>%
  group_by(species, method) %>%
  summarise(cor = cor(log2FC_mean, ISG_score_diff, method = "pearson", use = "complete.obs"), .groups = "drop")

method_order <- cor_df %>%
  group_by(method) %>%
  summarise(mean_cor = mean(cor, na.rm = TRUE)) %>%
  arrange(mean_cor) %>%
  pull(method)

# 2. factor 順序を設定（平均相関の高い順）
cor_df$method <- factor(cor_df$method, levels = method_order)

# 1. Inductionごとにラベルを作成（Positive → Negative の順で並べたい）
species_order <- c("Gallus gallus", "Homo sapiens", "Rattus norvegicus", "Myotis lucifugus",
                   "Pteropus vampyrus","Sus scrofa", "Ovis aries", "Bos taurus",
                   "Equus caballus", "Canis lupus")


cor_heatmap <- ggplot(cor_df, aes(x = species, y = method, fill = cor)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(option = "D", name = "Pearson\nCorrelation", limits = c(0, 1)) +
  theme_classic(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.title = element_blank()
  ) +
  labs(
    title = "Correlation between log2FC and ISG Score",
    subtitle = "Across species and tools"
  )

print(cor_heatmap)

output_dir <- "."
file_name <- paste0(output_dir, "/", "tools5_comparison_260131.pdf")
#ggsave(file_name, plot = cor_heatmap, width = 5, height = 3)

