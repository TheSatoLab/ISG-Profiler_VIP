rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library(emmeans)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)
library(caret)
library(ggtree)
library(broom)
library(ape)
library(patchwork)
library(ggtext)
library(ggplot2)
library(multcompView)
library(scales)
library(car)
library(effectsize)


#import data
ISG_meta <- read.table("Aves_Mam_mbio_250613_ISGscore_logan.txt", sep = "\t", header = T)
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

ISG_virus_neg_mod <- ISG_virus_neg %>%
  separate(col = taxonomy, into=c("Virus", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";") %>%
  mutate(
    Virus = ifelse(Infection == "Negative" & is.na(Virus), "Negative", Virus),
    Realm = ifelse(Infection == "Negative" & is.na(Realm), "Negative", Realm),
    Kingdom = ifelse(Infection == "Negative" & is.na(Kingdom), "Negative", Kingdom),
    Phylum = ifelse(Infection == "Negative" & is.na(Phylum), "Negative", Phylum),
    Class = ifelse(Infection == "Negative" & is.na(Class), "Negative", Class),
    Order = ifelse(Infection == "Negative" & is.na(Order), "Negative", Order),
    Family = ifelse(Infection == "Negative" & is.na(Family), "Negative", Family),
    Family = ifelse(Infection == "Positive" & nucleotide == "DNA" & Family == "", "DNA Unclassified", Family),
    Family = ifelse(Infection == "Positive" & nucleotide == "RNA" & Family == "", "RNA Unclassified", Family)
  ) %>%
  ungroup()

#factorでNegativeを0にする
ISG_virus_neg$taxonomy <- relevel(factor(ISG_virus_neg$taxonomy), ref = "Negative")
lm_result <- lm(ISG_mean ~ taxonomy + species, data = ISG_virus_neg)
summary(lm_result)

vir_tax <- ISG_virus_neg %>%
  ungroup() %>%
  select(taxonomy, nucleotide, Infection) %>%
  unique()

lm_summary <- lm_result %>%
  tidy() %>%
  mutate(p_adj = p.adjust(p.value, method = "BH")) %>%
  rename(taxonomy = term) %>%
  mutate(taxonomy = str_replace(taxonomy, "taxonomy", "")) %>%
  inner_join(vir_tax, by = "taxonomy") %>%
  mutate(taxonomy1 = taxonomy) %>%
  separate(col = taxonomy, into=c("Virus", "Realm", "Kingdom", "Phylum", "Class", "Order", "Family"), sep = ";") %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  mutate(
    Virus = ifelse(Infection == "Negative" & is.na(Virus), "Negative", Virus),
    Realm = ifelse(Infection == "Negative" & is.na(Realm), "Negative", Realm),
    Kingdom = ifelse(Infection == "Negative" & is.na(Kingdom), "Negative", Kingdom),
    Phylum = ifelse(Infection == "Negative" & is.na(Phylum), "Negative", Phylum),
    Class = ifelse(Infection == "Negative" & is.na(Class), "Negative", Class),
    Order = ifelse(Infection == "Negative" & is.na(Order), "Negative", Order),
    Family = ifelse(Infection == "Negative" & is.na(Family), "Negative", Family),
    Family = ifelse(Infection == "Positive" & nucleotide == "DNA" & is.na(Family), "DNA Unclassified", Family),
    Family = ifelse(Infection == "Positive" & nucleotide == "RNA" & is.na(Family), "RNA Unclassified", Family)
  ) %>%
  ungroup() %>%
  filter(Family != "RNA Unclassified" & Family != "DNA Unclassified")

#write.table(lm_summary, "lm_virus_250522.list", row.names = F, sep = "\t", quote = F)

# 並び順を作成（estimate で並べる）
dot_family_order <- lm_summary %>%
  arrange(estimate) %>%
  pull(Family)

# 並び順に従って factor に変換
lm_summary <- lm_summary %>%
  mutate(Family = factor(Family, levels = dot_family_order))

em_plot <- ggplot(lm_summary, aes(x = estimate, y = Family, color = nucleotide)) +
  geom_point(size = 3) +
  geom_errorbarh(aes(xmin = estimate - std.error,
                     xmax = estimate + std.error), height = 0.2) +
  coord_flip() +
  theme_classic() +
  labs(
    title = "Estimated ISG score by Virus Family (adjusted for Host species), linear model",
    x = "Estimated ISG_score",
    y = NULL
  ) +
  scale_color_manual(values = c("DNA" = "#619CFF", "RNA" = "salmon", "Negative" = "gray")) +
  theme(
    panel.grid = element_blank(),
    axis.text.x = element_blank(),
  )
plot(em_plot)

inf_num <- ISG_virus_neg_mod %>%
  ungroup() %>%
  group_by(species, order, Family, nucleotide) %>%
  summarise(n_sample = n()) %>%
  ungroup() %>%
  group_by(Family, nucleotide) %>%
  summarise(n_sp = n(), n_sample_total = sum(n_sample)) %>%
  ungroup() %>%
  mutate(log10_sample_number = log10(n_sample_total), log10_species_number = log10(n_sp)) %>%
  filter(Family %in% lm_summary$Family) %>%
  mutate(Family = factor(Family, levels = dot_family_order))


p_species <- ggplot(inf_num, aes(x = Family, y = "log10_species_number", fill = log10_species_number)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  scale_fill_gradient(low = "white", high = "blue", name = "log10 species", limits = c(-0.3, 3)) +
  theme_minimal(base_size = 10) +
  labs(x = NULL, y = NULL) +
  theme(
    axis.text.x = element_blank(),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

p_samples <- ggplot(inf_num, aes(x = Family, y = "log10_sample_number", fill = log10_sample_number)) +
  geom_tile(color = "white", width = 0.9, height = 0.9) +
  scale_fill_gradient(low = "white", high = "red", name = "log10 samples", limits = c(-0.4, 5.3)) +
  theme_minimal(base_size = 10) +
  labs(x = "Virus Family", y = NULL) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    axis.ticks = element_blank(),
    panel.grid = element_blank()
  )

p_species
p_samples

em_hm_plot <- em_plot + p_species + p_samples +
  plot_layout(heights = c(8, 1, 1), guides = "collect")
plot(em_hm_plot)

output_dir <- "./"
em_hm_name <- paste0(output_dir, "/", "lm_family_hm_260110.pdf")
#ggsave(em_hm_name, plot = em_hm_plot, width = 8, height = 5)


lm_summary_filt <- lm_summary %>%
  filter(p_adj < 0.01 & estimate > 0.1)
virus_filt_list <- lm_summary_filt %>%
  select(taxonomy1) %>%
  rename(taxonomy = taxonomy1)

lm_summary_rm <- lm_summary %>%
  filter(p_adj >= 0.01 | estimate <= 0.1) %>%
  select(taxonomy1) %>%
  rename(taxonomy = taxonomy1)
#write.table(lm_summary_rm, "low_ISG_virus_260110.txt", row.names = F, sep = "\t", quote = F)

#per virus ISG expression analysis
target_families <- c("Bornaviridae", "Togaviridae", "Flaviviridae", "Peribunyaviridae",
                     "Parvoviridae", "Anelloviridae", "Smacoviridae")

# 2. Positiveサンプルを取得
positive_samples <- ISG_virus_neg_mod %>%
  filter(Infection == "Positive", Family %in% target_families)

# 3. Positiveサンプルから、各ウイルスファミリーとspeciesの組み合わせを取得
family_species_pairs <- positive_samples %>%
  select(Family, species) %>%
  distinct()

# 4. Negativeサンプルから、上記speciesのデータを抽出し、Familyごとに割り当て（= speciesごとにfamilyを複製）
negative_samples <- ISG_virus_neg_mod %>%
  filter(Infection == "Negative") %>%
  select(-Family)

# 5. Family情報と結合して、正しいNegativeペアを作る
negative_expanded <- family_species_pairs %>%
  left_join(negative_samples, by = "species") %>%
  mutate(Family = Family, Infection = "Negative") %>%
  select(colnames(positive_samples))  # 列順を positive_samples に揃える

# 6. PositiveとNegativeの結合
combined_expanded <- bind_rows(positive_samples, negative_expanded) %>%
  filter(!is.na(ISG_mean))

my_order <- c("Togaviridae", "Flaviviridae", "Peribunyaviridae", "Bornaviridae", "Parvoviridae", "Smacoviridae", "Anelloviridae")

combined_expanded$Family <- factor(combined_expanded$Family, levels = my_order)

low_high_boxplot <- ggplot(combined_expanded, aes(x = Infection, y = ISG_mean, fill = Infection)) +
  geom_boxplot(outlier.shape = NA) +
  #geom_jitter(width = 0.2, alpha = 0.5, size = 0.8) +
  facet_wrap(~ Family, scales = "fixed", nrow = 1) +  # 横一列に並べる
  stat_compare_means(method = "t.test", label = "p.signif") +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) +
  theme_classic(base_size = 12) +
  labs(
    title = "ISG Score Comparison by Virus Family",
    x = "Infection Status",
    y = "ISG Score"
  ) +
  coord_cartesian(ylim = c(-2, 2))

output_dir <- "./"
low_high_boxplot_name <- paste0(output_dir, "/", "low_high_virus_260110.pdf")
#ggsave(low_high_boxplot_name, plot = low_high_boxplot, width = 6, height = 4)

#line plot with top 10 species in each Family
top10_species_per_family <- combined_expanded %>%
  filter(Infection == "Positive") %>%
  group_by(Family, species) %>%
  summarise(n = n(), .groups = "drop") %>%
  arrange(Family, desc(n)) %>%
  group_by(Family) %>%
  slice_head(n = 10) # 各ファミリーで上位10種を抽出

#loopでplot分ける
# 2. Familyリストを取得
target_families <- unique(top10_species_per_family$Family)

library(lme4)
library(lmerTest)  # p-value 出すために便利

# 3. Familyごとにプロットを作成
plot_list <- list()
lmer_results_df <- data.frame()
all_plot_data <- data.frame()
for (fam in target_families) {
  
  ## 1. Family × species × Infection のサンプル数
  species_counts <- combined_expanded %>%
    filter(Family == fam) %>%
    count(species, Infection) %>%
    tidyr::pivot_wider(
      names_from = Infection,
      values_from = n,
      values_fill = 0
    )
  
  ## 2. Positive >= 3 & Negative >= 3 を満たす species
  valid_species <- species_counts %>%
    filter(Positive >= 3, Negative >= 3)
  
  if (nrow(valid_species) == 0) next  # 念のため
  
  ## 3. その中で Positive 数が多い Top10 species
  top10_species <- valid_species %>%
    filter(!str_ends(species, "_sp")) %>%
    arrange(desc(Positive)) %>%
    slice_head(n = 10) %>%
    pull(species)
  
  ## 4. プロット用データ作成
  data_fam <- combined_expanded %>%
    filter(
      Family == fam,
      species %in% top10_species
    ) %>%
    group_by(species, Infection) %>%
    summarise(
      ISG_mean_m = mean(ISG_mean, na.rm = TRUE),
      n = n(),
      .groups = "drop"
    ) %>%
    mutate(log_num = log10(n), Family = fam)
  
  data_fam_lmer <- combined_expanded %>%
    filter(
      Family == fam,
      species %in% top10_species
    )
  
  # 線形混合モデル：ISG_score ~ Infection + (1 | species)
  model <- lmer(ISG_mean ~ Infection + (1 | species), data = data_fam_lmer)
  summary_model <- summary(model)
  infection_result <- summary_model$coefficients["InfectionPositive", ]
  
  # モデル結果を表に追加
  result_row <- data.frame(
    Family = fam,
    Estimate = infection_result["Estimate"],
    Std_Error = infection_result["Std. Error"],
    DF = infection_result["df"],
    t_value = infection_result["t value"],
    p_value = infection_result["Pr(>|t|)"]
  )
  lmer_results_df <- rbind(lmer_results_df, result_row)
  
  # アスタリスクラベル作成
  pval <- infection_result["Pr(>|t|)"]
  label <- ifelse(pval < 0.001, "***",
                  ifelse(pval < 0.01, "**",
                         ifelse(pval < 0.05, "*", "")))
  data_fam$p_label <- label
  
  # 統合
  all_plot_data <- bind_rows(all_plot_data, data_fam)
  
  ## 5. プロット
  p <- ggplot(
    data_fam,
    aes(x = Infection, y = ISG_mean_m, group = species, color = species)
  ) +
    geom_line(linewidth = 1) +
    geom_point(aes(size = log_num), alpha = 1.0) +  
    scale_color_viridis_d(option = "D") +
    scale_size_continuous(range = c(1, 5), limits = c(0.4, 4.5), name = "Sample Size") +  
    theme_classic(base_size = 12) +
    labs(
      title = paste0("ISG Score by Infection (Top10 species) - ", fam),
      x = "Infection Status",
      y = "Mean ISG Score"
    ) +
    coord_cartesian(ylim = c(-1.2, 1.2)) +
    annotate("text", x = 1.5, y = 1.2, label = label, size = 6)
  
  plot_list[[fam]] <- p
}

# 4. 表示
plot_list[["Togaviridae"]]
plot_list[["Flaviviridae"]]
plot_list[["Peribunyaviridae"]]
plot_list[["Bornaviridae"]]
plot_list[["Parvoviridae"]]
plot_list[["Smacoviridae"]]
plot_list[["Anelloviridae"]]

for (fam in names(plot_list)) {
  #ggsave(paste0("ISG_lineplot_", fam, "_260110.pdf"), plot_list[[fam]], width = 4, height = 3)
  #ggsave(paste0("ISG_lineplot_", fam, "_260112.pdf"), plot_list[[fam]], width = 4, height = 3)
  #ggsave(paste0("ISG_lineplot_", fam, "_260113.pdf"), plot_list[[fam]], width = 4, height = 3)
  #ggsave(paste0("ISG_lineplot_", fam, "_260205.pdf"), plot_list[[fam]], width = 2, height = 4)
}
print(lmer_results_df)

wrap <- wrap_plots(plot_list, ncol = 1)
wrap
#ggsave("combined_plots.pdf", plot = wrap, width = 5, height = 25)

