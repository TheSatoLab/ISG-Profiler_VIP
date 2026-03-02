rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ggplot2")
library(psych)
library(scales)

all_filt <- read.table("modified_humref_ori_260217.txt", header = T, sep = "\t")
ISG_score <- read.table("modified_ISGcore_humref_260217.txt", header = T, sep = "\t")

#ISG scoreをとって発現量全体の傾向を見るようにする
ref_matrix <- tribble(
  ~group,         ~original, ~human, ~homininae, ~hominidae, ~hominoidea, ~catarrhini, ~simiiformes, ~haplorrhini, ~primates,
  "human",            1,        0,        0,           0,            0,           0,           0,         0,           0,
  "homininae",        1,        1,        0,           0,            0,           0,           0,         0,           0,
  "hominidae",        1,        1,        1,           0,            0,           0,           0,         0,           0,
  "hominoidea",       1,        1,        1,           1,            0,           0,           0,         0,           0,
  "catarrhini",       1,        1,        1,           1,            1,           0,           0,         0,           0,
  "simiiformes",      1,        1,        1,           1,            1,           1,           0,         0,           0,
  "haplorrhini",      1,        1,        1,           1,            1,           1,           1,         0,           0,
  "primates",         1,        1,        1,           1,            1,           1,           1,         1,           0
)
# 並び順を指定
order_levels <- c("original", "human", "homininae", "hominidae", "hominoidea", "catarrhini", "simiiformes", "haplorrhini", "primates")
order_levels_y <- c("human", "homininae", "hominidae", "hominoidea", "catarrhini", "simiiformes", "haplorrhini", "primates")

# ロング形式へ
plot_df <- ref_matrix %>%
  pivot_longer(-group, names_to = "ISG_DB", values_to = "included") %>%
  mutate(
    ISG_DB = factor(ISG_DB, levels = order_levels),
    group = factor(group, levels = rev(order_levels)))

# 描画
status <- ggplot(plot_df, aes(x = ISG_DB, y = group)) +
  geom_point(aes(fill = factor(included)), shape = 21, size = 5, color = "gray50") +
  scale_fill_manual(
    values = c("0" = "gray20", "1" = "gray80"),
    labels = c("Absence", "Presence"),
    name = "Reference"
  ) +
  theme_minimal(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
    #axis.text.x = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(
    #x = "ISG DB version",
    #y = "Phylogenetic group"
    x = NULL,
    y = NULL
  )
status


ISG_score$noRef <- factor(ISG_score$noRef, levels = order_levels)

#original panel の Positive / Negative の中央値を計算
medians <- ISG_score %>%
  filter(noRef == "original") %>%
  group_by(con) %>%
  summarise(med = median(ISG_mean, na.rm = TRUE)) %>%
  deframe()

#プロットに追加
noref_box <- ggplot(ISG_score, aes(x = con, y = ISG_mean, fill = con)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +
  geom_jitter(aes(color = con), shape = 21, color = "black", width = 0.2, alpha = 1.0, size = 1.5) +
  stat_compare_means(aes(group = con), method = "t.test", label = "p.signif") +
  facet_grid(. ~ noRef) +
  theme_classic() +
  theme(axis.text.x = element_text()) +
  labs(title = "ISG score of human samples, Reference Availability in ISG DB Versions",
       y = "Mean ISG Score", x = "", fill = "Induction") +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) +
  scale_x_discrete(
    limits = c("Negative", "Positive"),
    labels = c("Negative" = "-", "Positive" = "+")
  ) +
  geom_hline(yintercept = medians["Positive"], linetype = "dashed", color = "#E41A1C") +
  geom_hline(yintercept = medians["Negative"], linetype = "dashed", color = "#377EB8")

plot(noref_box)


#z-score difference HM
pos_spp_zmed <- ISG_score %>%
  filter(con == "Positive") %>%
  group_by(noRef) %>%
  summarise(mean_pos = mean(ISG_mean, na.rm = TRUE),
            n_pos    = n(), .groups = "drop")

neg_spp_zmed <- ISG_score %>%
  filter(con == "Negative") %>%
  group_by(noRef) %>%
  summarise(mean_neg = mean(ISG_mean, na.rm = TRUE),
            n_neg    = n(), .groups = "drop")

log2fc_species_zmed <- inner_join(pos_spp_zmed, neg_spp_zmed, by = c("noRef")) %>%
  mutate(pos_neg = mean_pos - mean_neg) %>%
  mutate(n_total = n_pos + n_neg) %>%
  mutate(pos_rate = n_pos / n_total)

log2fc_species_zmed$noRef <- factor(log2fc_species_zmed$noRef, levels = order_levels)

heatmap_ref_effect <- ggplot(log2fc_species_zmed,
                             aes(x = noRef, y = "ISG difference", fill = pos_neg)) +
  geom_tile(color = "grey80", height = 0.9) +  
  scale_fill_gradient(
    low = "white",
    high = "green4", 
    limits = c(0, 1.25),
    oob = scales::squish,
    name = "pos - neg"
  ) +
  theme_minimal() +
  theme(
    axis.title.y = element_blank(),
    axis.title.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.text.x = element_blank())
heatmap_ref_effect

final_plot <-  noref_box + heatmap_ref_effect + status +
  plot_layout(heights = c(5, 1, 3), guides = "collect")
plot(final_plot)

output_dir <- "."
hm_tree_ISG_file <- paste0(output_dir, "/", "robust_human_modISGDB_260115.pdf")
#ggsave(hm_tree_ISG_file, plot = final_plot, height=9, width = 8)
