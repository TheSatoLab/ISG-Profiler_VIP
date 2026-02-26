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
library(scales)

#import data
ISG_meta <- read.table("Aves_Mam_mbio_250613_ISGscore_logan.txt", sep = "\t", header = T)
Aves_Mam_mbio_ISG_norm.list <- read.table("Aves_Mam_mbio_250613.ISGcntl_norm_logan.txt", sep = "\t", header = T)
Aves_Mam_mbio_genomad_filt <- read.table("Aves_Mam_mbio_genomad_filt_ML_250630.txt", sep = "\t", header = T)
Aves_Mam_mbio_meta_filt <- read.table("Aves_Mam_mbio_metadata_ISG_logan_ML_250616_2.txt", sep = "\t", header = T)

tmp_gen <- Aves_Mam_mbio_genomad_filt %>%
  select(-n_contig, -ISG_level)
ISG_virus <- ISG_meta %>%
  left_join(tmp_gen, by = 'ID', relationship = "many-to-many")

ISG_virus_mod <- ISG_virus %>%
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

ISG_virus_meta <- ISG_virus_mod %>%
  group_by(Infection, Family, order, class, genus, species) %>%
  summarise(
    ISG_mean_m = mean(ISG_mean, na.rm = TRUE),  
    sample_size = n(),  
    .groups = "drop"
  ) %>%
  ungroup()


heatmap_species <- unique(ISG_virus$species)

ordered_families <- c("DNA Unclassified", "Herpesviridae", "Papillomaviridae", "Polyomaviridae",
                      "Parvoviridae", "Circoviridae", "Smacoviridae",
                      "Anelloviridae", "Poxviridae", "Adintoviridae",
                      "Adenoviridae", "Hepadnaviridae", "Birnaviridae",
                      "Sedoreoviridae", "Spinareoviridae", "Coronaviridae",
                      "Arteriviridae", "Tobaniviridae", "Picornaviridae",
                      "Caliciviridae", "Astroviridae", "Flaviviridae",
                      "Hepeviridae", "Togaviridae", "Nodaviridae",
                      "Phenuiviridae", "Peribunyaviridae", "Arenaviridae",
                      "Orthomyxoviridae", "Paramyxoviridae", "Bornaviridae",
                      "Rhabdoviridae", "Pneumoviridae", "Filoviridae",
                      "RNA Unclassified", "Negative")

ISG_virus_meta <- ISG_virus_meta %>%
  mutate(Family = factor(Family, levels = ordered_families))

# 1. データのフィルタリング。先に平均をとった値でASR
neg_data <- ISG_virus_meta %>%
  filter(Family == "Negative") %>%
  filter(!is.na(ISG_mean_m)) %>% 
  filter(!is.na(species)) %>%
  filter(sample_size >= 5) %>%
  mutate(ISG_scale = as.numeric(scale(ISG_mean_m)))

tmp <- neg_data %>%
  select(species) %>%
  unique()

tmp2 <- neg_data %>%
  select(order, class, genus, species, ISG_mean_m) %>%
  unique()

virus_num <- Aves_Mam_mbio_meta_filt %>%
  inner_join(tmp, by = "species") %>%
  group_by(species, Infection) %>%
  summarise(inf_num = n()) %>%
  ungroup() %>%
  group_by(species) %>%
  summarise(
    total_num = sum(inf_num, na.rm = TRUE),
    pos_num   = sum(inf_num[Infection == "Positive"], na.rm = TRUE),
    neg_num   = sum(inf_num[Infection == "Negative"], na.rm = TRUE),
    infection_ratio = pos_num / total_num
  ) %>%
  ungroup() %>%
  inner_join(tmp2, by = "species")

order_rm <- virus_num %>%
  group_by(order) %>%
  summarise(n = n()) %>%
  filter(n>=3)

virus_num_filt <- virus_num %>%
  filter(order %in% order_rm$order)

tree <- read.tree("timetree_250910.nwk")
# 使いたいorder（virus_num_filtに出てくるorder）
target_orders <- unique(virus_num_filt$order)

# species→orderの対応表（treeに存在するspeciesだけ）
sp_order <- virus_num_filt %>%
  filter(order %in% target_orders,
         species %in% tree$tip.label) %>%
  distinct(species, order) %>%
  group_by(order) %>%
  slice(1) %>%                    
  ungroup()

# treeにいる全部のspecies
all_tips <- tree$tip.label

# 落としたいtip = 代表に選ばれていないspecies
tips_to_drop <- setdiff(all_tips, sp_order$species)

# 枝刈り
tree_order <- drop.tip(tree, tips_to_drop)
# species → order に変換する lookup
match_idx <- match(tree_order$tip.label, sp_order$species)
tree_order$tip.label <- sp_order$order[match_idx]

library(ggtree)
tree_plot <- ggtree(tree_order, layout = "rectangular") +
  coord_flip() +
  theme_tree2() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),       
    panel.grid.major = element_blank(),  
    panel.grid.minor = element_blank(),  
    panel.background = element_blank(), 
  )
plot(tree_plot)

tree_data <- ggtree(tree_order)$data

# ターミナルノード（葉）を抽出し、樹形順にラベルを取得
tip_labels_in_order <- tree_data %>%
  filter(isTip) %>% 
  arrange(y) %>%     
  pull(label)        

virus_num_filt <- virus_num_filt %>%
  mutate(order = factor(order, levels = tip_labels_in_order))

Negative_host <- ggplot(virus_num_filt, aes(x=order, y=ISG_mean_m, fill = class)) +
  geom_boxplot() +
  geom_jitter(width = 0.2, size = 2, alpha = 0.7) +
  #facet_wrap(~ BioProject_ID) +
  theme_classic() +
  theme(axis.text.x = element_blank(),
        axis.title.x = element_blank()) +
  labs(title = "Basal ISG expression level on each animal order",
       x = "", y = "ISG_score")

plot(Negative_host)

#STAT1とIRF7のHMを結合させる
stat_irf <- data.frame(hum_symbol = c("STAT1", "IRF7"))
stat_irf_norm <- Aves_Mam_mbio_ISG_norm.list %>%
  inner_join(stat_irf, by = "hum_symbol")

ISG_virus_mod_neg <- ISG_virus_mod %>%
  filter(Family == "Negative") %>%
  inner_join(tmp, by = "species") %>%
  select(ID) %>%
  unique()

stat_irf_norm_neg <- stat_irf_norm %>%
  inner_join(ISG_virus_mod_neg, by = "ID")

stat_irf_norm_neg_sp <- stat_irf_norm_neg %>%
  ungroup() %>%
  group_by(species, hum_symbol, order, class) %>%
  summarise(ISG_mean_sp = mean(ISG_score)) %>%
  ungroup()

stat_irf_norm_neg_ord <- stat_irf_norm_neg_sp %>%
  group_by(hum_symbol, order, class) %>%
  summarise(ISG_mean_ord = mean(ISG_mean_sp), n = n()) %>%
  ungroup() %>%
  filter(n>=3)

stat_irf_norm_neg_ord <- stat_irf_norm_neg_ord %>%
  mutate(order = factor(order, levels = tip_labels_in_order))

p_heatmap <- ggplot(stat_irf_norm_neg_ord, aes(x = order, y = hum_symbol, fill = ISG_mean_ord)) +
  geom_tile(color = "white") +
  scale_fill_gradient(
    low = "white",       
    high = "darkgreen",      
    limits = c(-1.0, 0.3),   
    oob = scales::squish,
    name = "Expression level"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = NULL)

p_heatmap <- ggplot(stat_irf_norm_neg_ord, aes(x = order, y = hum_symbol, fill = ISG_mean_ord)) +
  geom_tile(color = "white") +
  scale_fill_viridis_c(
    option = "D",           
    limits = c(-0.5, 0.1),  
    oob = scales::squish,
    name = "Expression level"
  ) +
  theme_minimal(base_size = 11) +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title = element_blank(),
    panel.grid = element_blank()
  )
print(p_heatmap)


tree_box_p <- Negative_host + p_heatmap + tree_plot +
  plot_layout(heights = c(8, 1, 1), guides = "collect")
plot(tree_box_p)

output_dir <- "./"
tree_box_file_name <- paste0(output_dir, "/", "basalISG_order_ISGHM_260110.pdf")
#ggsave(tree_box_file_name, plot = tree_box_p, width = 6, height = 6)

isg_basal <- neg_data %>%
  filter(Infection == "Negative") %>%
  filter(!is.na(order), !is.na(ISG_mean_m))  

isg_basal <- isg_basal %>%
  mutate(order = as.factor(order))

lm_order <- lm(ISG_mean_m ~ order, data = isg_basal)
summary(lm_order)
anova(lm_order)

library(multcomp)
tukey_res <- glht(lm_order, linfct = mcp(order = "Tukey"))
summary(tukey_res)

library(ggpubr)
ggboxplot(isg_basal, x = "order", y = "ISG_mean_m",
          color = "order", palette = "jco", add = "jitter") +
  stat_compare_means(method = "anova")

filtered_data <- ISG_virus_mod %>%
  filter(Infection == "Negative", species %in% neg_data$species)

filtered_data$order <- as.factor(filtered_data$order)

lm_order <- lm(ISG_mean ~ order, data = filtered_data)
summary(lm_order)

# 多重比較
library(multcomp)
tukey_res <- glht(lm_order, linfct = mcp(order = "Tukey"))
summary(tukey_res, test = adjusted("none"))

# Orderごとに平均と標準誤差を計算
order_summary <- filtered_data %>%
  group_by(order) %>%
  summarise(
    mean_ISG = mean(ISG_mean, na.rm = TRUE),
    se_ISG = sd(ISG_mean, na.rm = TRUE) / sqrt(n()),
    n = n()
  )

ggplot(order_summary, aes(x = reorder(order, mean_ISG), y = mean_ISG)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  geom_errorbar(aes(ymin = mean_ISG - se_ISG, ymax = mean_ISG + se_ISG), width = 0.3) +
  coord_flip() +
  theme_classic(base_size = 14) +
  labs(
    title = "Basal ISG Score by Animal Order",
    x = "Order",
    y = "Mean ISG Score (± SE)"
  )

plot(tukey_res, las = 1)  # 横軸が Order の組み合わせ、縦軸が差

ggplot(filtered_data, aes(x = order, y = ISG_mean, fill = order)) +
  geom_boxplot(outlier.shape = NA) +
  geom_jitter(width = 0.2, alpha = 0.3, size = 0.7) +
  theme_classic(base_size = 13) +
  coord_flip() +
  labs(
    title = "Distribution of Basal ISG Scores by Order",
    x = "Order",
    y = "ISG Score"
  ) +
  theme(legend.position = "none")


# モデルの係数（estimate, conf.intなど）を取り出す
coef_df <- broom::tidy(lm_order, conf.int = TRUE)

tmp3 <- virus_num %>%
  dplyr::select(order, class) %>%
  unique()

# Intercept 以外を抽出
coef_df_filtered <- coef_df %>%
  filter(term != "(Intercept)") %>%
  mutate(term = str_replace(term, "order", "")) %>%
  rename(order = term) %>%
  inner_join(tmp3, by = "order") %>%
  mutate(signif_label = ifelse(p.value < 0.05, "*", ""))

# プロット
lm_order_p <- ggplot(coef_df_filtered, aes(x = estimate, y = reorder(order, estimate), color = class)) +
    geom_point(size = 3) +
    geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2) +
  　geom_text(
    aes(x = estimate + 0.4, label = signif_label),  
    color = "black",
    size = 5,
    show.legend = FALSE
  　) +
    coord_flip() +
    theme_classic(base_size = 12) +
    labs(
      x = "Effect size on ISG score",
      title = "Effect of Taxonomic Order on ISG Score"
    ) +
    theme(
      axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
      axis.title.x = element_blank(),
      panel.grid = element_blank()
    )
lm_order_p

output_dir <- "./"
lm_order_name <- paste0(output_dir, "/", "basalISG_order_lm_260113.pdf")
#ggsave(lm_order_name, plot = lm_order_p, width = 8, height = 5)
