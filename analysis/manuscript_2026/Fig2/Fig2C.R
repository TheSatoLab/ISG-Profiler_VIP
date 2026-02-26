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
library(castor)
library(cowplot)

setwd("~/OneDrive/ドキュメント/解析結果/Virome/virome/geNomad/Mammals_Aves_20k/")

#import data
ISG_meta <- read.table("Aves_Mam_mbio_250613_ISGscore_logan.txt", sep = "\t", header = T)
Aves_Mam_mbio_meta_filt <- read.table("Aves_Mam_mbio_metadata_ISG_logan_ML_250616_2.txt", sep = "\t", header = T)
Aves_Mam_mbio_genomad_filt <- read.table("Aves_Mam_mbio_genomad_filt_ML_250630.txt", sep = "\t", header = T)

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
    ISG_median = median(ISG_mean, na.rm = TRUE),  # 🔹 ISGスコアの中央値
    sample_size = n(),  # 🔹 サンプルサイズ
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

#ASR解析
# 1. データのフィルタリング。先に平均をとった値でASR
neg_data <- ISG_virus_meta %>%
  filter(Family == "Negative") %>%
  filter(!is.na(ISG_median)) %>%
  filter(!is.na(species)) %>%
  filter(sample_size >= 5) %>%
  mutate(ISG_scale = as.numeric(scale(ISG_median)))

tmp <- neg_data %>%
  select(species) %>%
  unique()

tmp2 <- neg_data %>%
  select(order, class, genus, species, ISG_median) %>%
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

# スコアをベクトル化（species名を名前に）
score_vec <- deframe(neg_data %>% distinct(species, ISG_scale))

tree <- read.tree("timetree_250910.nwk")
# pruning（スコアがある種だけに絞る）
matched_species <- intersect(tree$tip.label, names(score_vec))

pruned_tree2 <- drop.tip(tree, setdiff(tree$tip.label, matched_species))
score_vec2 <- score_vec[pruned_tree2$tip.label]
names(score_vec2) <- pruned_tree2$tip.label

# ASR実行
asr <- asr_squared_change_parsimony(
  tree = pruned_tree2,
  tip_states = score_vec2
)

#tip/nodeのstateをggtreeのデータに合体
td <- ggtree(pruned_tree2)$data
# tip
td$state <- NA_real_

#tips
tip_idx    <- which(td$isTip)
tip_labels <- td$label[tip_idx]
tip_vals   <- score_vec2[tip_labels]

# 参照できなかったラベル（名前不一致）チェック
if (anyNA(tip_vals)) {
  cat("Unmatched tip labels:\n")
  print(tip_labels[is.na(tip_vals)])
}
td$state[tip_idx] <- tip_vals

#node
node_idx <- which(!td$isTip)

# castor の出力に名前を付ける（必須）
node_ids <- (Ntip(pruned_tree2) + 1):(Ntip(pruned_tree2) + pruned_tree2$Nnode)
names(asr$ancestral_states) <- as.character(node_ids)

node_keys <- as.character(td$node[node_idx])   # ggtreeの node 番号（文字列に）
node_vals <- asr$ancestral_states[node_keys]

if (anyNA(node_vals)) {
  cat("Unmatched internal nodes (unexpected):\n")
  print(node_keys[is.na(node_vals)])
}
td$state[node_idx] <- node_vals

# 0 白の発色、外れ値はクリップ
min_val <- -2.0; max_val <- 2.0
td$state_clip <- pmin(pmax(td$state, min_val), max_val)

p_base <- ggtree(pruned_tree2, layout = "rectangular")

#ASR plot
p_asr <- p_base +
  geom_tree() +
  geom_point(
    data = td[!td$isTip, ],
    aes(x = x, y = y, color = state_clip),
    size = 3,
    inherit.aes = FALSE
  ) +
  geom_point(
    data = td[td$isTip, ],
    aes(x = x, y = y, color = state_clip),
    size = 3,
    inherit.aes = FALSE
  ) +
  scale_color_gradient2(low = "blue", mid = "white", high = "red",
                        midpoint = 0, name = "ASR (z)") +
  coord_flip() +
  theme_tree2() +
  theme(
    text = element_text(size = 10),
    axis.line.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.title.x = element_blank(),
    panel.grid = element_blank()
  )


## 1) orderごとのMRCAノードを作る（ツリーに実在するtipのみ）
ord_groups <- neg_data %>%
  filter(!is.na(order)) %>%
  filter(species %in% pruned_tree2$tip.label) %>%
  distinct(species, order) %>%
  group_by(order) %>%
  summarise(tips = list(unique(species)),
            n_tips = length(tips[[1]]),
            .groups = "drop")

mrca_df <- ord_groups %>%
  mutate(
    node = map2_int(tips, n_tips, ~{
      tv <- intersect(.x, pruned_tree2$tip.label)
      if (length(tv) >= 2) {
        res <- suppressWarnings(ape::getMRCA(pruned_tree2, tv))
        if (is.null(res) || length(res) == 0) NA_integer_ else as.integer(res)
      } else if (length(tv) == 1) {
        as.integer(match(tv, pruned_tree2$tip.label))
      } else {
        NA_integer_
      }
    })
  ) %>%
  filter(!is.na(node))

# 既存のp_asrに「orderラベル」だけを追加
# mrca_df: 列に node（整数ノードID）, order（ラベル文字列）がある前提
p_asr <- p_asr +
  ggtree::geom_cladelab(
    data = mrca_df,
    mapping = aes(node = node, label = order),
    align = TRUE,
    offset = 0.5,
    barsize = 0.3,
    fontsize = 2.6,
    angle = 90,                 # 縦向きツリーに合わせて回転
    inherit.aes = FALSE
  )

## 1) ツリーのtip順（描画順）を取得
tip_order <- p_asr$data %>%
  filter(isTip) %>%
  arrange(y) %>%
  pull(label)

#HMにする
# 2) 指標を作って並びを固定
virus_num_sc <- virus_num %>%
  filter(species %in% tip_order) %>%
  mutate(
    species = factor(species, levels = tip_order),
    total_num_log10 = log10(pmax(total_num, 1)),                # Nはlog10
    ISG_median_clip = scales::squish(ISG_median, c(-1.5, 1.5))  # 表示範囲をクリップ
  )

# ASRのggtreeオブジェクトから tip の順番を取る
tip_order <- p_asr$data %>%
  dplyr::filter(isTip) %>%
  dplyr::arrange(y) %>%
  dplyr::pull(label)

# HMとdotplotに同じ順序を強制
virus_num_sc <- virus_num_sc %>%
  dplyr::mutate(species = factor(species, levels = tip_order))

p_hm <- ggplot(virus_num_sc, aes(x = species)) +
  # line 1: log10(N)
  geom_tile(aes(y = "sample number, log10", fill = total_num_log10), height = 1) +
  scale_fill_viridis_c(
    name = "log10(N)",
    guide = guide_colorbar(
      barheight = unit(3, "mm"),
      barwidth  = unit(18, "mm"),
      title.position = "top",
      ticks = FALSE
    )
  ) +
  new_scale_fill() +

  # line 2: 陽性割合（mid=0.5を白）
  geom_tile(aes(y = "Virus positive ratio", fill = infection_ratio), height = 1) +
  scale_fill_gradient(
    low = "white",
    high = "tomato",
    limits = c(0, 1),
    name = "Positive ratio",
    guide = guide_colorbar(
      barheight = unit(3, "mm"),
      barwidth  = unit(18, "mm"),
      title.position = "top",
      ticks = FALSE
    )
  ) +
  new_scale_fill() +

  coord_cartesian(expand = FALSE) +
  scale_x_discrete(drop = FALSE) +
  labs(x = NULL, y = NULL) +
  theme_minimal(base_size = 10) +
  theme(
    axis.text.x  = element_blank(), axis.ticks.x = element_blank(),
    panel.grid   = element_blank(),
    plot.margin  = margin(2, 6, 2, 6),
    legend.position   = "top",
    legend.box        = "horizontal",
    legend.title      = element_text(size = 8),
    legend.text       = element_text(size = 7),
    legend.key.width  = unit(6, "mm"),
    legend.key.height = unit(3, "mm"),
    legend.margin     = margin(0,0,0,0)
  )

p_isg_dot <- ggplot(virus_num_sc, aes(x = species, y = ISG_median)) +
  geom_hline(yintercept = 0, linetype = 2, linewidth = 0.3) +
  geom_point(size = 1.8) +  # 必要なら aes(color = infection_ratio, size = total_num_log10) など足してもOK
  labs(x = NULL, y = "ISG median") +
  theme_classic(base_size = 10) +
  theme(
    axis.text.x  = element_blank(),
    axis.ticks.x = element_blank(),
    panel.grid.minor = element_blank(),
    plot.margin  = margin(2, 6, 2, 6)
  )
p_isg_dot

# 4 上下にスタック（横幅そろえて縦方向で整列）
p_asr <- p_asr + coord_flip(clip = "off")

p_asr2 <- p_asr +
  scale_x_continuous(expand = expansion(mult = c(0.02, 0.1))) +
  theme(plot.margin = margin(5.5, 50, 5.5, 5.5))
fina_p <- plot_grid(p_hm, p_asr2, ncol = 1, align = "v", rel_heights = c(1, 5))
fina_p
#ggsave("asr_tree_HM_251112.png", fina_p, width = 40, height = 10, units = "in", dpi = 300)

p_all <- plot_grid(
  p_hm,
  p_isg_dot,
  p_asr2,
  ncol = 1,
  align = "v",
  axis  = "lr",
  rel_heights = c(1, 1, 3)
)

p_all
#ggsave("asr_tree_HM_251225.pdf", p_all, width = 40, height = 5, units = "in", dpi = 300)

