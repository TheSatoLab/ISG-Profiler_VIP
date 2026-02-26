rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ggplot2")
library(ggtree)
library(patchwork)
library(ggnewscale)

taxid_info <- read.table("Amniota398_sp_id.list", header = T, sep = "\t")
ISGcntl_filt_norm.list <- read.table("perspecies_ISG.txt", header = T, sep = "\t")

#get 398 Amniota taxonomy info
library(ape)
library(taxize)

Sys.setenv(ENTREZ_KEY = "hoge")

# tax_id ベクトルの抽出（NA除外＆重複除去）
tax_ids <- taxid_info$tax_id %>% na.omit() %>% unique()

# 一括でNCBIから分類階層を取得
classifications_list <- classification(tax_ids, db = "ncbi")
names(classifications_list) <- as.character(names(classifications_list)) 

# 指定した階級名（例: order）を取り出す関数
extract_rank <- function(class_list, tid, rank) {
  tid <- as.character(tid)
  if (tid %in% names(class_list)) {
    cl <- class_list[[tid]]
    if (!is.null(cl) && inherits(cl, "data.frame")) {
      row <- cl %>% filter(rank == !!rank)
      if (nrow(row) > 0) return(row$name[1])
    }
  }
  return(NA_character_)
}

# 階級情報を taxid_info に付け加える
taxid_info_mod <- taxid_info %>%
  mutate(
    tax_id_chr = as.character(tax_id),
    phylum = map_chr(tax_id_chr, ~ extract_rank(classifications_list, .x, "phylum")),
    class  = map_chr(tax_id_chr, ~ extract_rank(classifications_list, .x, "class")),
    order  = map_chr(tax_id_chr, ~ extract_rank(classifications_list, .x, "order")),
    family = map_chr(tax_id_chr, ~ extract_rank(classifications_list, .x, "family")),
    genus  = map_chr(tax_id_chr, ~ extract_rank(classifications_list, .x, "genus"))
  ) %>%
  dplyr::select(-tax_id_chr)

#write.table(taxid_info_mod, "../control/Amniota398_sp_id_info.list", row.names = F, sep = "\t", quote = F)


# 系統樹のターミナルラベルを tax_id から species 名に変換
tree <- read.tree("../control/tree_Amni_398_mod.nwk")
tree$tip.label <- taxid_info_mod$species[match(tree$tip.label, taxid_info_mod$tax_id)]
print(tree)

# 系統樹をデータフレームに変換
tree_data <- ggtree(tree)$data

# species と order の対応表（重複削除）
species_order_df <- taxid_info_mod %>%
  filter(species %in% tree$tip.label) %>%
  dplyr::select(species, order) %>%
  distinct()

library(purrr)
library(dplyr)

ord_groups <- species_order_df %>%
  group_by(order) %>%
  summarise(
    tips = list(species),
    n_tips = length(tips[[1]]),
    .groups = "drop"
  )

mrca_df <- ord_groups %>%
  mutate(
    node = map2_int(tips, n_tips, ~{
      tips_in_tree <- intersect(.x, tree$tip.label)
      if (length(tips_in_tree) >= 2) {
        ape::getMRCA(tree, tips_in_tree)
      } else if (length(tips_in_tree) == 1) {
        match(tips_in_tree, tree$tip.label)
      } else {
        NA_integer_
      }
    })
  ) %>%
  filter(!is.na(node))

p_tree <- ggtree(tree, layout = "rectangular") +
  coord_flip() +
  theme_tree2() +
  theme(
    axis.title.x = element_blank(),
    axis.text.x = element_blank(),
    axis.ticks.x = element_blank(),
    axis.line.x = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    panel.background = element_blank()
  )

# order ラベルを追加
p_tree <- p_tree +
  geom_cladelab(
    data = mrca_df,
    mapping = aes(node = node, label = order),
    align = TRUE,
    offset = 0.5,
    barsize = 0.3,
    fontsize = 2.6,
    angle = 90,
    inherit.aes = FALSE
  )

p_tree

tip_labels_ordered <- tree_data %>%
  filter(isTip) %>%   
  arrange(y) %>%       
  pull(label)           

ISG_filt_norm.list <- ISGcntl_filt_norm.list %>%
  filter(type == "ISG")

all_combinations <- expand.grid(
  species_db = unique(ISG_filt_norm.list$species_db),
  hum_symbol = unique(ISG_filt_norm.list$hum_symbol),
  stringsAsFactors = FALSE
)

#ループでサンプルごとのheatmapを作成する
row.line <- length(unique(ISGcntl_filt_norm.list$ID))
for (i in 1:row.line) {
  id <- unique(ISGcntl_filt_norm.list$ID)[i]
  # idのデバッグ用出力
  message(paste("Processing SRA:", id))
  
  #サンプル読み込み&系統順に並び替え
  sample.i <- ISGcntl_filt_norm.list %>%
    dplyr::filter(ID == !!id, type == "ISG") %>%
    mutate(species = factor(species, levels = tip_labels_ordered)) %>%
    right_join(all_combinations, by = c("species_db", "hum_symbol"))
  if (!all(levels(sample.i$species) %in% tip_labels_ordered)) {
    warning("ヒートマップのspeciesとtreeのtip.labelに不一致があります。")
  }
  if (nrow(sample.i) == 0) {
    message(paste("No data for SRA:", id, "- skipping."))
    next
  }
  
  # NAのみの場合のチェック
  if (all(is.na(sample.i$ISG_score))) {
    message(paste("Skipping", id, "- all ISG_score values are NA."))
    next
  }
  
  sample.i <- sample.i %>%
    dplyr::mutate(species_db = factor(species_db, levels = tip_labels_ordered))
  
  # 4. ヒートマップの作成
  heatmap_plot <- ggplot(sample.i, aes(x = species_db, y = hum_symbol)) +
    geom_tile(aes(fill = ISG_score)) +
    scale_fill_viridis_c(
      option = "D",
      oob = scales::squish,
      limits = c(-3, 2),
      name = "Expression level",
      na.value = "gray50"  
    ) +
    theme_minimal() +
    theme(
      panel.grid = element_blank(),
      axis.text.x = element_text(angle = 90, hjust = 1, size = 5, vjust = 0.5),
      axis.text.y = element_text(size = 5),
      axis.title.y = element_blank()
    ) +
    labs(
      title = "ISGs Gene Expression Heatmap",
      x = NULL,
      fill = "Expression"
    ) +
    scale_x_discrete(expand = c(0, 0)) +
    scale_y_discrete(expand = c(0, 0))
  
  # 表示
  heatmap_plot
  
  # barplot of reads number for each species
  reads_num_isg <- sample.i %>%
    select(species_db, species_mean_ISG) %>% unique()
  dotplot_isg_plot <- ggplot(reads_num_isg,
                             aes(x = species_db, y = species_mean_ISG)) +
    geom_point(
      color = "gray30",
      size = 1.8,
      alpha = 0.8,
      na.rm = TRUE
    ) +
    theme_minimal() +
    theme(
      axis.line = element_line(color = "black", linewidth = 0.5),
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.ticks.x = element_blank(),
      axis.title.x = element_blank(),
      axis.title.y = element_text(size = 6)
    ) +
    labs(
      y = "ISG score",
      title = paste("ISG Expression for", id)
    )
  
  dotplot_isg_plot
  
  # Adjust layout
  final_plot <- dotplot_isg_plot + heatmap_plot + p_tree +
    plot_layout(heights = c(1, 6, 1))
  
  plot(final_plot)
  output_dir <- "./"
  if (!dir.exists(output_dir)) {
    dir.create(output_dir)
  }
  
  file_name <- file.path(output_dir, paste0(id, "_sp_plots_260120.pdf"))
  #ggsave(file_name, plot = final_plot, width = 25, height = 10)
  
}
