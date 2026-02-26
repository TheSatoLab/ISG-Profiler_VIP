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

#遺伝子の保存状況を調べる
gene2ref <- read.table("gene2refseq_Amniota_ISGcntl_241003.list", header = T, sep = "\t")
tax_comb_org <- read.table("Mammals_Aves_20k/taxon_info_mod_250610.txt", sep = "\t", header = T)
tax_comb <- tax_comb_org %>%
  select(-Host_species) %>%
  unique()

gene2ref_mod <- gene2ref %>%
  select(-Isoform) %>%
  unique()

gene2ref_tax <- gene2ref_mod %>%
  inner_join(tax_comb, by = "tax_id")

gene2ref_tax_gene <- gene2ref_tax %>%
  select(tax_id, species, order, type, hum_symbol) %>%
  group_by(order, type, hum_symbol) %>%
  summarise(n_sp = n())

gene2ref_tax_count <- gene2ref_tax %>%
  group_by(tax_id, species, order, type) %>%
  summarise(n_genes = n())

gene2ref_tax_mod <- gene2ref_tax_count %>%
  group_by(order, type) %>%
  summarise(n_genes_mean = mean(n_genes))

#write.table(gene2ref_tax_mod, "gene_n_250421.txt", row.names = F, sep = "\t", quote = F)

hum_symbol_type <- gene2ref_tax %>%
  select(hum_symbol, type) %>%
  distinct()

# 全 order × hum_symbol の組み合わせ
full_grid <- expand_grid(
  order = unique(gene2ref_tax$order),
  hum_symbol = unique(gene2ref_tax$hum_symbol)
) %>%
  left_join(hum_symbol_type, by = "hum_symbol")  

# 1. orderごとの生物種数（分母）
order_gene_n <- gene2ref_tax %>%
  distinct(order, species, hum_symbol, type) %>%
  count(order, hum_symbol, type, name = "n_species_with_gene")

# 組み合わせにマージしてNA→0
order_gene_coverage <- full_grid %>%
  left_join(order_gene_n, by = c("order", "hum_symbol", "type")) %>%
  mutate(n_species_with_gene = replace_na(n_species_with_gene, 0))

# 分母（orderごとの生物種数）
order_species_n <- gene2ref_tax %>%
  select(species, order) %>%
  distinct() %>%
  count(order, name = "n_species")

# 保存率計算
order_gene_coverage <- order_gene_coverage %>%
  left_join(order_species_n, by = "order") %>%
  mutate(coverage = n_species_with_gene / n_species * 100)

low_cov <- order_gene_coverage %>%
  filter(coverage<50)

order_gene_coverage_ISG <- order_gene_coverage %>%
  filter(type == "ISG")

#tree
taxid_info <- tax_comb_org %>%
  select(Host_species, tax_id)
tree <- read.tree("../virome/geNomad/Mammals_Aves_20k/timetree_250910.nwk")
target_orders <- unique(order_gene_coverage_ISG$order)

# species→orderの対応表（treeに存在するspeciesだけ）
sp_order <- tax_comb %>%
  filter(order %in% target_orders,
         species %in% tree$tip.label) %>%
  distinct(species, order) %>%
  group_by(order) %>%
  slice(1) %>%                     
  ungroup()

stable_species_order <- readLines("stable_species_order.txt")
# stable順に species をフィルタし、それに対応する order の順番を取得
ordered_orders <- sp_order %>%
  filter(species %in% stable_species_order) %>%
  mutate(species = factor(species, levels = stable_species_order)) %>%
  arrange(species) %>%
  pull(order) %>%
  unique()  

# treeにいる全部のspecies
all_tips <- tree$tip.label

# 落としたいtip = 代表に選ばれていないspecies
tips_to_drop <- setdiff(all_tips, sp_order$species)

# 枝刈り
tree_order <- drop.tip(tree, tips_to_drop)
# species → order に変換する lookup
match_idx <- match(tree_order$tip.label, sp_order$species)
tree_order$tip.label <- sp_order$order[match_idx]

# `fortify(tree)` を使ってデータフレーム化
tree_data <- fortify(tree_order)

# 系統樹をデータフレームに変換
tree_data <- ggtree(tree_order)$data

# ターミナルノード（葉）を抽出し、樹形順にラベルを取得
tip_labels_in_order <- tree_data %>%
  filter(isTip) %>%  
  arrange(y) %>%     
  pull(label)        

# 系統樹の描画
tree_plot <- ggtree(tree_data, layout = "rectangular") +
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

order_gene_coverage_ISG <- order_gene_coverage_ISG %>%
  mutate(order = factor(order_gene_coverage_ISG$order, levels = tip_labels_in_order))


# ① wide 形式にして、クラスタリング対象の行列を作成
heatmap_mat <- order_gene_coverage_ISG %>%
  select(order, hum_symbol, coverage) %>%
  pivot_wider(names_from = order, values_from = coverage, values_fill = 0) %>%
  column_to_rownames("hum_symbol")

# ② クラスタリング
dist_mat <- dist(heatmap_mat)
clust <- hclust(dist_mat)
clust_ward <- hclust(dist_mat, method = "ward.D2")
clust_avg  <- hclust(dist_mat, method = "average")
clust_mc   <- hclust(dist_mat, method = "mcquitty")
ordered_symbols <- clust$labels[clust$order] 
ordered_symbolsw <- clust$labels[clust_ward$order]  
ordered_symbolsa <- clust$labels[clust_avg$order]  
ordered_symbolsm <- clust$labels[clust_mc$order]  

# ③ hum_symbolの因子順をクラスタ順に並べる
order_gene_coverage_ISG <- order_gene_coverage_ISG %>%
  mutate(hum_symbol = factor(hum_symbol, levels = ordered_symbolsw))

# heatmapデータのorderカラムのfactor順を明示的に指定
order_gene_coverage_ISG <- order_gene_coverage_ISG %>%
  mutate(order = factor(order, levels = ordered_orders))

# 指定した順番
prioritized <- c("ZC3HAV1", "OAS1", "TAP2", "MX1", "TRIM21", "IFIT3", "IFIT2")

# 残りのクラスタ順
remaining <- ordered_symbolsw[!ordered_symbolsw %in% prioritized]

# 優先順を最後に（上に）したいので、後ろに足す
final_order <- c(remaining, prioritized)

# 因子順を更新（上に来てほしい遺伝子がY軸上部に来る）
order_gene_coverage_ISG <- order_gene_coverage_ISG %>%
  mutate(hum_symbol = factor(hum_symbol, levels = final_order))

# ④ ggplotで描画（先ほどのコードと同じ）
coverage_HM_ISG <- ggplot(order_gene_coverage_ISG, aes(x = order, y = hum_symbol, fill = coverage)) +
  geom_tile(color = "white") +
  scale_fill_gradientn(
    colours = c("white", "skyblue", "blue4"),
    values = scales::rescale(c(0, 0.01, 50, 100)),
    limits = c(0, 100),
    # na.value = "grey80",
    guide = guide_colorbar(title = "Coverage (%)")
  ) +
  facet_grid(rows = vars(type), scales = "free_y", space = "free_y") +
  scale_x_discrete(drop = FALSE) +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1),
    axis.title.x = element_blank(),
    axis.title.y = element_blank(),
    panel.grid = element_blank()
  ) +
  labs(title = "Gene Conservation Rate by Order")

plot(coverage_HM_ISG)

output_dir <- "."
hm_tree_ISG_file <- paste0(output_dir, "/", "coreISG_conserve_260126.pdf")
#ggsave(hm_tree_ISG_file, plot = coverage_HM_ISG, height=4, width = 4)
