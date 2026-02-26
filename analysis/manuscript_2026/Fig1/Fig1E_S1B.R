rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ggplot2")
library("dplyr")
library(pryr)
library("caret")
library(ggtree)
library(ape)
library(ggtext)
library(ggpubr)
library(patchwork)

Game_fur_filt <- read.table("Game_fur_ISG_260217.txt",  sep = "\t", header = T)
ISG_score <- read.table("Game_fur_ISGscore_nometa_260217.txt",  sep = "\t", header = T)
virus_ISG_status <- read.table("ISG_virus_family_251111.txt",  sep = "\t", header = T) %>%
  separate(col = taxonomy, into=c("Virus_v", "Realm_v", "Kingdom_v", "Phylum_v", "Class_v", "Order_v", "Family_v"), sep = ";") %>%
  mutate(across(everything(), ~na_if(.x, ""))) %>%
  ungroup()
ISG_meta_update_level <- read.table("Game_Fur_ISG_meta_update.txt",  sep = "\t", header = T)

#get viral taxonomy info
library(ape)
library(taxize)

Sys.setenv(ENTREZ_KEY = "hoge")

# #get additional tax info
# test <- species_df_mod_join %>% filter(is.na(tax_id_v)) %>% select(Viral_species) %>% unique()
# # NCBIからtax_idを取得
# uids <- get_uid(test$Viral_species, names = TRUE)
# uids_named <- setNames(uids, test$Viral_species)
# 
# # 分類情報を取得
# classifications <- classification(uids, db = "ncbi")
# 
# # familyを抽出してDataFrameにする
# species_df <- tibble(
#   species = names(uids_named),
#   tax_id  = as.character(uids_named)
# )
# #write.table(species_df, "Game_Fur_virus_additional_taxid.txt", row.names = F, sep = "\t", quote = F)
# #modified the above file
# 
# # 種名をとる
# species_names <- unique(Game_fur_virus$Viral_species)
# species_df <- read.table("Game_Fur_virus_taxid_mod_260205.txt", sep = "\t", header = T)
# tax_ids <- species_df$tax_id %>% na.omit()
# 
# # 一括で分類情報を取得
# classifications2 <- classification(tax_ids, db = "ncbi")
# names(classifications2) <- as.character(names(classifications2))
# 
# species_df_mod <- species_df %>%
#   filter(!is.na(tax_id)) %>%
#   mutate(
#     family = map_chr(tax_id, function(tid) {
#       tid <- as.character(tid)  
#       if (tid %in% names(classifications2)) {
#         cl <- classifications2[[tid]]
#         if (!is.null(cl) && inherits(cl, "data.frame")) {
#           fam_row <- cl %>% filter(rank == "family")
#           if (nrow(fam_row) > 0) fam_row$name[1] else NA_character_
#         } else {
#           NA_character_
#         }
#       } else {
#         NA_character_
#       }
#     })
#   )
# 
# # 必要な分類階級
# target_ranks <- c("phylum", "class", "order", "family", "genus")
# 
# extract_rank <- function(class_list, tid, rank) {
#   tid <- as.character(tid)
#   if (tid %in% names(class_list)) {
#     cl <- class_list[[tid]]
#     if (!is.null(cl) && inherits(cl, "data.frame")) {
#       row <- cl %>% filter(rank == !!rank)
#       if (nrow(row) > 0) return(row$name[1])
#     }
#   }
#   return(NA_character_)
# }
# 
# # すべての分類階級を抽出するコード
# species_df_mod <- species_df %>%
#   filter(!is.na(tax_id)) %>%
#   mutate(
#     tax_id = as.character(tax_id),  
#     across(.cols = everything(), .fns = identity), 
#     phylum = map_chr(tax_id, ~ extract_rank(classifications2, .x, "phylum")),
#     class = map_chr(tax_id, ~ extract_rank(classifications2, .x, "class")),
#     order = map_chr(tax_id, ~ extract_rank(classifications2, .x, "order")),
#     family = map_chr(tax_id, ~ extract_rank(classifications2, .x, "family")),
#     genus = map_chr(tax_id, ~ extract_rank(classifications2, .x, "genus"))
#   ) %>%
#   unique()
# 
# species_df_mod_join <- species_df_mod %>%
#   rename(Viral_species = species, Phylum_v = phylum, Class_v = class, Order_v = order, Family_v = family, Genus_v = genus, tax_id_v = tax_id) %>%
#   full_join(Game_fur_virus, by = "Viral_species") %>%
#   unique()
# 
# #write.table(species_df_mod_join, "Game_Fur_meta_virus.txt", row.names = F, sep = "\t", quote = F)

species_df_mod_join <- read.table("Game_Fur_meta_virus.txt", header = T, sep = "\t")

#z-score difference, lung
pos_spp_zmed_lung_inf <- ISG_meta_update_level %>%
  filter(tissue == "Lung") %>%
  filter(Infection == "Positive") %>%
  group_by(species, Family_v) %>%
  summarise(mean_pos = mean(ISG_mean, na.rm = TRUE),
            n_pos    = n(), .groups = "drop") %>%
  filter(n_pos >= 2)

neg_spp_zmed_lung_inf <- ISG_meta_update_level %>%
  filter(tissue == "Lung") %>%
  filter(Infection == "Negative") %>%
  group_by(species) %>%
  summarise(mean_neg = mean(ISG_mean, na.rm = TRUE),
            n_neg    = n(), .groups = "drop") %>%
  filter(n_neg >= 2)

log2fc_species_zmed_lung_inf <- inner_join(pos_spp_zmed_lung_inf, neg_spp_zmed_lung_inf, by = c("species")) %>%
  filter(n_pos >= 2, n_neg >= 2) %>%
  mutate(pos_neg = mean_pos - mean_neg)  %>%
  mutate(n_total = n_pos + n_neg) %>%
  mutate(pos_rate = n_pos / n_total)

full_mat_z_lung_inf <- log2fc_species_zmed_lung_inf %>%
  select(species, Family_v, pos_neg) %>%
  pivot_wider(
    names_from = Family_v,
    values_from = pos_neg
  )
hm_df_z_lung_inf <- full_mat_z_lung_inf %>%
  pivot_longer(
    cols = -species,
    names_to = "Family_v",
    values_to = "pos_neg"
  ) %>%
  left_join(log2fc_species_zmed_lung_inf, by = c("species", "Family_v", "pos_neg"))

tip_order <- c("Vulpes_lagopus", "Oryctolagus_cuniculus", "Marmota_himalayana", "Cavia_porcellus")
hm_df_z_lung_inf <- hm_df_z_lung_inf %>%
  mutate(species = factor(species, levels = tip_order)) %>%
  mutate(Family_b = fct_reorder(Family_v, pos_neg, .fun = median, .desc = TRUE))


hm_log2fc_z_lung_inf <- ggplot(hm_df_z_lung_inf, aes(x = species, y = Family_v, fill = pos_neg)) +
  geom_point(aes(size = n_pos, fill = pos_neg), shape = 21, color = "grey0", alpha = 0.8) +
  #geom_tile(color = "grey90") +
  scale_fill_distiller(
    palette = "RdBu",         
    direction = -1,           
    limits = c(-0.5, 0.5),   
    oob = scales::squish,
    name = "pos-neg"
  ) +
  scale_size(range = c(3, 10)) +
  labs(x = "", y = "Virus family", title = "Positive − Negative ISG score, gaming/fur animal, lung") +
  theme_minimal(base_size = 10) +
  theme(
    axis.line = element_line(color = "black"),
    axis.text.x = element_blank(),
    panel.grid.major = element_line(color = "grey85", linetype = "dotted", linewidth = 0.2),  
    panel.grid.minor = element_blank()
  )

hm_log2fc_z_lung_inf

output_dir <- "."
hm_log2fc_z_lung_inf_name <- paste0(output_dir, "/", "gamefur_posneg_Infection_lung_HM_260205.pdf")
#ggsave(hm_log2fc_z_lung_inf_name, plot = hm_log2fc_z_lung_inf, height=5, width = 4)

#t-test
df_inf <- hm_df_z_lung_inf %>%
  filter(!is.na(pos_neg)) %>%
  mutate(x = "All")
test_res_inf <- t.test(df_inf$pos_neg, mu = 0, alternative = "greater")

gf_one_l_inf <- ggplot(df_inf, aes(x = "", y = pos_neg)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +  
  geom_jitter(aes(size = n_pos, color = species), width = 0.1, alpha = 1) + 
  geom_boxplot(fill = "white", outlier.shape = NA, alpha = 0.4) +
  annotate("text",
           x = 1,
           y = max(df_inf$pos_neg, na.rm = TRUE) + 0.1,
           label = paste0("p = ", signif(test_res_inf$p.value, 3)),
           size = 4) +
  scale_size_continuous(name = "Sample count", range = c(2, 6)) +
  theme_classic() +
  labs(
    title = "Game/fur animal lung, One-sample t-test (H1:mu> 0)",
    x = NULL,
    y = "ISG score difference (pos - neg)",
    color = "Host species"
  )

gf_one_l_inf
output_dir <- "."
gf_one_l_file_inf <- paste0(output_dir, "/", "gamefur_posneg_lung_onetest_infection_260205.pdf")
#ggsave(gf_one_l_file_inf, plot = gf_one_l_inf, height=5, width = 5)

