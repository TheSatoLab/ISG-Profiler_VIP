rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ggplot2")
library("dplyr")


#import data
ISG_meta <- read.table("Aves_Mam_mbio_250613_ISGscore_logan.txt", header = T, sep = "\t")
amniota_list <- read.table("Amniota399_sp_id.list", sep = "\t", header = T)
library(taxize)
Sys.setenv(ENTREZ_KEY = "hoge")

taxids <- unique(amniota_list$tax_id)

# NCBI から階層情報を取得
tax_class <- classification(taxids, db = "ncbi")

# リストをデータフレーム化して genus, family, order だけ抽出
tax_df <- bind_rows(lapply(names(tax_class), function(id) {
  df <- tax_class[[id]]
  if (is.null(df)) return(data.frame(tax_id = id, genus = NA, family = NA, order = NA))
  tibble(
    tax_id = as.integer(id),
    genus  = df$name[df$rank == "genus"]  %||% NA,
    family = df$name[df$rank == "family"] %||% NA,
    order = df$name[df$rank == "order"] %||% NA
  )
})) %>%
  select(-tax_id) %>%
  unique()

Amni_ord <- tax_df %>%
  select(order) %>%
  unique() %>%
  mutate(DB_ord = "YES")

Amni_fam <- tax_df %>%
  select(family) %>%
  unique() %>%
  mutate(DB_fam = "YES")

Amni_gen <- tax_df %>%
  select(genus) %>%
  unique() %>%
  mutate(DB_gen = "YES")

ISG_meta_nonref <- ISG_meta %>%
  anti_join(amniota_list, by = "tax_id")

nonref_num <- ISG_meta_nonref %>%
  group_by(species, Infection, genus, family, order) %>%
  summarise(inf_num = n()) %>%
  ungroup() %>%
  group_by(species, genus, family, order) %>%
  summarise(
    total_num = sum(inf_num, na.rm = TRUE),
    pos_num   = sum(inf_num[Infection == "Positive"], na.rm = TRUE),
    neg_num   = sum(inf_num[Infection == "Negative"], na.rm = TRUE),
    infection_ratio = pos_num / total_num
  ) %>%
  ungroup()

amniota_list_mod <- amniota_list %>%
  mutate(species = sub("^([^_]+_[^_]+).*", "\\1", species))

include_list <- nonref_num %>%
  filter(!str_detect(species, "_sp$")) %>%
  filter(!str_detect(species, "_x_")) %>%
  filter(!str_detect(species, "_sp._")) %>%
  anti_join(amniota_list_mod, by = "species") %>%
  mutate(genus = sub("^([^_]+).*", "\\1", species)) %>%
  left_join(Amni_ord, by = "order") %>%
  left_join(Amni_fam, by = "family") %>%
  left_join(Amni_gen, by = "genus") %>%
  mutate(DB_gen = replace_na(DB_gen, "NO"), DB_fam = replace_na(DB_fam, "NO")) %>%
  filter(pos_num >= 2 & neg_num >= 2)

ISG_meta_nonref_many <- ISG_meta_nonref %>%
  filter(species %in% include_list$species)


PRJNA816878_meta <- read.table("filereport_read_run_PRJNA816878_tsv.txt", header = T, sep = "\t") %>%
  mutate(exp_inf = case_when(
    str_detect(sample_alias, regex("healthy", ignore_case = TRUE)) ~ "Negative",
    str_detect(sample_alias, regex("RSV", ignore_case = TRUE)) ~ "Positive",
    TRUE ~ NA_character_ 
  )) %>%
  mutate(tissue = case_when(
    str_detect(sample_alias, regex("intestin", ignore_case = TRUE)) ~ "Intestin",
    str_detect(sample_alias, regex("lung", ignore_case = TRUE)) ~ "Lung",
    str_detect(sample_alias, regex("kidney", ignore_case = TRUE)) ~ "kidney",
    str_detect(sample_alias, regex("heart", ignore_case = TRUE)) ~ "heart",
    str_detect(sample_alias, regex("spleen", ignore_case = TRUE)) ~ "spleen",
    TRUE ~ NA_character_  
  )) %>%
  select(run_accession, exp_inf, tissue) %>%
  rename(ID = run_accession)

PRJNA816878 <- ISG_meta %>%
  filter(BioProject_ID == "PRJNA816878") %>%
  inner_join(PRJNA816878_meta, by = "ID")  %>%
  mutate(species = if_else(species == "Sigmodon_fulviventer",
                           "Sigmodon_fulviventer_RSV",
                           species)) %>%
  mutate(species = if_else(species == "Sigmodon_hispidus",
                           "Sigmodon_hispidus_RSV",
                           species)) %>%
  mutate(ISG_DB = "Family") 

ggplot(PRJNA816878, aes(x = exp_inf, y = ISG_mean, fill = exp_inf)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.4) +
  facet_grid(Host_species ~ tissue) +
  scale_fill_manual(values = c("Negative" = "#0072B2", "Positive" = "#E69F00")) +
  theme_bw(base_size = 11) +
  labs(
    title = "ISG_mean by Infection Status, PRJNA816878",
    x = "Infection Status",
    y = "ISG_mean"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

#chrolocebus
PRJNA396021_meta <- read.table("filereport_read_run_PRJNA396021.tsv", header = T, sep = "\t") %>%
  mutate(exp_inf = case_when(
    str_detect(sample_title, regex("Uninfected", ignore_case = TRUE)) ~ "Negative",
    str_detect(sample_title, regex("negative", ignore_case = TRUE)) ~ "Negative",
    str_detect(sample_title, regex("positive", ignore_case = TRUE)) ~ "Positive",
    TRUE ~ NA_character_   
  )) %>%
  select(run_accession, exp_inf, library_strategy) %>%
  rename(ID = run_accession)

PRJNA396021 <- ISG_meta %>%
  filter(BioProject_ID == "PRJNA396021") %>%
  inner_join(PRJNA396021_meta, by = "ID")  %>%
  mutate(species = if_else(species == "Chlorocebus_aethiops",
                           "Chlorocebus_aethiops_HSV",
                           species))  %>%
  mutate(ISG_DB = "Genus")

ggplot(PRJNA396021, aes(x = exp_inf, y = ISG_mean, fill = exp_inf)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.4) +
  facet_wrap(~ species) +
  scale_fill_manual(values = c("Negative" = "#0072B2", "Positive" = "#E69F00")) +
  theme_bw(base_size = 11) +
  labs(
    title = "ISG_mean by Infection Status, PRJNA981676",
    x = "Infection Status",
    y = "ISG_mean"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )


#chrolocebus
PRJNA668354_meta <- read.table("filereport_read_run_PRJNA668354.tsv", header = T, sep = "\t") %>%
  mutate(exp_inf = case_when(
    str_detect(sample_title, regex("day 0", ignore_case = TRUE)) ~ "Negative",
    str_detect(sample_title, regex("day 1", ignore_case = TRUE)) ~ "Positive",
    str_detect(sample_title, regex("day 2", ignore_case = TRUE)) ~ "Positive",
    TRUE ~ NA_character_   
  )) %>%
  select(run_accession, exp_inf, library_strategy) %>%
  rename(ID = run_accession)

PRJNA668354 <- ISG_meta %>%
  filter(BioProject_ID == "PRJNA668354") %>%
  inner_join(PRJNA668354_meta, by = "ID") %>%
  mutate(species = if_else(species == "Chlorocebus_aethiops",
                                  "Chlorocebus_aethiops_SCV2",
                                  species)) %>%
  mutate(ISG_DB = "Genus")

ggplot(PRJNA668354, aes(x = exp_inf, y = ISG_mean, fill = exp_inf)) +
  geom_boxplot(alpha = 0.8, outlier.shape = NA) +
  geom_jitter(width = 0.2, size = 1.2, alpha = 0.4) +
  facet_wrap(~ species) +
  scale_fill_manual(values = c("Negative" = "#0072B2", "Positive" = "#E69F00")) +
  theme_bw(base_size = 11) +
  labs(
    title = "ISG_mean by Infection Status, PRJNA981676",
    x = "Infection Status",
    y = "ISG_mean"
  ) +
  theme(
    legend.position = "none",
    axis.text.x = element_text(angle = 45, hjust = 1)
  )

PRJNA816878_mod <- PRJNA816878 %>%
  select(-tissue)

PRJNA668354_mod <- PRJNA668354 %>%
  select(-library_strategy)

PRJNA396021_mod <- PRJNA396021 %>%
  select(-library_strategy)


ISG_PRJ_comb <- rbind(PRJNA816878_mod, PRJNA668354_mod, PRJNA396021_mod) %>%
  filter(!is.na(species))

ISG_PRJ_comb_mod <- ISG_PRJ_comb %>%
  ungroup() %>%
  group_by(ID, BioProject_ID, tax_id, kingdom, phylum, class, order, family, genus, species, Host_species, exp_inf, ISG_DB) %>%
  summarise(mean_ISG = mean(ISG_mean))

#write.table(ISG_PRJ_comb_mod, "Experimental_infection_noISGDB.txt", row.names = F, sep = "\t", quote = F)

wilcox_results <- ISG_PRJ_comb_mod %>%
  group_by(species) %>%
  summarise(
    p_value = wilcox.test(mean_ISG ~ exp_inf)$p.value,
    .groups = "drop"
  )

library(ggpubr)
ISG_DB_list <- ISG_PRJ_comb_mod %>%
  ungroup() %>%
  select(species, ISG_DB) %>%
  unique()

# 1列のheatmapとして描画
heat_annot <- ggplot(ISG_DB_list, aes(x = species, y = 1, fill = ISG_DB)) +
  geom_tile(color = "white") +
  scale_fill_manual(values = c("Genus" = "darkgreen", "Family" = "orange")) +
  theme_void() +
  theme(
    legend.position = "bottom",
    axis.text.x = element_blank(),
    plot.margin = margin(0, 5, 0, 5)
  ) +
  labs(fill = "ISG_DB Level")

score_box <- ggplot(ISG_PRJ_comb_mod, aes(x = exp_inf, y = mean_ISG, fill = exp_inf)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  
  geom_jitter(aes(color = exp_inf), shape = 21, color = "black", width = 0.2, alpha = 1.0, size = 1.5) +  
  stat_compare_means(aes(group = exp_inf), method = "t.test", label = "p.signif") +
  facet_grid(. ~ species) +
  xlab("") +
  ylab("Z-score of normalized reads count") +
  labs(color="condition", fill = "Infection") +
  theme_classic() +
  scale_x_discrete(limit = c("Negative", "Positive")) +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) +
  scale_x_discrete(
    limits = c("Negative", "Positive"),
    labels = c("Negative" = "-", "Positive" = "+")
  )

final_box <-  score_box + heat_annot +
  plot_layout(heights = c(8, 1), guides = "collect")
plot(final_box)

score_box
output_dir <- "."
file_name <- paste0(output_dir, "/", "noISGDB_expinf_boxplot_251224.pdf")
#ggsave(file_name, plot = score_box, width = 5, height = 5)
