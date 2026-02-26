rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ggplot2")
library("dplyr")


#import data
Aves_Mam_mbio_meta_filt <- read.table("Aves_Mam_mbio_metadata_ISG_logan_ML_250616_2.txt", sep = "\t", header = T)
Aves_Mam_mbio_genomad_summary_filt <- read.table("Aves_Mam_mbio_genomad_filt_ML_250630.txt", sep = "\t", header = T)
ISG_meta <- read.table("Aves_Mam_mbio_250613_ISGscore_logan.txt", header = T, sep = "\t")
ISG_virus <- ISG_meta %>%
  filter(Infection == "Positive") %>%
  select(-ISG_level) %>%
  inner_join(Aves_Mam_mbio_genomad_summary_filt, by = c("ID"), relationship = "many-to-many") %>%
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


ISG_summary <- ISG_virus_neg_mod %>%
  mutate(ISG_level = tolower(ISG_level)) %>%
  group_by(ID) %>%
  summarise(
    final_ISG_level = case_when(
      any(ISG_level == "high") ~ "high",
      any(ISG_level %in% c("low", "negative")) ~ "low",
      TRUE ~ NA_character_
    ),
    .groups = "drop"
  ) %>%
  mutate(
    Inf_ISG = case_when(
      final_ISG_level == "high" ~ "Positive",
      final_ISG_level %in% c("low", "negative") ~ "Negative",
      TRUE ~ NA_character_
    )
  )

# 元のデータに結合
ISG_virus_neg_mod_annot <- ISG_virus_neg_mod %>%
  left_join(ISG_summary, by = "ID")

ISG_virus_neg_mod_ID <- ISG_virus_neg_mod_annot %>%
  select(-Virus, -Realm, -Kingdom, -Phylum, -Class, -Order, -Family, -nucleotide, -ISG_level, -final_ISG_level) %>%
  unique()

###ISG Database animal taxonomy###
amniota_list <- read.table("Amniota399_sp_id.list", sep = "\t", header = T)
library(taxize)
Sys.setenv(ENTREZ_KEY = "hoge")

# taxid のリスト
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

Amni_sp <- amniota_list %>%
  select(species) %>%
  unique() %>%
  mutate(DB_sp = "YES")

all_num <- ISG_virus_neg_mod_ID %>%
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
  mutate(species = sub("^([^_]+_[^_]+).*", "\\1", species)) %>%
  mutate(DB_sp = "YES") %>%
  select(-tax_id) %>%
  unique()

all_list <- all_num %>%
  filter(!str_detect(species, "_sp$")) %>%
  filter(!str_detect(species, "_x_")) %>%
  filter(!str_detect(species, "_sp._")) %>%
  full_join(amniota_list_mod, by = "species") %>%
  mutate(genus = sub("^([^_]+).*", "\\1", species)) %>%
  left_join(Amni_ord, by = "order") %>%
  left_join(Amni_fam, by = "family") %>%
  left_join(Amni_gen, by = "genus") %>%
  mutate(DB_gen = replace_na(DB_gen, "NO"), DB_fam = replace_na(DB_fam, "NO"), DB_sp = replace_na(DB_sp, "NO"), DB_ord = replace_na(DB_ord, "NO")) %>%
  filter(!is.na(total_num)) %>%
  mutate(
    ISG_DB = case_when(
      DB_sp  == "YES" ~ "species",
      DB_gen == "YES" ~ "genus",
      DB_fam == "YES" ~ "family",
      DB_ord == "YES" ~ "order",
      TRUE ~ "class"
    )
  )

all_part_list <- all_num %>%
  filter(str_detect(species, "_sp$") | str_detect(species, "_x_") | str_detect(species, "_sp._")) %>%
  mutate(DB_sp = "NO") %>%
  left_join(Amni_ord, by = "order") %>%
  left_join(Amni_fam, by = "family") %>%
  left_join(Amni_gen, by = "genus") %>%
  mutate(DB_gen = replace_na(DB_gen, "NO"), DB_fam = replace_na(DB_fam, "NO"), DB_sp = replace_na(DB_sp, "NO"), DB_ord = replace_na(DB_ord, "NO")) %>%
  mutate(
    ISG_DB = case_when(
      DB_sp  == "YES" ~ "species",
      DB_gen == "YES" ~ "genus",
      DB_fam == "YES" ~ "family",
      DB_ord == "YES" ~ "order",
      TRUE ~ "class"
    )
  ) %>%
  mutate(DB_sp = if_else(
    str_detect(species, "Bos_indicus_x_Bos_taurus"),
    "YES",
    DB_sp
  ))

#種レベルでISG_DBにどの階層であるのか
species_level <- rbind(all_list, all_part_list) %>%
  distinct(species, ISG_DB) %>%     
  count(ISG_DB) %>%
  mutate(
    ISG_DB = factor(ISG_DB, levels = c("class","order","family","genus","species")),
    pct = n / sum(n) * 100
  )

sp_level_p <- ggplot(species_level, aes(x = "", y = pct, fill = ISG_DB)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  scale_fill_brewer(palette = "BuGn", direction = -1) +
  labs(
    title = "Reference genome level availability per species (n = 915)",
    x = NULL,
    y = "Percentage of species",
    fill = "Reference Level"
  ) +
  theme_classic(base_size = 12)

sp_level_p
output_dir <- "."
sp_level_p_name <- paste0(output_dir, "/", "ISGDB_species_260122.pdf")
#ggsave(sp_level_p_name, plot = sp_level_p, height=3, width = 4)

#サンプル数レベルで
sample_level <- rbind(all_list, all_part_list) %>%
  group_by(ISG_DB) %>%
  summarise(sample_n = sum(total_num)) %>%
  mutate(
    ISG_DB = factor(ISG_DB, levels = c("class","order","family","genus","species")),
    pct = sample_n / sum(sample_n) * 100
  )

sam_level_p <- ggplot(sample_level, aes(x = "", y = pct, fill = ISG_DB)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  scale_fill_brewer(palette = "BuGn",  direction = -1) +
  labs(
    title = "Reference genome level availability per sample",
    y = "Percentage of samples",
    x = NULL,
    fill = "Reference Level"
  ) +
  theme_classic(base_size = 12)
sam_level_p
output_dir <- "."
sam_level_p_name <- paste0(output_dir, "/", "ISGDB_sample_260122.pdf")
#ggsave(sam_level_p_name, plot = sam_level_p, height=3, width = 4)

###ISG Database animal taxonomy###
taxonomy_2506 <- read.csv("amniota_40905_taxrank.list", header = F) %>%
  select(-V1, -V2, -V3) %>%
  rename(species = V8, genus = V7, family = V6, order = V5, class = V4)

tax_398_2506 <- read.csv("amniota_398_taxrank.list", header = F) %>%
  select(-V1, -V2, -V3) %>%
  rename(species = V8, genus = V7, family = V6, order = V5, class = V4)

test <- tax_398_2506 %>% select(species) %>% left_join(taxonomy_2506) %>% select(species) %>% unique() %>% mutate(exist = "YES")
test2 <- left_join(tax_398_2506, test, by = "species")


# ユニークなカバー対象
species_set <- unique(tax_398_2506$species)
genus_set   <- unique(tax_398_2506$genus)
family_set  <- unique(tax_398_2506$family)
order_set   <- unique(tax_398_2506$order)

# 各階層で一致する行だけ抽出
species_cov <- taxonomy_2506 %>% filter(species %in% species_set)
genus_cov   <- taxonomy_2506 %>% filter(!(species %in% species_set) & genus %in% genus_set)
family_cov  <- taxonomy_2506 %>% filter(!(species %in% species_set) &
                                          !(genus %in% genus_set) &
                                          family %in% family_set)
order_cov   <- taxonomy_2506 %>% filter(!(species %in% species_set) &
                                          !(genus %in% genus_set) &
                                          !(family %in% family_set) &
                                          order %in% order_set)

# 未カバー種数
covered_species <- nrow(species_cov)
covered_genus   <- nrow(genus_cov)
covered_family  <- nrow(family_cov)
covered_order   <- nrow(order_cov)
uncovered       <- nrow(taxonomy_2506) - (covered_species + covered_genus + covered_family + covered_order)

# 結果のデータフレーム作成
coverage_df <- tibble(
  level = c("Species", "Genus", "Family", "Order", "Class"),
  count = c(covered_species, covered_genus, covered_family, covered_order, uncovered)
) %>%
  mutate(percent = count / nrow(taxonomy_2506))

coverage_df <- coverage_df %>%
  mutate(level = factor(level, c("Class", "Order", "Family", "Genus", "Species")))

# 積み上げ棒グラフ
ggplot(coverage_df, aes(y = level, x = percent, fill = level)) +
  geom_bar(stat = "identity") +
  scale_x_continuous(labels = scales::percent_format(accuracy = 1)) +
  scale_fill_brewer(palette = "Set2", name = "Taxonomic Level") +
  theme_minimal(base_size = 12) +
  labs(
    x = "Coverage (%)",
    y = NULL,
    title = "Taxonomic Coverage of RefSeq Amniota by ISG DB"
  ) +
  theme(
    axis.text.y = element_text(size = 10),
    legend.position = "none"  
  )

sam_level_p <- ggplot(coverage_df, aes(x = "", y = percent, fill = level)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  scale_fill_brewer(palette = "BuGn",  direction = -1) +
  labs(
    title = "Reference genome level availability per sample",
    y = "Percentage of samples",
    x = NULL,
    fill = "Reference Level"
  ) +
  theme_classic(base_size = 12)
sam_level_p

output_dir <- "."
sam_level_p_name <- paste0(output_dir, "/", "ISGDB_refseq_260120.pdf")
#ggsave(sam_level_p_name, plot = sam_level_p, height=3, width = 4)

taxonomy_2506_norep <- taxonomy_2506 %>%
  filter(class != "Lepidosauria" & order != "Crocodylia" & order != "Lepidosauria" & order != "Testudines" & family != "Hadrosauridae" & family != "Tyrannosauridae")

# 各階層で一致する行だけ抽出
species_cov_nr <- taxonomy_2506_norep %>% filter(species %in% species_set)
genus_cov_nr   <- taxonomy_2506_norep %>% filter(!(species %in% species_set) & genus %in% genus_set)
family_cov_nr  <- taxonomy_2506_norep %>% filter(!(species %in% species_set) &
                                                   !(genus %in% genus_set) &
                                                   family %in% family_set)
order_cov_nr   <- taxonomy_2506_norep %>% filter(!(species %in% species_set) &
                                                   !(genus %in% genus_set) &
                                                   !(family %in% family_set) &
                                                   order %in% order_set)

# 未カバー種数
covered_species_nr <- nrow(species_cov_nr)
covered_genus_nr   <- nrow(genus_cov_nr)
covered_family_nr  <- nrow(family_cov_nr)
covered_order_nr   <- nrow(order_cov_nr)
uncovered_nr       <- nrow(taxonomy_2506_norep) - (covered_species_nr + covered_genus_nr + covered_family_nr + covered_order_nr)

# 結果のデータフレーム作成
coverage_df_nr <- tibble(
  level = c("Species", "Genus", "Family", "Order", "Class"),
  count = c(covered_species_nr, covered_genus_nr, covered_family_nr, covered_order_nr, uncovered_nr)
) %>%
  mutate(percent = count / nrow(taxonomy_2506_norep))

coverage_df_nr <- coverage_df_nr %>%
  mutate(level = factor(level, c("Class", "Order", "Family", "Genus", "Species")))

# 積み上げ棒グラフ
sam_level_p_nr <- ggplot(coverage_df_nr, aes(x = "", y = percent, fill = level)) +
  geom_bar(stat = "identity", width = 0.6) +
  coord_flip() +
  scale_fill_brewer(palette = "BuGn",  direction = -1) +
  labs(
    title = "Reference genome level availability per sample",
    y = "Percentage of samples",
    x = NULL,
    fill = "Reference Level"
  ) +
  theme_classic(base_size = 12)
sam_level_p_nr

output_dir <- "."
sam_level_p_nr_name <- paste0(output_dir, "/", "ISGDB_refseq_norep_260120.pdf")
#ggsave(sam_level_p_nr_name, plot = sam_level_p_nr, height=3, width = 4)
