rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ggplot2")
library(tidytext)

setwd("~/OneDrive/ドキュメント/解析結果/Virome/virome/geNomad/Mammals_Aves_20k/")
#import data
Aves_Mam_mbio_genomad_filt <- read.table("Aves_Mam_mbio_genomad_filt_250630_for168438.txt", sep = "\t", header = T) %>%
  separate(taxonomy, into = c("Virus_g", "Realm_g", "Kingdom_g", "Phylum_g", "Class_g", "Order_g", "Family_g"), sep = ";", fill = "right") %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  select(-Virus_g)
ISG_meta <- read.table("Aves_Mam_mbio_250613_ISGscore_logan.txt", sep = "\t", header = T)

blastx <- read_tsv("blast/tophit_lin_mod.blastx.txt") %>%
  separate(tax_lin, into = c("Realm_b", "Kingdom_b", "Phylum_b", "Class_b", "Order_b", "Family_b", "Genus_b", "Species_b"), sep = ";", fill = "right") %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  group_by(seq_name) %>%
  slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  rename(bl_len = length) %>%
  mutate(
    Acc_ID = if_else(
      str_detect(full_lin, "unclassified_entries|Metazoa|Bacillati|Pseudomonadati|Viridiplantae|other_entries|cellular_organisms"),
      NA_character_,
      Acc_ID
    )
  )

genomad_seq <- Aves_Mam_mbio_genomad_filt %>%
  select(seq_name) %>%
  unique()
blastx_out <- read_tsv("blast/tophit_lin_mod.blastx.txt") %>%
  separate(tax_lin, into = c("Realm_b", "Kingdom_b", "Phylum_b", "Class_b", "Order_b", "Family_b", "Genus_b", "Species_b"), sep = ";", fill = "right") %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  separate(
    col = seq_name,
    into = c("run_ID", "number"),
    sep = "_",         
    remove = FALSE     
  ) %>%
  select(-number) %>%
  select(run_ID, seq_name, everything()) %>%
  inner_join(genomad_seq, by = "seq_name")

genomad_blastx <- Aves_Mam_mbio_genomad_filt %>%
  left_join(blastx, by = c("seq_name")) %>%
  mutate(across(where(is.character), ~ na_if(., "")))

herpes <- genomad_blastx %>%
  filter(Family_g == "Orthoherpesviridae")

herpes_nohit <- herpes %>%
  filter(is.na(Acc_ID))
herpes_nohit_seq <- herpes_nohit %>%
  select(seq_name) %>%
  unique()

herpes <- genomad_blastx %>%
  filter(Family_g == "Orthoherpesviridae")

herpes_nohit <- herpes %>%
  filter(is.na(Acc_ID))

test <- herpes_nohit %>%
  count(Species_b) %>%
  mutate(ratio = n / sum(n)) %>%
  mutate(classification = case_when(
    is.na(Species_b) ~ NA_character_,  # NAはそのまま
    str_detect(Species_b, regex("synthetic", ignore_case = TRUE)) ~ "synthetic",
    str_detect(Species_b, regex("vector|cloning", ignore_case = TRUE)) ~ "vector",
    TRUE ~ "cellular"
  ))

test2 <- test %>%
  group_by(classification) %>%
  summarise(cum_ratio = sum(ratio), cum_n = sum(n))

herpes_nohit_contig <- herpes_nohit %>%
  select(seq_name) %>%
  unique()

blastx_2602 <- read_tsv("blast/noherpes.blastx_lin.txt") %>%
  separate(tax_lin, into = c("Realm_b", "Kingdom_b", "Phylum_b", "Class_b", "Order_b", "Family_b", "Genus_b", "Species_b"), sep = ";", fill = "right") %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  group_by(seq_name) %>%
  slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%
  ungroup()

blastx_2602_mod <- herpes_nohit_contig %>%
  left_join(blastx_2602, by = "seq_name")

test_x <- blastx_2602_mod %>%
  count(Species_b) %>%
  mutate(ratio = n / sum(n)) %>%
  mutate(classification = case_when(
    is.na(Species_b) ~ NA_character_,  # NAはそのまま
    str_detect(Species_b, regex("synthetic", ignore_case = TRUE)) ~ "synthetic",
    str_detect(Species_b, regex("vector|cloning", ignore_case = TRUE)) ~ "vector",
    TRUE ~ "cellular"
  ))

test2_x <- test_x %>%
  group_by(classification) %>%
  summarise(cum_ratio = sum(ratio), cum_n = sum(n))

kingdom_x <- blastx_2602_mod %>%
  count(Kingdom_b) %>%
  mutate(ratio = n / sum(n))

blastx_2602_mod_sp <- blastx_2602_mod %>%
  separate(seq_name, into = c("ID", "number"), sep = "_", remove = FALSE) %>%
  left_join(ISG_meta, by = "ID")

tmp <- blastx_2602_mod_sp %>%
  filter(Kingdom_b == "Metazoa") %>%
  filter(Order_b == order)

order_x <- blastx_2602_mod %>%
  count(Kingdom_b, Order_b) %>%
  mutate(ratio = n / sum(n))

#write.table(order_x, "herpes_fp_taxcount.txt", row.names = F, sep = "\t", quote = F)

blast_herp <- blastx_2602_mod_sp %>%
  filter(Kingdom_b == "Heunggongvirae")

blastx_ori <- read_tsv("blast/noherpes.blastx.txt", col_names = FALSE)
colnames(blastx_ori) <- c(
  "seq_name","sseqid","pident","length","mismatch","gapopen",
  "qstart","qend","sstart","send","evalue","bitscore","strand", "taxid", "protein_name", "na1", "na2", "na3", "na4"
)
# 2. qseqidごとに bitscore 最大1件（tiesなし）
blastx_top <- blastx_ori %>%
  group_by(seq_name) %>%
  slice_max(order_by = bitscore, n = 1, with_ties = FALSE) %>%
  ungroup()

# 3. stitle加工（awk + cut + sed 相当）
blastx_protein <- blastx_top %>%
  mutate(
    protein = str_remove(protein_name, "^[^ ]+ "),
    protein = str_remove(protein, "\\[.*\\]"),
    protein = str_trim(protein),
    protein = str_replace_all(protein, " ", "_")
  ) %>%
  count(protein, sort = TRUE)

total <- 4752
# kingdom_x から数を取得
metazoa_n  <- kingdom_x %>% filter(Kingdom_b == "Metazoa") %>% pull(n)
herpes_n   <- kingdom_x %>% filter(Kingdom_b == "Heunggongvirae") %>% pull(n)
bacteria_n <- kingdom_x %>% filter(Kingdom_b == "Bacillati") %>% pull(n)
bacteria_plant_n <- kingdom_x %>%
  filter(Kingdom_b %in% c("Bacillati",
                          "Pseudomonadati",
                          "Viridiplantae")) %>%
  summarise(total = sum(n)) %>%
  pull(total)
nohit_n    <- kingdom_x %>% filter(is.na(Kingdom_b) | Kingdom_b == "NA") %>% pull(n)

collagen_n <- blastx_protein %>%
  filter(str_detect(protein, regex("collagen.*alpha.*3", ignore_case = TRUE))) %>%
  summarise(sum(n)) %>% pull()

bcl2_n <- blastx_protein %>%
  filter(str_detect(protein, regex("bcl", ignore_case = TRUE))) %>%
  summarise(sum(n)) %>% pull()

pabp_n <- blastx_protein %>%
  filter(str_detect(protein, regex("polyadenylate", ignore_case = TRUE))) %>%
  summarise(sum(n)) %>% pull()

metazoa_others_n <- metazoa_n - collagen_n - bcl2_n - pabp_n

pie_df <- tibble(
  category = c(
    "Metazoa: collagen alpha 3",
    "Metazoa: Bcl2",
    "Metazoa: polyadenylate binding protein",
    "Metazoa: Other protein",
    "Orthoerpesviridae",
    "Others",
    "No hit"
  ),
  n = c(
    collagen_n,
    bcl2_n,
    pabp_n,
    metazoa_others_n,
    herpes_n,
    bacteria_n,
    nohit_n
  )
) %>%
  mutate(
    perc = 100 * n / sum(n),
    category = fct_reorder(category, n, .desc = TRUE)
  )

fill_colors <- c(
  "Metazoa: collagen alpha 3" = "#4C72B0",
  "Metazoa: Bcl2" = "#55A868",
  "Metazoa: polyadenylate binding protein" = "#C44E52",
  "Metazoa: Other protein" = "#8172B3",
  "Orthoerpesviridae" = "#CCB974",
  "Others" = "grey75",   # 少し濃いめグレー
  "No hit" = "grey40"
)

pie_df <- pie_df %>%
  mutate(
    category = fct_reorder(category, n, .desc = TRUE),
    category = fct_relevel(category, "No hit", after = Inf)
  )

herpes_fp_pie <- ggplot(pie_df, aes(x = "", y = n, fill = category)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y", direction = -1) +
  scale_fill_manual(values = fill_colors) +
  theme_void() +
  theme(legend.position = "right") +
  labs(
    title = "Composition of 4,752 False-Positive Contigs",
    fill = NULL
  )
herpes_fp_pie
output_dir <- "."
herpes_fp_pie_name <- paste0(output_dir, "/", "herpes_fp_pie_260303.pdf")
#ggsave(herpes_fp_pie_name, plot = herpes_fp_pie, height=4, width = 4)
