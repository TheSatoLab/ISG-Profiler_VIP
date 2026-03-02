rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ggplot2")
library(patchwork)
library(ggtext)
library(ggpubr)
library(multcompView)
library(tidytext)
library(forcats)

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

genomad_blastx <- Aves_Mam_mbio_genomad_filt %>%
  left_join(blastx, by = c("seq_name")) %>%
  mutate(across(where(is.character), ~ na_if(., "")))

# 1. Identity bin の定義関数
identity_bin <- function(x) {
  case_when(
    x >= 95 ~ "≥95%",
    x >= 90 ~ "90–95%",
    x >= 70 ~ "70–90%",
    TRUE    ~ "<70%"
  )
}

blast_labelled <- genomad_blastx %>%
  mutate(
    ## detection level
    detection = if_else(is.na(Acc_ID), "Undetected", "Detected"),
    matched_status = case_when(
      detection == "Undetected" ~ "Undetected",
      Family_g == Family_b & !is.na(Family_b) ~ "Matched",
      TRUE ~ "Unmatched"
    ),
    ## identity_group: only for Matched
    identity_group = case_when(
      detection == "Undetected" ~ "Undetected",
      matched_status == "Unmatched" ~ "Unmatched",
      matched_status == "Matched" & Identity >= 95 ~ "≥95%",
      matched_status == "Matched" & Identity >= 90 ~ "90–95%",
      matched_status == "Matched" & Identity >= 70 ~ "70–90%",
      matched_status == "Matched" & Identity < 70  ~ "<70%",
      TRUE ~ NA_character_
    ),
    identity_group = factor(identity_group,
                            levels = c("≥95%", "90–95%", "70–90%", "<70%", "Unmatched", "Undetected")
    )
  )

# 3. seq_name × Family_g × identity_group で unique count
family_id_count <- blast_labelled %>%
  distinct(seq_name, Family_g, identity_group) %>%
  group_by(Family_g, identity_group) %>%
  summarise(n = n(), .groups = "drop")

# 4. Top 15 の Family_g 抽出
top20_family <- family_id_count %>%
  group_by(Family_g) %>%
  summarise(total = sum(n)) %>%
  filter(!is.na(Family_g)) %>%
  arrange(desc(total)) %>%
  slice_head(n = 15) %>%
  pull(Family_g)

# 5. プロット用データ作成
plot_family_num <- family_id_count %>%
  filter(Family_g %in% top20_family) %>%
  mutate(Family_g = fct_reorder(Family_g, n, .fun = sum))  

# 6. カラーパレット設定
fill_colors <- c(
  "≥95%" = "#171c61",
  "90–95%" = "#006d2c",
  "70–90%" = "#a5c675",
  "<70%" = "#fde725",
  "Undetected" = "#898989",
  "Unmatched" = "#dcdddd"
)

plot_family_num <- plot_family_num %>%
  mutate(identity_group = factor(
    identity_group,
    levels = c("≥95%", "90–95%", "70–90%", "<70%", "Unmatched", "Undetected")
  ))
# 7. 描画
bar_family_num <- ggplot(plot_family_num, aes(x = n, y = Family_g, fill = identity_group)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = fill_colors, name = "Identity") +
  labs(
    x = "Number of viral contigs",
    y = "geNomad Viral Family",
    title = "Top 15 Viral Families by Blastx Identity Bin",
    fill = "Identity"
  ) +
  theme_classic(base_size = 12)

bar_family_num

output_dir <- "."
bar_family_num_name <- paste0(output_dir, "/", "virfamily_barplot_percontig_260225.pdf")
#ggsave(bar_family_num_name, plot = bar_family_num, height=7, width = 7)

#coronavirus
Coronaviridae <- genomad_blastx %>%
  filter(Family_b == "Coronaviridae")

Corona_Ident <- Coronaviridae %>%
  group_by(run_ID, Species_b) %>%
  summarise(mean_ident = mean(Identity), bit_max = max(bitscore)) %>%
  ungroup()

tmp2 <- Corona_Ident %>%
  filter(mean_ident >= 95) %>%
  unique()


# データ集計
pie_data <- tmp2 %>%
  count(Species_b) %>%
  arrange(desc(n)) %>%
  mutate(
    perc = round(100 * n / sum(n), 1),
    label = paste0(Species_b, " (", perc, "%)")
  )

pie_data_top <- pie_data %>%
  mutate(Species_b = ifelse(row_number() > 5, "Other", Species_b)) %>%
  group_by(Species_b) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(
    Species_b = fct_reorder(Species_b, n, .desc = TRUE)
  ) %>%
  mutate(
    Species_b = fct_relevel(Species_b, "Other", after = Inf)
  ) %>%
  mutate(
    perc = round(100 * n / sum(n), 1),
    label = paste0(Species_b, " (", perc, "%)")
  )

# カラーパレット
fill_colors <- c(
  "#0072B2",  # blue
  "#D55E00",  # vermillion
  "#009E73",  # bluish green
  "#CC79A7",  # reddish purple
  "#E69F00",  # orange
  "grey70"    # Other
)

# Pie chart
virus_pie <- ggplot(pie_data_top, aes(x = "", y = -n, fill = Species_b)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = fill_colors) +
  labs(title = "Composition of Species", fill = "Species") +
  theme_classic() +
  theme(legend.position = "right")
virus_pie

output_dir <- "."
virus_pie_name <- paste0(output_dir, "/", "corona_vir_pie_260227.pdf")
#ggsave(virus_pie_name, plot = virus_pie, height=4, width = 4)

tmp3 <- Corona_Ident %>%
  filter(mean_ident >= 95) %>%
  select(run_ID) %>%
  unique()

tmp4 <- Corona_Ident %>%
  filter(mean_ident >= 95) %>%
  filter(Species_b == "Betacoronavirus_pandemicum") %>%
  select(run_ID) %>%
  unique()

#betacoronavirus pandemicum infected host
host <- tmp4 %>%
  rename(ID = run_ID) %>%
  inner_join(ISG_meta, by = "ID")

host_pie_data <- host %>%
  filter(!is.na(species)) %>%
  count(species, sort = TRUE) %>%  
  mutate(
    species = ifelse(row_number() > 5, "Other", species)
  ) %>%
  group_by(species) %>%
  summarise(n = sum(n), .groups = "drop") %>%
  mutate(
    species = fct_reorder(species, n, .desc = TRUE)
  ) %>%
  mutate(
    species = fct_relevel(species, "Other", after = Inf)
  ) %>%
  mutate(
    perc = round(100 * n / sum(n), 1),
    label = paste0(species, " (", perc, "%)")
  )

# カラーパレット
fill_colors_host <- c(
  "#4C72B0",  # muted blue
  "#55A868",  # muted green
  "#C44E52",  # muted red
  "#8172B3",  # muted purple
  "#CCB974",  # muted yellow-brown
  "grey75"    # Other
)

# Pie chart
host_pie <- ggplot(host_pie_data, aes(x = "", y = -n, fill = species)) +
  geom_col(width = 1, color = "white") +
  coord_polar(theta = "y") +
  scale_fill_manual(values = fill_colors_host) +
  labs(title = "Composition of Host species", fill = "Species") +
  theme_classic() +
  theme(legend.position = "right")

host_pie

output_dir <- "."
host_pie_name <- paste0(output_dir, "/", "betacovp_host_pie_260227.pdf")
#ggsave(host_pie_name, plot = host_pie, height=4, width = 4)
