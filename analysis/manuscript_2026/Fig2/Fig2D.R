rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ggplot2")
library(patchwork)
library(ggtext)
library(ggpubr)
library(multcompView)
library(tidytext)


#import data
Aves_Mam_mbio_genomad_filt <- read.table("Aves_Mam_mbio_genomad_filt_250630.txt", sep = "\t", header = T) %>%
  separate(taxonomy, into = c("Virus_g", "Realm_g", "Kingdom_g", "Phylum_g", "Class_g", "Order_g", "Family_g"), sep = ";", fill = "right") %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  select(-Virus_g)

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

library(forcats)
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
      TRUE ~ "Unassigned"
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
  "Undetected" = "#dcdddd",
  "Unmatched" = "#898989"
  #"Unassigned" = "grey70"
)

plot_family_num <- plot_family_num %>%
  mutate(identity_group = factor(
    identity_group,
    levels = c("≥95%", "90–95%", "70–90%", "<70%", "Undetected","Unmatched")
    #levels = c("≥95%", "90–95%", "70–90%", "<70%", "Unassigned")
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
bar_family_num_name <- paste0(output_dir, "/", "virfamily_barplot_260213.pdf")
#ggsave(bar_family_num_name, plot = bar_family_num, height=7, width = 7)
