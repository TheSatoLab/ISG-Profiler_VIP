rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library(dplyr)
library(tidyr)
library(ggplot2)

alignment <- read.csv("Alignment_coverage_20210519.txt", header = T, sep = "\t") %>%
  rename(ID = sra)
alignment_family <- alignment %>%
  filter(coverage >= 20) %>%
  select(ID, family) %>%
  unique() %>%
  rename(Family_k = family)

#count taxonomy
genomad_result <- read.table("all_44879.virus_cons_250312.tsv", sep = "\t", header = T) %>%
  separate(taxonomy, into = c("Virus_g", "Realm_g", "Kingdom_g", "Phylum_g", "Class_g", "Order_g", "Family_g"), sep = ";", fill = "right") %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  mutate(Family_g = if_else(is.na(Family_g), "Unassigned", Family_g)) %>%
  mutate(across(where(is.character), ~ replace_na(., "Unclassified"))) %>%
  select(-Virus_g)
metadata <- read.csv("mbio_Kawasaki_all_metadata_250123.txt", sep = "\t", header = T)

virus_all <- genomad_result %>%
  inner_join(metadata, by = "ID")

blast_sum <- metadata %>%
  ungroup() %>%
  group_by(blast_status) %>%
  summarise(sum_blast = n()) %>%
  ungroup()

virus_all_count <- virus_all %>%
  group_by(Family_g, blast_status) %>%
  summarise(sum_contig = sum(n_contig), sum_sample = n()) %>%
  inner_join(blast_sum, by = "blast_status") %>%
  ungroup() %>%
  mutate(rate_sample = sum_sample/sum_blast*100) %>%
  ungroup()

genomad_mod <- virus_all %>%
  filter(nucleotide == "RNA") %>%
  select(Family_g, ID, nucleotide) %>%
  unique()

## Kawasaki側: ID × Family
kawa_pairs <- alignment_family %>%
  filter(!is.na(Family_k)) %>%
  distinct(ID, Family = Family_k)

kawa_ID <- kawa_pairs %>%
  select(ID) %>%
  unique() %>%
  mutate(kawasaki = 1)

## geNomad側: ID × Family (+ nucleotideは取っておくと後で便利)
geno_pairs <- genomad_mod %>%
  filter(nucleotide == "RNA") %>%
  filter(!is.na(Family_g)) %>%
  distinct(ID, Family = Family_g, nucleotide)

#####RNA virus only metadata#####
geno_RNA <- geno_pairs %>%
  select(ID) %>%
  unique() %>%
  mutate(RNA_virus = 1)

metadata_mod_rna <- metadata %>%
  left_join(geno_RNA, by = "ID") %>%
  mutate(RNA_virus = if_else(is.na(RNA_virus), 0, RNA_virus))

# カテゴリ付け：factorに明示的に変換
metadata_plot <- metadata_mod_rna %>%
  mutate(
    Kawasaki = if_else(blast_status == "over20", "Kawasaki_Positive", "Kawasaki_Negative"),
    geNomad = if_else(RNA_virus == 1, "Positive", "Negative")
  ) %>%
  mutate(
    Kawasaki = factor(Kawasaki, levels = c("Kawasaki_Positive", "Kawasaki_Negative")),
    geNomad  = factor(geNomad, levels = c("Negative", "Positive"))  
  )

# 集計
summary_df <- metadata_plot %>%
  group_by(Kawasaki, geNomad) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(Kawasaki) %>%
  mutate(
    total = sum(n),
    percentage = round(n / total * 100, 1),
    label = paste0(n, "/", total, " (", percentage, "%)")
  )

# 可視化
geno_kawa <- ggplot(summary_df, aes(x = Kawasaki, y = percentage, fill = geNomad)) +
  geom_col() +
  scale_fill_manual(values = c("Positive" = "salmon", "Negative" = "grey80")) +
  labs(
    x = "Kawasaki et al. classification",
    y = "Positive ratio",
    fill = "geNomad result",
    title = "Comparison of RNA Virus Detection: geNomad vs Kawasaki et al."
  ) +
  theme_classic(base_size = 13)
geno_kawa
output_dir <- "."
geno_kawa_file_name <- paste0(output_dir, "/", "geno_kawa_bar_260204.pdf")
#ggsave(geno_kawa_file_name, plot = geno_kawa, width = 3, height = 3)
