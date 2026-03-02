rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ggplot2")

#import data
Aves_Mam_mbio_genomad_filt_sum <- read.table("Aves_Mam_mbio_genomad_filt_ML_250630.txt", sep = "\t", header = T) %>%
  separate(taxonomy, into = c("Virus_g", "Realm_g", "Kingdom_g", "Phylum_g", "Class_g", "Order_g", "Family_g"), sep = ";", fill = "right") %>%
  mutate(across(where(is.character), ~ na_if(., ""))) %>%
  select(-Virus_g)

# 4. Top 15 の Family_g 抽出
top15_df <- Aves_Mam_mbio_genomad_filt_sum %>%
  filter(!is.na(Family_g)) %>%
  count(Family_g, sort = TRUE) %>%
  slice_head(n = 15)

plot_family_num <- Aves_Mam_mbio_genomad_filt_sum %>%
  filter(Family_g %in% top15_df$Family_g) %>%
  count(Family_g, nucleotide)

plot_family_num <- plot_family_num %>%
  mutate(Family_g = factor(Family_g, levels = rev(top15_df$Family_g)))

# 7. 描画
bar_family_num <- ggplot(plot_family_num, aes(x = n, y = Family_g, fill = nucleotide)) +
  geom_col(position = "stack") +
  labs(
    x = "Number of samples",
    y = "geNomad Viral Family",
    title = "Top 15 Viral Families",
    fill = "Virus type"
  ) +
  scale_fill_manual(values = c("DNA" = "#619CFF", "RNA" = "salmon")) +
  theme_classic(base_size = 12)

bar_family_num

output_dir <- "."
bar_family_num_name <- paste0(output_dir, "/", "virfamily_barplot_persample_simple_260227.pdf")
#ggsave(bar_family_num_name, plot = bar_family_num, height=7, width = 7)

