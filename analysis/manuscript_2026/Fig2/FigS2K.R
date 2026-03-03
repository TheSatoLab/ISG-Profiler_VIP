rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library(dplyr)
library(ggplot2)
library(tidyverse)
library(ggpubr)

#import data
ISG_meta <- read.table("Aves_Mam_mbio_250613_ISGscore_logan.txt", sep = "\t", header = T)
Aves_Mam_mbio_genomad_filt_sum <- read.table("Aves_Mam_mbio_genomad_filt_ML_250630_for168438.txt", sep = "\t", header = T) %>% rename(ID = run_ID)

ISG_virus <- ISG_meta %>%
  filter(Infection == "Positive") %>%
  select(-ISG_level) %>%
  inner_join(Aves_Mam_mbio_genomad_filt_sum, by = c("ID"), relationship = "many-to-many") %>%
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

Herpes <- ISG_virus_neg_mod %>%
  filter(Family == "Orthoherpesviridae")

Herpes_host <- Herpes %>%
  select(species) %>%
  unique()

Herpes_neg <- ISG_virus_neg_mod %>%
  filter(Family == "Negative" & species %in% Herpes_host$species)

Herpes_pos_neg <- rbind(Herpes, Herpes_neg)

# 最終プロット
posneg_plot<- ggplot(Herpes_pos_neg, aes(x = Infection, y = ISG_mean, fill = Infection)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  coord_cartesian(ylim = c(-2, 2)) +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) +
  theme_classic() +
  labs(title = "Herpesviridae ISG score comparison with geNomad label",
       x = "Infection condition", y = "ISG score")

plot(posneg_plot)
output_dir <- "."
file_name <- paste0(output_dir, "/", "Herpes_box_all_260227.pdf")
#ggsave(file_name, plot = posneg_plot, width = 4, height = 4)


#import data
Aves_Mam_mbio_genomad_filt <- read.table("Aves_Mam_mbio_genomad_filt_250630_for168438.txt", sep = "\t", header = T) %>% rename(ID = run_ID)

genomad_fp <- read.table("blastx_herp_undetected.txt", sep = "\t", header = T)
Aves_Mam_mbio_genomad_filt_und <- Aves_Mam_mbio_genomad_filt %>%
  anti_join(genomad_fp, by = "seq_name")

lowISG_virus <- read.table("low_ISG_virus_260110.txt", sep = "\t", header = T) %>%
  mutate(ISG_level = "low")

Aves_Mam_mbio_genomad_filt_und_sum <- Aves_Mam_mbio_genomad_filt_und %>%
  group_by(ID, taxonomy, nucleotide) %>%
  summarise(n_contig = n(), .groups = "drop") %>%
  ungroup() %>%
  left_join(lowISG_virus, by = "taxonomy") %>%
  mutate(ISG_level = if_else(is.na(ISG_level), "high", ISG_level))

genomad_tp <- Aves_Mam_mbio_genomad_filt_und_sum %>%
  select(ID) %>%
  unique() %>%
  mutate(gen_tp = "Positive")

ISG_meta_tp <- ISG_meta %>%
  left_join(genomad_tp, by = "ID") %>%
  mutate(gen_tp = ifelse(is.na(gen_tp), "Negative", gen_tp))

ISG_virus_tp <- ISG_meta_tp %>%
  filter(gen_tp == "Positive") %>%
  select(-ISG_level) %>%
  inner_join(Aves_Mam_mbio_genomad_filt_und_sum, by = c("ID"), relationship = "many-to-many") %>%
  select(-n_contig)

ISG_neg_tp <- ISG_meta_tp %>%
  filter(gen_tp == "Negative") %>%
  mutate(taxonomy = "Negative", nucleotide = "Negative")

ISG_virus_neg_tp <- rbind(ISG_virus_tp, ISG_neg_tp) %>%
  select(-Infection) %>%
  rename(Infection = gen_tp)

ISG_virus_neg_mod_tp <- ISG_virus_neg_tp %>%
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

Herpes_tp <- ISG_virus_neg_mod_tp %>%
  filter(Family == "Orthoherpesviridae")

Herpes_host_tp <- Herpes_tp %>%
  select(species) %>%
  unique()

Herpes_neg_tp <- ISG_virus_neg_mod_tp %>%
  filter(Family == "Negative" & species %in% Herpes_host_tp$species)

Herpes_pos_neg_tp <- rbind(Herpes_tp, Herpes_neg_tp)

# 最終プロット
posneg_plot_tp <- ggplot(Herpes_pos_neg_tp, aes(x = Infection, y = ISG_mean, fill = Infection)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.7) +
  coord_cartesian(ylim = c(-2, 2)) +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8")) +
  theme_classic() +
  labs(title = "Herpesviridae ISG score comparison with geNomad & blastx label",
       x = "Infection condition", y = "ISG score")

plot(posneg_plot_tp)
output_dir <- "."
file_name_tp <- paste0(output_dir, "/", "Herpes_box_tp_260227.pdf")
#ggsave(file_name_tp, plot = posneg_plot_tp, width = 4, height = 4)
