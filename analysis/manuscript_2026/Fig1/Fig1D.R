rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library(emmeans)
library(dplyr)
library(tidyr)
library(ggpubr)
library(stringr)
library(caret)
library(ggtree)

###this script see the ISG expression differences on experimental infection samples
###you should download the SRA detailed descriptions

ISG_meta <- read.table("exp_inf_260216.txt", sep = "\t", header = T)

#get specific BioProjectID data
#PRJEB9318
PRJEB9318_list <- read.table("infection/filereport_read_run_PRJEB9318_tsv.txt", sep = "\t", header = T) %>%
  rename(ID = run_accession) %>%
  mutate(
    virus = case_when(
      str_detect(sample_title, "inf") ~ "IBDV",
      str_detect(sample_title, "control") ~ "Negative",
      TRUE ~ NA_character_
    ),
    status = case_when(
      str_detect(sample_title, "inf") ~ "Positive",
      str_detect(sample_title, "control") ~ "Negative",
      TRUE ~ NA_character_  # その他はNA
    )
  ) %>%
  select(ID, virus, status, sample_title)

ISG_PRJEB9318_list <- PRJEB9318_list %>%
  inner_join(ISG_meta, by = "ID")

#PRJNA290579
PRJNA290579_list <- read.table("infection/filereport_read_run_PRJNA290579_tsv.txt", sep = "\t", header = T) %>%
  rename(ID = run_accession) %>%
  mutate(
    virus = case_when(
      str_detect(sample_title, "MEC_") ~ "PRRSV",
      str_detect(sample_title, "^INF_") ~ "PRRSV",
      str_detect(sample_title, "UNINF_") ~ "Negative",
      str_detect(sample_title, "CON_") ~ "Negative",
      TRUE ~ NA_character_
    ),
    status = case_when(
      str_detect(sample_title, "MEC_") ~ "Positive",
      str_detect(sample_title, "^INF_") ~ "Positive",
      str_detect(sample_title, "UNINF_") ~ "Negative",
      str_detect(sample_title, "CON_") ~ "Negative",
      TRUE ~ NA_character_  # その他はNA
    )
  ) %>%
  select(ID, virus, status, sample_title)

ISG_PRJNA290579_list <- PRJNA290579_list %>%
  inner_join(ISG_meta, by = "ID")

#PRJEB29321
PRJEB29321_list <- read.table("infection/filereport_read_run_PRJEB29321_tsv.txt", sep = "\t", header = T) %>%
  rename(ID = run_accession) %>%
  mutate(
    virus = case_when(
      str_detect(library_name, "Infected") ~ "MDV",
      str_detect(library_name, "Control") ~ "Negative",
      TRUE ~ NA_character_  # その他はNA
    ),
    status = case_when(
      str_detect(library_name, "Infected") ~ "Positive",
      str_detect(library_name, "Control") ~ "Negative",
      TRUE ~ NA_character_  # その他はNA
    )
  ) %>%
  select(ID, virus, status, library_name) %>%
  rename(sample_title = library_name)

ISG_PRJEB29321_list <- PRJEB29321_list %>%
  inner_join(ISG_meta, by = "ID")

#PRJNA303166
PRJNA303166_list <- read.table("infection/filereport_read_run_PRJNA303166_tsv.txt", sep = "\t", header = T) %>%
  rename(ID = run_accession) %>%
  mutate(
    virus = case_when(
      !str_detect(sample_title, "_D0_") ~ "PRRSV",
      str_detect(sample_title, "_D0_") ~ "Negative",
      TRUE ~ NA_character_  # その他はNA
    ),
    status = case_when(
      !str_detect(sample_title, "_D0_") ~ "Positive",
      str_detect(sample_title, "_D0_") ~ "Negative",
      TRUE ~ NA_character_  # その他はNA
    )
  ) %>%
  select(ID, virus, status, sample_title)

ISG_PRJNA303166_list <- PRJNA303166_list %>%
  inner_join(ISG_meta, by = "ID")

#PRJNA395043
PRJNA395043_list <- read.table("infection/filereport_read_run_PRJNA395043_tsv.txt", sep = "\t", header = T) %>%
  rename(ID = run_accession) %>%
  mutate(
    virus = case_when(
      str_detect(sample_title, "POWV") ~ "POWV",
      str_detect(sample_title, "control") ~ "Negative",
      TRUE ~ NA_character_
    ),
    status = case_when(
      str_detect(sample_title, "POWV") ~ "Positive",
      str_detect(sample_title, "control") ~ "Negative",
      TRUE ~ NA_character_  # その他はNA
    )
  ) %>%
  select(ID, virus, status, sample_title)

ISG_PRJNA395043_list <- PRJNA395043_list %>%
  inner_join(ISG_meta, by = "ID")

#PRJNA481895
PRJNA481895_list <- read.table("infection/filereport_read_run_PRJNA481895_tsv.txt", sep = "\t", header = T) %>%
  rename(ID = run_accession) %>%
  mutate(
    virus = case_when(
      str_detect(sample_title, "MARV") ~ "MARV",
      str_detect(sample_title, "SeV") ~ "SeV",
      str_detect(sample_title, "Mock") ~ "Negative",
      TRUE ~ NA_character_
    ),
    status = case_when(
      virus %in% c("MARV", "SeV") ~ "Positive",
      virus == "Negative" ~ "Negative",
      TRUE ~ NA_character_
    )
  ) %>%
  select(ID, virus, status, sample_title)

PRJNA481895_MARV_list <- PRJNA481895_list %>%
  inner_join(ISG_meta, by = "ID") %>%
  filter(virus == "MARV"| virus == "Negative") %>%
  mutate(BioProject_ID = if_else(BioProject_ID == "PRJNA481895",
                                 "PRJNA481895_MARV",
                                 BioProject_ID))


PRJNA481895_SeV_list <- PRJNA481895_list %>%
  filter(virus == "SeV"| virus == "Negative")

PRJNA481895_SeV_list <- PRJNA481895_list %>%
  inner_join(ISG_meta, by = "ID") %>%
  filter(virus == "SeV"| virus == "Negative") %>%
  mutate(BioProject_ID = if_else(BioProject_ID == "PRJNA481895",
                                 "PRJNA481895_SeV",
                                 BioProject_ID))


ISG_PRJNA481895_list <- PRJNA481895_list %>%
  inner_join(ISG_meta, by = "ID")

ISG_PRJ_comb <- rbind(ISG_PRJEB29321_list, ISG_PRJEB9318_list, ISG_PRJNA303166_list, PRJNA481895_MARV_list, PRJNA481895_SeV_list) %>%
  filter(!is.na(species))

ISG_PRJ_comb_mod <- ISG_PRJ_comb %>%
  ungroup() %>%
  group_by(ID, BioProject_ID, tax_id, kindom, phylum, class, order, family, genus, species, Host_species, virus, status) %>%
  summarise(mean_ISG = mean(ISG_score))

wilcox_results <- ISG_PRJ_comb_mod %>%
  group_by(BioProject_ID) %>%
  summarise(
    p_value = wilcox.test(mean_ISG ~ status)$p.value,
    .groups = "drop"
  )

score_box <- ggplot(ISG_PRJ_comb_mod, aes(x = status, y = mean_ISG, fill = status)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.5) +  # boxplot を薄く、線なしに
  geom_jitter(aes(color = status), shape = 21, color = "black", width = 0.2, alpha = 1.0, size = 1.5) +  # 点に色をつける
  stat_compare_means(aes(group = status), method = "t.test", label = "p.signif") +
  facet_grid(. ~ BioProject_ID) +
  xlab("") +
  ylab("ISG score") +
  labs(color="condition") +
  theme_classic() +
  scale_x_discrete(
    limits = c("Negative", "Positive"),
    labels = c("Negative" = "-", "Positive" = "+")  # ← ここで置換！
  ) +
  scale_fill_manual(values = c("Positive" = "#E41A1C", "Negative" = "#377EB8"))
score_box
output_dir <- "."
file_name <- paste0(output_dir, "/", "experiment_virus_251223.pdf")
#ggsave(file_name, plot = score_box, width = 5, height = 5)
