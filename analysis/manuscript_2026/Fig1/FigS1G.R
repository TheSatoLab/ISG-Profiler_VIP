rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ComplexHeatmap")
library("circlize")
library("ggplot2")
library(psych)

gene2ref <- read.table("gene2refseq_Amniota_ISGcntl_260205.list", sep = "\t", header = T) %>%
  select(-Isoform) %>%
  unique()

RPM_ISGcntl.list <- read.table("ISG_cntl_analyses/RPM_mean_comb_241003.list", header = T, sep="\t") %>%
  rename(tax_id = taxid) %>%
  inner_join(gene2ref, by = c("tax_id", "Other_symbol")) %>%
  select(ID, RPM, hum_symbol) %>%
  rename(RPM_log = RPM)

sp10_ISGcntl_meta_mod <- read.table("Aves_Mam_mbio_260205ref_ISGcntl_forin.txt", header = T, sep = "\t")

ISG_RPM <- sp10_ISGcntl_meta_mod %>%
  full_join(RPM_ISGcntl.list, c("ID", "hum_symbol"))

ISG_RPM_mean <- ISG_RPM  %>%
  filter(type == "ISG") %>%
  group_by(hum_symbol, species, Induction) %>%
  summarise(RPM_mean = mean(RPM_log), sal_mean = mean(norm_salmon_log))
# アウトライアー遺伝子を指定
outliers <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "EIF2AK2", "RNF19B")
species_order <- c("Gallus gallus", "Homo sapiens", "Rattus norvegicus", "Myotis lucifugus",
                   "Pteropus vampyrus","Sus scrofa", "Ovis aries", "Bos taurus",
                   "Equus caballus", "Canis lupus")

# outlierかどうかの列を追加
ISG_RPM_mean <- ISG_RPM_mean %>%
  mutate(outlier_flag = ifelse(hum_symbol %in% outliers, "Outlier", "Normal"))

ISG_RPM_mean$species <- factor(ISG_RPM_mean$species, levels = species_order)
# Induction も Positive → Negative の順に
ISG_RPM_mean$Induction <- factor(ISG_RPM_mean$Induction, levels = c("Positive", "Negative"))

salmon_data_mod <- ISG_RPM_mean %>%
  mutate(gene_group = case_when(
    hum_symbol %in% c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G") ~ "HLA",
    hum_symbol == "EIF2AK2" ~ "EIF2AK2",
    hum_symbol == "RNF19B" ~ "RNF19B",
    TRUE ~ "Normal"
  ))

color_values <- c(
  "Normal"   = "grey20",
  "HLA"      = "#D55E00",  # vermillion
  "EIF2AK2"  = "#0072B2",  # blue
  "RNF19B"   = "#009E73"   # green
)

salmon_data_mod <- salmon_data_mod %>%
  mutate(gene_group = factor(gene_group,
                             levels = c("HLA","EIF2AK2","RNF19B", "Normal")))

scatter_plot_mod <- ggplot() +
  geom_point(
    data = salmon_data_mod %>% filter(gene_group == "Normal"),
    aes(x = RPM_mean, y = sal_mean),
    color = "grey50", alpha = 0.5, size = 2
  ) +
  geom_point(
    data = salmon_data_mod %>% filter(gene_group != "Normal"),
    aes(x = RPM_mean, y = sal_mean, color = gene_group),
    size = 2
  ) +
  geom_smooth(
    data = salmon_data_mod %>% filter(gene_group == "Normal"),
    aes(x = RPM_mean, y = sal_mean),
    method = "lm",
    se = FALSE,
    color = "black"
  ) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 20)) +
  scale_color_manual(values = color_values) +
  facet_grid(rows = vars(Induction), cols = vars(species)) +
  theme_classic() +
  labs(
    title = "Scatter plot: STAR_log2_RPM_mean vs. salmon_normalized_count_mean",
    x = "STAR_log2_RPM_mean",
    y = "salmon_normalized_count_mean"
  )

print(scatter_plot_mod)

output_dir <- "."
scatter_plot_mod_file_name <- paste0(output_dir, "/", "outlier_260205_260218.pdf")
#ggsave(scatter_plot_mod_file_name, plot = scatter_plot_mod, width = 10, height = 3)
