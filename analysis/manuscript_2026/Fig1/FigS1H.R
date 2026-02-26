rm( list=ls(all=TRUE) ) # clean up R workspace
##load package
library("tidyverse")
library("ggplot2")
library(MASS)
library(viridis)
library(pheatmap) 
library(reshape2)
library(ggpubr)

ISGcntl_sal_RPM <- read.table("sp10_240801ref_ISGcntl_forin.txt", header = T, sep = "\t")

# ロング形式に変換
ISG_long <- ISGcntl_sal_RPM %>%
  filter(type == "ISG") %>%
  pivot_longer(cols = c(norm_kallisto, norm_kma, norm_salmon, norm_star_ort, norm_bowtie), names_to = "method", values_to = "normalized_value") %>%
  mutate(method = str_remove(method, "^norm_")) %>%
  mutate(reads_log = log2(normalized_value*(10E+5)+1))

# プロット（facetでツールごとに表示）
p <- ggplot(ISG_long, aes(x = RPM_log, y = reads_log)) +
  geom_point(aes(color = Induction), alpha = 0.7) +
  geom_smooth(method = "lm", se = FALSE, linetype = "dashed", color = "gray40") +
  stat_cor(method = "spearman", label.x = 0, label.y = 0.9, size = 3) +
  facet_grid(species ~ method + Induction, scales = "free") +
  theme_classic(base_size = 12) +
  labs(
    title = "Comparison of RPM_mean and normalized values across tools",
    x = "RPM_mean (reference)",
    y = "Normalized value"
  )
print(p)


# 対象のデータをフィルタリング（method == "salmon"）
salmon_data <- ISG_long %>%
  ungroup() %>%
  filter(method == "salmon") %>%
  group_by(hum_symbol, species, Induction) %>%
  summarise(RPM_mean = mean(RPM_log), sal_mean = mean(reads_log))

# アウトライアー遺伝子を指定
outliers <- c("HLA-A", "HLA-B", "HLA-C", "HLA-E", "HLA-F", "HLA-G", "EIF2AK2", "RNF19B")

# outlierかどうかの列を追加
salmon_data <- salmon_data %>%
  mutate(outlier_flag = ifelse(hum_symbol %in% outliers, "Outlier", "Normal"))

species_order <- c("Gallus gallus", "Homo sapiens", "Rattus norvegicus", "Myotis lucifugus",
                   "Pteropus vampyrus","Sus scrofa", "Ovis aries", "Bos taurus",
                   "Equus caballus", "Canis lupus")

salmon_data$species <- factor(salmon_data$species, levels = species_order)
# Induction も Positive → Negative の順に
salmon_data$Induction <- factor(salmon_data$Induction, levels = c("Positive", "Negative"))

salmon_data_mod <- salmon_data %>%
  mutate(gene_group = case_when(
    hum_symbol %in% c("HLA-A","HLA-B","HLA-C","HLA-E","HLA-F","HLA-G") ~ "HLA",
    hum_symbol == "EIF2AK2" ~ "EIF2AK2",
    hum_symbol == "RNF19B" ~ "RNF19B",
    TRUE ~ "Normal"
  )) %>%
  filter(species == "Rattus norvegicus")

color_values <- c(
  "Normal"   = "grey20",
  "HLA"      = "#D55E00",  # vermillion
  "EIF2AK2"  = "#0072B2",  # blue
  "RNF19B"   = "#009E73"   # green
)

salmon_data_mod <- salmon_data_mod %>%
  mutate(gene_group = factor(gene_group,
                             levels = c("Normal","HLA","EIF2AK2","RNF19B")))

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
  scale_color_manual(values = color_values) +
  coord_cartesian(xlim = c(0, 15), ylim = c(0, 20)) +
  facet_grid(rows = vars(Induction), cols = vars(species)) +
  theme_classic() +
  labs(
    title = "Scatter plot: STAR_log2_RPM_mean vs. salmon_normalized_count_mean",
    x = "STAR_log2_RPM_mean",
    y = "salmon_normalized_count_mean"
  )


print(scatter_plot_mod)

output_dir <- "."
scatter_plot_mod_file_name <- paste0(output_dir, "/", "outlier_240801rat_260218.pdf")
#ggsave(scatter_plot_mod_file_name, plot = scatter_plot_mod, width = 2.7, height = 3)


#make a scatter plot; star (with references) vs other methods (with custom DB)

ISG_reads_sp_mean_mod.list <- ISGcntl_sal_RPM %>%
  dplyr::select(-sum_kallisto, -sum_kma, -sum_salmon, -sum_star_ort, -sum_bowtie) %>%
  pivot_longer(cols = c(norm_kallisto, norm_kma, norm_salmon, norm_star_ort, norm_bowtie), names_to = "method", values_to = "norm_reads") %>%
  filter(type == "ISG")

ISG_reads_sp_mean_mod2.list <- ISG_reads_sp_mean_mod.list %>%
  mutate(Species = species) %>%
  unite(col=sp_con, Species, Induction, sep = "_")

ISG_reads_sp_mean_mod22.list <- ISG_reads_sp_mean_mod2.list %>%
  mutate(reads_log = log10(norm_reads*(10E+5)+1))


#get outlier information
ISG_reads_sp_mean_mod3.list <- ISG_reads_sp_mean_mod22.list %>%
  unite(col=sp_con_met, sp_con, method, sep = "_") %>%
  ungroup()
row_s.line <- length(unique(ISG_reads_sp_mean_mod3.list$sp_con_met))
outlier.list <- data.frame()
outlier_hs.list <- data.frame()
for (i in 1:row_s.line) {
  Species <- unique(ISG_reads_sp_mean_mod3.list$sp_con_met)[i]
  sample.i <- ISG_reads_sp_mean_mod3.list %>%
    dplyr::filter(sp_con_met == Species)
  sample.filt <- sample.i %>%
    filter(RPM_log >= 1) 
  x = sample.filt$RPM_log
  y = sample.filt$reads_log
  sample.out <- sample.filt %>%
    mutate(hatval_lm = hatvalues(lm(y ~ x)))
  sample.out <- sample.out %>%
    mutate(hatval_rlm = hatvalues(rlm(y ~ x)))
  influence <- lm.influence(lm(y ~ x))
  cooks_d <- cooks.distance(lm(y ~ x))
  leverage <- influence$hat
  sample.out$CooksD <- cooks_d
  sample.out$Leverage <- leverage
  sample.out$StudentizedResiduals <- studres(lm(y ~ x))
  robust_cov <- cov.rob(sample.out[, c("reads_log", "RPM_log")])
  
  ModelLM <- lm(y ~ x, data = sample.filt)
  fitted_intercept <- coef(ModelLM)[1]
  fitted_slope <- coef(ModelLM)[2]
  sample.out$finterc <- fitted_intercept
  sample.out$fslope <- fitted_slope
  
  # Mahalanobis距離の計算
  mahal_dist <- mahalanobis(sample.out[, c("reads_log", "RPM_log")], center = robust_cov$center, cov = robust_cov$cov)
  sample.out$MahalanobisDistance <- mahal_dist
  sample.out <- sample.out %>%
    dplyr::select(hum_symbol, hatval_lm, hatval_rlm, CooksD, Leverage, StudentizedResiduals, MahalanobisDistance, finterc, fslope)
  
  sample.i <- sample.i %>%
    left_join(sample.out, by = "hum_symbol")
  outlier_hs.list <- rbind(outlier_hs.list, sample.i)
  
}

outlier_mod.list <- outlier_hs.list %>%
  dplyr::select(hum_symbol, hatval_lm, hatval_rlm, sp_con_met, CooksD, Leverage, StudentizedResiduals, MahalanobisDistance, finterc, fslope) %>%
  dplyr::mutate(sp_con_met=gsub(sp_con_met,pattern="star_ort",replacement = "starort", ignore.case = TRUE)) %>%
  separate(col = sp_con_met, into = c("species", "condition", "norm", "met"), sep = "_") %>%
  unite(col= sp_con, species, condition, sep="_") %>%
  unite(col= method, norm, met, sep="_")

sal_out <- outlier_mod.list %>% filter(method == "norm_salmon") %>% group_by(hum_symbol) %>% summarise(stre_mean = mean(StudentizedResiduals, na.rm = TRUE ))
