library(readr)
library(dplyr)
library(tidyr)
library(ggplot2)
library(patchwork)
library(RColorBrewer)

# change directory
setwd("/Users/kyokokurihara/iLab/itolab_backup/backup-latest/Lab/projects/2507blastx/output/250909_4474_samples/final_test/")

# read data
rate_df  <- read_tsv("infection_events_three_cols_genomad_positive_rate.tsv",  show_col_types = FALSE)
count_df <- read_tsv("total_count_of_infection_events.tsv", show_col_types = FALSE)
rate_cols <- c("infection_events", "infection_events_200aa", "infection_events_200aa_<95%")
genus_levels <- unique(rate_df$genus)

# long format
rate_long <- rate_df %>%
  select(genus, all_of(rate_cols)) %>%
  pivot_longer(cols = all_of(rate_cols), names_to = "event", values_to = "rate") %>%
  mutate(
    rate  = ifelse(is.nan(rate), NA_real_, as.numeric(rate)),
    genus = factor(genus, levels = rev(genus_levels)),
    event = factor(event, levels = rate_cols)
  )

count_long <- rate_df %>%
  distinct(genus) %>%
  left_join(count_df %>% select(genus, total_count), by = "genus") %>%
  mutate(
    total_count = ifelse(is.nan(total_count), NA_real_, as.numeric(total_count)),
    genus = factor(genus, levels = rev(genus_levels)),
    colname = "infection_events"
  )

# plot
rmax <- max(rate_long$rate, na.rm = TRUE)
cmax <- max(count_long$total_count, na.rm = TRUE)
pal_rate  <- c("white", RColorBrewer::brewer.pal(9, "BuGn")[-1])
pal_count <- c("white", RColorBrewer::brewer.pal(9, "YlGnBu")[-1])

base_theme <- theme_minimal(base_size = 11) +
  theme(
    text = element_text(color = "black"),
    axis.text = element_text(color = "black"),
    axis.title = element_blank(),
    panel.grid = element_blank(),

    panel.border = element_blank(),
    axis.line.x.bottom = element_line(color = "black", linewidth = 0.4),
    axis.line.y.left   = element_line(color = "black", linewidth = 0.4),
    axis.line.x.top    = element_blank(),
    axis.line.y.right  = element_blank(),

    plot.background  = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA)
  )

p_rate <- ggplot(rate_long, aes(x = event, y = genus, fill = rate)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(colours = pal_rate, limits = c(0, rmax), na.value = "white") +
  coord_fixed(ratio = 1) +
  base_theme +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

p_count <- ggplot(count_long, aes(x = colname, y = genus, fill = total_count)) +
  geom_tile(color = "white", linewidth = 0.3) +
  scale_fill_gradientn(colours = pal_count, limits = c(0, cmax), na.value = "white") +
  coord_fixed(ratio = 1) +
  base_theme +
  theme(
    axis.text.y  = element_blank(),
    axis.ticks.y = element_blank(),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1)
  )

p <- p_rate + p_count
p

ggsave("heatmap_genomad_positive_rate_plus_total_count_three_cols_genomadv1.11.2.pdf", p, width = 5, height = 5, units = "in")