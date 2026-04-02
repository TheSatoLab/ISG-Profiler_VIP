library(ComplexHeatmap)
library(circlize)
grid::unit

# change working directory
setwd("/Users/kyokokurihara/iLab/itolab_backup/backup-latest/Lab/projects/2507blastx/output/250909_4474_samples/final_test/chaphama_260206_MSigDB_Hallmark_2020_stat_blastx_filtered/heatmap_chaphama/")

# read tables
mat_df  <- read.csv("gene_project_log2FC_matrix.csv", row.names = 1, check.names = FALSE)
term_bin_df <- read.csv("term_anno.csv", row.names = 1, check.names = FALSE)

# format tables
mat <- as.matrix(mat_df)
term_bin <- as.matrix(term_bin_df[rownames(mat), , drop=FALSE])

# plot term
term_col <- c("0"="white", "1"="#2B2D42")
ht_term <- Heatmap(
  term_bin,
  name = "term",
  col = term_col,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = FALSE,
  show_column_names = TRUE,
  column_names_rot = 90,
  width = unit(ncol(term_bin) * 4, "mm"),
  show_heatmap_legend = FALSE
)

# plot Log2FC
col_fun <- colorRamp2(c(-6, 0, 6), c("#1B9E9E", "white", "#C51B7D"))
ht_main <- Heatmap(
  mat,
  name = "log2FC",
  col = col_fun,
  # left_annotation = left_anno,
  cluster_rows = FALSE,
  cluster_columns = FALSE,
  show_row_names = TRUE,
  row_names_side = "right",
  column_names_rot = 90,
)

# draw and save
ht <- ht_term + ht_main
png("stat_log2FC_heatmap.png", width = 4, height = 11.2, units = "in", res = 300)
draw(ht)
dev.off()

pdf("stat_log2FC_heatmap.pdf", width = 4, height = 11.2)
draw(ht)
dev.off()