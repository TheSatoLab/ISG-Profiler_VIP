rm( list=ls(all=TRUE) )
library("tidyverse")
library("DESeq2")
packageVersion("DESeq2")

# output
setwd("/Users/kyokokurihara/iLab/itolab_backup/backup-latest/Lab/projects/2507blastx/data/251215_rna_seq/To_kurihara_downloaded260209_chaphama")

# input
sample_directoy <- "/Users/kyokokurihara/iLab/itolab_backup/backup-latest/Lab/projects/2507blastx/output/250909_4474_samples/group_analysis/Chaphamaparvovirus_bioproject_260206_blastx_filtered/"

for (i in c("PRJNA577590", "PRJNA612882","PRJNA622813")) {
  count <- read.table(paste("featureCount/all_counts_",i,"_ver2.txt", sep = ""),
                      sep = "\t",  header = T, row.names = 1)
  
  sample <- read.csv(paste(sample_directoy, "SraRunTable_", i, ".csv", sep = ""),
                     header = TRUE, stringsAsFactors = FALSE)
  
  sample <- sample %>%
    filter((!is.na(Positive) & Positive == 1) | (!is.na(Negative) & Negative == 1))
  
  sample <- sample %>%
    mutate(con = ifelse(!is.na(Positive) & Positive == 1, "Positive",
                        ifelse(!is.na(Negative) & Negative == 1, "Negative", NA))) %>%
    filter(!is.na(con))
  
  sample <- sample %>% arrange(con)
  
  out <- sample %>%
    group_by(con) %>%
    summarise(n = n(), .groups = "drop") %>%
    mutate(ratio = 100 * n / sum(n))
  print(i)
  print(out)
  
  ### Positive vs Negative
  
  count <- count %>% select(one_of(sample$Run))
  count <- as.matrix(count)
  
  group <- data.frame(con = factor(sample$con, levels = c("Negative", "Positive")))
  rownames(group) <- sample$Run
  
  dds <- DESeqDataSetFromMatrix(countData = count, colData = group, design = ~ con)
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("con", "Positive", "Negative"))
  head(res)
  
  out.pdfname <- paste("output/plotMA_", i, "_Positive_vs_Negative.pdf", sep = "")
  pdf(out.pdfname, width = 5, height = 5)
  plotMA(res, alpha = 0.01)
  dev.off()
  
  f.name <- paste("output/result_",i, "_Positive_vs_Negative.txt", sep = "")
  write.table(res, file = f.name, row.names = T, col.names = T, sep = "\t", quote=F)
  
  ortholog <- read.table("go/gene_orthologs_hs_gallus_table.txt", sep = "\t",  header = T)
  
  res <- res %>% as.data.frame() %>% mutate(Symbol.2 = rownames(res))
  res.hs.Id <- res %>% left_join(ortholog, by = "Symbol.2")
  res.hs.Id <- res.hs.Id %>% select(Symbol.1, Symbol.2, 1:6)
  res.hs.Id <- res.hs.Id %>% dplyr::rename(gene_name_hs = Symbol.1, gene_name_gallus = Symbol.2)
  
  f.name <- paste("output/result_",i, "_Positive_vs_Negative_hs_gallus_and_Gene.txt", sep = "")
  write.table(res.hs.Id, file = f.name, row.names = F, col.names = T, sep = "\t", quote=F)
}