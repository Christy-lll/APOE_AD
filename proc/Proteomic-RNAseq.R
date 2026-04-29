# Set up ----
suppressMessages(library(tidyverse))
library(here)

output_dir <- here("output", "Proteomic-RNAseq")

# iA results
iA_limma <- read.csv(here("output", "tables", "ROSMAP", "iAstrocytes_limma.csv"), header = T)
iA_mi <- read.csv(here("output", "tables", "ROSMAP", "iAstrocytes_MI.csv"), header = T)

# iN results
iN_limma <- read.csv(here("output", "tables", "ROSMAP", "iNeurons_limma.csv"), header = T)
iN_mi <- read.csv(here("output", "tables", "ROSMAP", "iNeurons_MI.csv"), header = T)

# dlPFC results 
dlPFC_limma <- read.csv(here("output", "tables", "DiverseCohorts", "DC_frontal_all_limma.csv"), header = T)

# RNA-seq results
RNA_dlPFC_deseq2 <- read.csv(here("output", "tables", "ROSMAP", "dlPFC_deseq2.csv"), header = T)
RNA_hCN_deseq2 <- read.csv(here("output", "tables", "ROSMAP", "hCN_deseq2.csv"), header = T)
RNA_PCC_deseq2 <- read.csv(here("output", "tables", "ROSMAP", "PCC_deseq2.csv"), header = T)

# Preprocessing ----
## Clean and deduplicate limma/mi results ----
iA_limma <- iA_limma %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  arrange(P.Value) %>%
  distinct(GeneSymbol, .keep_all = TRUE) %>%
  select(GeneSymbol, logFC, P.Value, adj.P.Val) 

iA_mi <- iA_mi %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  arrange(desc(importance)) %>%
  distinct(GeneSymbol, .keep_all = TRUE) %>%
  select(GeneSymbol, importance)

iN_limma <- iN_limma %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  arrange(P.Value) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  select(GeneSymbol = SYMBOL, logFC, P.Value, adj.P.Val)

iN_mi <- iN_mi %>%
  filter(!is.na(SYMBOL), SYMBOL != "") %>%
  arrange(desc(importance)) %>%
  distinct(SYMBOL, .keep_all = TRUE) %>%
  select(GeneSymbol = SYMBOL, importance)

dlPFC_limma <- dlPFC_limma %>%
  mutate(GeneSymbol = str_extract(Protein, "^[^|]+")) %>%
  filter(!is.na(GeneSymbol), GeneSymbol != "") %>%
  arrange(P.Value) %>%
  distinct(GeneSymbol, .keep_all = TRUE) %>%
  select(GeneSymbol, logFC, P.Value, adj.P.Val)

## Clean and deduplicate deseq2 results ----
clean_deseq2 <- function(df) {
  df %>%
    filter(!is.na(hgnc_symbol), hgnc_symbol != "") %>%
    arrange(pvalue) %>%
    distinct(hgnc_symbol, .keep_all = TRUE) %>%
    select(hgnc_symbol, log2FoldChange, pvalue, padj)
}

RNA_dlPFC_deseq2 <- clean_deseq2(RNA_dlPFC_deseq2)
RNA_hCN_deseq2 <- clean_deseq2(RNA_hCN_deseq2)
RNA_PCC_deseq2 <- clean_deseq2(RNA_PCC_deseq2)


# LIMMA & DESEQ2 result comparison ----
## Intersecting gene sets (limma & deseq2) ----
intersect_omics <- function(proteomics, transcriptomics) {
  inner_join(proteomics, transcriptomics,
             by = join_by(GeneSymbol == hgnc_symbol)) %>%
    select(gene = GeneSymbol,
           logFC_protein = logFC,
           logFC_rna = log2FoldChange,
           pval_protein = P.Value,
           padj_protein = adj.P.Val,
           pval_rna = pvalue,
           padj_rna = padj) %>%
    filter(!is.na(pval_protein), !is.na(pval_rna)) %>%
    mutate(rank_protein = sign(logFC_protein) * -log10(pval_protein),
           rank_rna = sign(logFC_rna) * -log10(pval_rna))
}

iA_dlPFC <- intersect_omics(iA_limma, RNA_dlPFC_deseq2)
iA_hCN <- intersect_omics(iA_limma, RNA_hCN_deseq2)
iA_PCC <- intersect_omics(iA_limma, RNA_PCC_deseq2)

iN_dlPFC <- intersect_omics(iN_limma, RNA_dlPFC_deseq2)
iN_hCN <- intersect_omics(iN_limma, RNA_hCN_deseq2)
iN_PCC <- intersect_omics(iN_limma, RNA_PCC_deseq2)

dlPFC_dlPFC <- intersect_omics(dlPFC_limma, RNA_dlPFC_deseq2)

## 2 way Venn diagram: DEGs vs DEPs (unadjusted p < 0.05) ----
library(ggVennDiagram)

plot_venn <- function(intersected_df, cell_type, region_name, output_dir) {
  dep <- intersected_df %>% filter(pval_protein < 0.05) %>% pull(gene)
  deg <- intersected_df %>% filter(pval_rna < 0.05) %>% pull(gene)
  
  venn_list <- list(DEP = dep, DEG = deg)
  overlap <- intersect(dep, deg)
  
  overlap_labels <- intersected_df %>%
    filter(gene %in% overlap) %>%
    mutate(direction = case_when(logFC_rna > 0 & logFC_protein > 0 ~ "↑↑",
                                 logFC_rna > 0 & logFC_protein < 0 ~ "↑↓",
                                 logFC_rna < 0 & logFC_protein < 0 ~ "↓↓",
                                 logFC_rna < 0 & logFC_protein > 0 ~ "↓↑"),
           label = paste0(gene, " ", direction)) %>%
    pull(label)
  
  overlap_label <- paste(overlap_labels, collapse = "\n")
  
  p <- ggVennDiagram(venn_list) +
    scale_fill_gradient(low = "white", high = "steelblue") +
    labs(title = paste(region_name, "Transcriptomics vs ", cell_type, "Proteomics" ),
         subtitle = "Threshold: unadjusted p < 0.05 | Direction: LogFC") +
    theme(plot.title = element_text(hjust = 0)) +
    annotate("text", x = 1.5, y = 1.9,
             label = overlap_label, size = 3)
  
  ggsave(file.path(output_dir, paste0("Venn_", cell_type, "_vs_", region_name, ".jpg")),
         plot = p, width = 9, height = 6, dpi = 300)
  
  cat(cell_type, "vs", region_name,
      "| Universe:", nrow(intersected_df),
      "| DEP:", length(dep),
      "| DEG:", length(deg),
      "| Overlap:", length(overlap), "\n")
  
  return(invisible(overlap_labels))
}

iA_venn_dlPFC <- plot_venn(iA_dlPFC, "iA", "dlPFC", output_dir)
iA_venn_hCN <- plot_venn(iA_hCN, "iA", "hCN", output_dir)
iA_venn_PCC <- plot_venn(iA_PCC, "iA", "PCC", output_dir)

iN_venn_dlPFC <- plot_venn(iN_dlPFC, "iN", "dlPFC", output_dir)
iN_venn_hCN <- plot_venn(iN_hCN, "iN", "hCN", output_dir)
iN_venn_PCC <- plot_venn(iN_PCC, "iN", "PCC", output_dir)

dlPFC_venn_dlPFC <- plot_venn(dlPFC_dlPFC, "dlPFC", "dlPFC", output_dir)

## 6 way Upset plot ----
upset_list <- list(iA_DEP = iA_limma %>% filter(P.Value < 0.05) %>% pull(GeneSymbol),
                   iN_DEP = iN_limma %>% filter(P.Value < 0.05) %>% pull(GeneSymbol),
                   dlPFC_DEP = dlPFC_limma %>% filter(P.Value < 0.05) %>% pull(GeneSymbol),
                   dlPFC_DEG = RNA_dlPFC_deseq2 %>% filter(pvalue < 0.05) %>% pull(hgnc_symbol),
                   hCN_DEG = RNA_hCN_deseq2 %>% filter(pvalue < 0.05) %>% pull(hgnc_symbol),
                   PCC_DEG = RNA_PCC_deseq2 %>% filter(pvalue < 0.05) %>% pull(hgnc_symbol))

png(file.path(output_dir, "UpSet_6way.png"),
    width = 2000, height = 1600, res = 300)

library(UpSetR)
upset(fromList(upset_list),
      nsets = 6,
      order.by = "freq",     
      decreasing = TRUE,
      mainbar.y.label = "Intersection size",
      sets.x.label = "Set size")

dev.off()

## Spearman rank correlation & Pearson correlation ----
run_correlation <- function(df, cell_type, region_name) {
  spearman <- cor.test(df$rank_protein, df$rank_rna, method = "spearman")
  pearson  <- cor.test(df$logFC_protein, df$logFC_rna, method = "pearson")
  
  tibble(cell_type = cell_type,
         region = region_name,
         spearman_rho = spearman$estimate,
         spearman_pval = spearman$p.value,
         pearson_r = pearson$estimate,
         pearson_pval = pearson$p.value,
         n = nrow(df))
}

correlation_results <- bind_rows(run_correlation(dlPFC_dlPFC, "dlPFC", "dlPFC"), # 8192 overlap
                                 run_correlation(iA_dlPFC, "iA", "dlPFC"), # 7161 overlap
                                 run_correlation(iA_hCN, "iA", "hCN"), # 7161 overlap
                                 run_correlation(iA_PCC, "iA", "PCC"), # 7161 overlap
                                 run_correlation(iN_dlPFC, "iN", "dlPFC"), # 1144 overlap
                                 run_correlation(iN_hCN, "iN", "hCN"), # 1144 overlap
                                 run_correlation(iN_PCC, "iN", "PCC")) # 1145 overlap

write.csv(correlation_results, file.path(output_dir, "PROT_RNA_Correlation.csv"), row.names = FALSE)


## Rank Rank Hypergeometric Overlap (RRHO) analysis ----
library(RRHO2)

run_rrho <- function(df, cell_type, region_name, output_dir) {
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  list1 <- df %>% select(gene, rank_protein) %>% as.data.frame()
  list2 <- df %>% select(gene, rank_rna) %>% as.data.frame()
  
  result <- RRHO2_initialize(list1, list2,
                             stepsize = 30,
                             labels = c(paste(cell_type, "Proteomics"), paste("RNA", region_name)),
                             log10.ind = TRUE)
  
  jpeg(file.path(output_dir, paste0("RRHO2_heatmap_", cell_type, "_vs_", region_name, ".jpg")),
       width = 800, height = 800)
  RRHO2_heatmap(result)
  dev.off()
  
  jpeg(file.path(output_dir, paste0("RRHO2_uu_venn_", cell_type, "_vs_", region_name, ".jpg")),
       width = 800, height = 400)
  RRHO2_vennDiagram(result, type = "uu")
  dev.off()
  
  jpeg(file.path(output_dir, paste0("RRHO2_dd_venn_", cell_type, "_vs_", region_name, ".jpg")),
       width = 800, height = 400)
  RRHO2_vennDiagram(result, type = "dd")
  dev.off()
  
  return(result)
}

# iA RRHO
iA_rrho_dlPFC <- run_rrho(iA_dlPFC, "iA", "dlPFC", output_dir)
iA_rrho_hCN <- run_rrho(iA_hCN,   "iA", "hCN",   output_dir)
iA_rrho_PCC <- run_rrho(iA_PCC,   "iA", "PCC",   output_dir)

# iN RRHO
iN_rrho_dlPFC <- run_rrho(iN_dlPFC, "iN", "dlPFC", output_dir)
iN_rrho_hCN <- run_rrho(iN_hCN,   "iN", "hCN",   output_dir)
iN_rrho_PCC <- run_rrho(iN_PCC,   "iN", "PCC",   output_dir)

# dlPFC RRHO
dlPFC_rrho_dlPFC <- run_rrho(dlPFC_dlPFC, "dlPFC", "dlPFC", output_dir)


# MI and DESEQ2 result comparison ----
## Intersecting gene sets (MI & deseq2) ---- (no direction)
intersect_omics_MI <- function(proteomics, transcriptomics) {
  inner_join(proteomics, transcriptomics,
             by = join_by(GeneSymbol == hgnc_symbol)) %>%
    select(gene = GeneSymbol,
           MI_protein = importance,
           logFC_rna = log2FoldChange,
           pval_rna = pvalue,
           padj_rna = padj) %>%
    filter(pval_rna < 0.05)
}

iA_mi_dlPFC <- intersect_omics_MI(iA_mi, RNA_dlPFC_deseq2) # LSP1
iA_mi_hCN <- intersect_omics_MI(iA_mi, RNA_hCN_deseq2) # 0 overlap
iA_mi_PCC <- intersect_omics_MI(iA_mi, RNA_PCC_deseq2) # LYRM1

iN_mi_dlPFC <- intersect_omics_MI(iN_mi, RNA_dlPFC_deseq2) # 0 overlap
iN_mi_hCN <- intersect_omics_MI(iN_mi, RNA_hCN_deseq2) # BIN1
iN_mi_PCC <- intersect_omics_MI(iN_mi, RNA_PCC_deseq2) # PDIA3



upset_list_2 <- list(iA_non_zero_MI_ = iA_mi %>% pull(GeneSymbol),
                     iN_non_zero_MI = iN_mi  %>% pull(GeneSymbol),
                     dlPFC_DEG = RNA_dlPFC_deseq2 %>% filter(pvalue < 0.05) %>% pull(hgnc_symbol),
                     hCN_DEG = RNA_hCN_deseq2 %>% filter(pvalue < 0.05) %>% pull(hgnc_symbol),
                     PCC_DEG = RNA_PCC_deseq2 %>% filter(pvalue < 0.05) %>% pull(hgnc_symbol))

png(file.path(output_dir, "UpSet_5way_withMI.png"),
    width = 2000, height = 1600, res = 300)

UpSetR::upset(fromList(upset_list_2),
              nsets = 5,
              order.by = "freq",     
              decreasing = TRUE,
              mainbar.y.label = "Intersection size",
              sets.x.label = "Set size")

dev.off()


