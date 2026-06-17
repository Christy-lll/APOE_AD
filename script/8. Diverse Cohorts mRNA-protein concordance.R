# Set up ----
suppressMessages(library(biomaRt))
suppressMessages(library(DESeq2))
suppressMessages(library(here))
suppressMessages(library(limma))
suppressMessages(library(tidyverse))

pic_output_dir <- here("output", "DiverseCohorts", "pics")
table_output_dir <- here("output", "DiverseCohorts", "tables")

dc_rna <- readRDS(here("raw", "DiverseCohorts", "DC-transcriptomics.rds"))
dc_prot <- readRDS(here("raw", "DiverseCohorts", "DC-proteomics.rds"))


# Preprocessing ----
# Filter to overlapping dlPFC cohort
ids_overlap <- intersect(dc_rna$meta$dlPFC$individualID, dc_prot$meta$individualID)
message("dlPFC | Transcriptomic: ", nrow(dc_rna$meta$dlPFC),
        " | Proteomic: ", nrow(dc_prot$meta),
        " | Overlap: ", length(ids_overlap))

rna_meta <- dc_rna$meta$dlPFC %>% filter(individualID %in% ids_overlap) %>% arrange(individualID)
prot_meta <- dc_prot$meta %>% filter(individualID %in% ids_overlap) %>% arrange(individualID)
rna_count <- dc_rna$count$dlPFC[, rna_meta$specimenID, drop = FALSE]
prot_expr <- dc_prot$expr[, prot_meta$specimenID, drop = FALSE]

# all(rna_meta$individualID == prot_meta$individualID)
rm(dc_rna, dc_prot)


# Helper functions ----
p_val_threshold <- 0.05

run_deseq2 <- function(count, meta, contrast) {
  meta <- as.data.frame(meta) %>% droplevels()
  rownames(meta) <- meta$specimenID
  
  count <- round(count[, meta$specimenID])
  mode(count) <- "integer"
  
  site_term <- if (nlevels(meta$dataGenerationSite) > 1) "+ dataGenerationSite" else ""
  design_formula <- as.formula(paste0("~ apoe4 + ADoutcome + sex + age_death_scaled", site_term))
  message("design: ", deparse(design_formula))
  
  dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = design_formula)
  
  min_n <- min(table(meta[[contrast[1]]]))
  n_before <- nrow(dds)
  dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ]
  message("Before: ", n_before, " | After: ", nrow(dds), " | Removed: ", n_before - nrow(dds))
  
  dds <- DESeq(dds)
  
  results(dds, contrast = contrast, alpha = p_val_threshold) %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_gene_id") %>%
    mutate(ensembl_gene_id = sub("\\..*", "", ensembl_gene_id),
           negLog10FDR = -log10(padj),
           direction = case_when(padj < p_val_threshold & log2FoldChange > 0 ~ "Upregulated",
                                 padj < p_val_threshold & log2FoldChange < 0 ~ "Downregulated",
                                 TRUE ~ "Not significant")) %>%
    arrange(padj)
}

run_limma <- function(expr, meta, coef) {
  meta <- as.data.frame(meta)
  design_mat <- model.matrix(~ apoe4 + ADoutcome + sex + age_death_num + batch, data = meta)
  fit <- lmFit(expr, design_mat) %>% eBayes()
  
  topTable(fit, coef = coef, number = Inf) %>%
    rownames_to_column("protein_id") %>%
    mutate(negLog10FDR = -log10(adj.P.Val),
           direction = case_when(adj.P.Val < p_val_threshold & logFC > 0 ~ "Upregulated",
                                 adj.P.Val < p_val_threshold & logFC < 0 ~ "Downregulated",
                                 TRUE ~ "Not significant"))
}

rna_to_gene <- function(ensembl_ids) {
  tryCatch({
    mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
    getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
          filters = "ensembl_gene_id",
          values = ensembl_ids,
          mart = mart) %>%
      distinct(ensembl_gene_id, .keep_all = TRUE) %>%
      mutate(hgnc_symbol = if_else(hgnc_symbol == "", ensembl_gene_id, hgnc_symbol))
  }, error = function(e) {
    message("Ensembl unavailable, falling back to org.Hs.eg.db: ", conditionMessage(e))
    suppressMessages(library(org.Hs.eg.db))
    AnnotationDbi::select(org.Hs.eg.db,
                          keys = ensembl_ids,
                          columns = "SYMBOL",
                          keytype = "ENSEMBL") %>%
      dplyr::rename(ensembl_gene_id = ENSEMBL, hgnc_symbol = SYMBOL) %>%
      distinct(ensembl_gene_id, .keep_all = TRUE) %>%
      mutate(hgnc_symbol = if_else(is.na(hgnc_symbol), ensembl_gene_id, hgnc_symbol))
  })
}

dedup_by_symbol <- function(df, expr_col) {
  df %>%
    filter(!is.na(hgnc_symbol)) %>%
    group_by(hgnc_symbol) %>%
    slice_max(!!sym(expr_col), n = 1, with_ties = FALSE) %>%
    ungroup()
}


# Concordance of APOE ε4-Associated RNA and Protein Change ----
## DEGs via DESeq2 ----
rna_res_apoe4 <- run_deseq2(rna_count, rna_meta, contrast = c("apoe4", "APOE4pos", "APOE4neg"))

## DAPs via limma ----
prot_res_apoe4 <- run_limma(prot_expr, prot_meta, coef = "apoe4APOE ε4+") %>%
  mutate(hgnc_symbol = sub("\\|.*", "", protein_id),
         hgnc_symbol = if_else(hgnc_symbol == "", NA_character_, hgnc_symbol))

## Gene ID mapping ----
id_map_apoe4 <- rna_to_gene(rna_res_apoe4$ensembl_gene_id)
rna_res_apoe4 <- left_join(rna_res_apoe4, id_map_apoe4, by = "ensembl_gene_id")

## log2FC scatter plot ----
rna_clean_apoe4 <- dedup_by_symbol(rna_res_apoe4, "baseMean") %>%
  dplyr::select(hgnc_symbol, rna_lfc = log2FoldChange, rna_dir = direction)

prot_clean_apoe4 <- dedup_by_symbol(prot_res_apoe4, "AveExpr") %>%
  dplyr::select(hgnc_symbol, prot_lfc = logFC, prot_dir = direction)

joined_apoe4 <- inner_join(rna_clean_apoe4, prot_clean_apoe4, by = "hgnc_symbol") %>%
  mutate(de_status = case_when(rna_dir != "Not significant" & prot_dir != "Not significant" ~ "DEG + DAP",
                               rna_dir != "Not significant" ~ "DEG only",
                               prot_dir != "Not significant" ~ "DAP only",
                               TRUE ~ "Not significant"),
         label = if_else(de_status != "Not significant", hgnc_symbol, NA_character_))

test_apoe4 <- cor.test(joined_apoe4$rna_lfc, joined_apoe4$prot_lfc, method = "spearman")
stats_label_apoe4 <- paste0("Spearman r = ", round(test_apoe4$estimate, 3),
                            ", p = ", format.pval(test_apoe4$p.value, digits = 3),
                            "\nn = ", nrow(joined_apoe4), " genes")

p_scatter_apoe4 <- ggplot(joined_apoe4, aes(x = rna_lfc, y = prot_lfc)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.4) +
  geom_point(size = 0.8, color = "grey", alpha = 0.8) +
  geom_smooth(method = "lm", color = "steelblue", linewidth = 0.6, se = TRUE) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = stats_label_apoe4, size = 3.2, fontface = "italic") +
  labs(title = "dlPFC RNA vs Protein log2FC (APOE ε4+ vs APOE ε4-)",
       x = "RNA log2 fold change", y = "Protein log2 fold change", color = NULL) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(pic_output_dir, "apoe4_dlpfc_logFC_scatter.tiff"),
       p_scatter_apoe4, width = 7, height = 6, dpi = 300, compression = "lzw")


# Concordance of AD-Associated RNA and Protein Change ----
## DEGs via DESeq2 ----
rna_res_ad <- run_deseq2(rna_count, rna_meta, contrast = c("ADoutcome", "AD", "Control"))

## DAPs via limma ----
prot_res_ad <- run_limma(prot_expr, prot_meta, coef = "ADoutcomeAD") %>%
  mutate(hgnc_symbol = sub("\\|.*", "", protein_id),
         hgnc_symbol = if_else(hgnc_symbol == "", NA_character_, hgnc_symbol))

## Gene ID mapping ----
id_map_ad <- rna_to_gene(rna_res_ad$ensembl_gene_id)
rna_res_ad <- left_join(rna_res_ad, id_map_ad, by = "ensembl_gene_id")

## Concordance summary ----
rna_sig_ad <- filter(rna_res_ad, direction != "Not significant", !is.na(hgnc_symbol))
prot_sig_ad <- filter(prot_res_ad, direction != "Not significant", !is.na(hgnc_symbol))

concordance_ad <- inner_join(
  dplyr::select(rna_sig_ad, hgnc_symbol, rna_direction = direction,
                rna_lfc = log2FoldChange, rna_padj = padj),
  dplyr::select(prot_sig_ad, hgnc_symbol, prot_direction = direction,
                prot_lfc = logFC, prot_padj = adj.P.Val),
  by = "hgnc_symbol"
) %>%
  mutate(concordant = rna_direction == prot_direction)

message("dlPFC AD | DEGs: ", nrow(rna_sig_ad), " | DAPs: ", nrow(prot_sig_ad),
        " | Overlap: ", nrow(concordance_ad), " | Concordant: ", sum(concordance_ad$concordant))

## log2FC scatter plot ----
rna_clean_ad <- dedup_by_symbol(rna_res_ad, "baseMean") %>%
  dplyr::select(hgnc_symbol, rna_lfc = log2FoldChange, rna_dir = direction)

prot_clean_ad <- dedup_by_symbol(prot_res_ad, "AveExpr") %>%
  dplyr::select(hgnc_symbol, prot_lfc = logFC, prot_dir = direction)

joined_ad <- inner_join(rna_clean_ad, prot_clean_ad, by = "hgnc_symbol") %>%
  mutate(de_status = case_when(
    rna_dir != "Not significant" & prot_dir != "Not significant" &
      rna_dir == prot_dir ~ "Concordant (DEG + DAP)",
    rna_dir != "Not significant" & prot_dir != "Not significant" &
      rna_dir != prot_dir ~ "Discordant (DEG + DAP)",
    rna_dir != "Not significant" ~ "DEG only",
    prot_dir != "Not significant" ~ "DAP only",
    TRUE ~ "Not significant"
  ))

test_ad <- cor.test(joined_ad$rna_lfc, joined_ad$prot_lfc, method = "spearman")
stats_label_ad <- paste0("Spearman r = ", round(test_ad$estimate, 3),
                         ", p = ", format.pval(test_ad$p.value, digits = 3),
                         "\nn = ", nrow(joined_ad), " genes")

p_scatter_ad <- ggplot(joined_ad, aes(x = rna_lfc, y = prot_lfc)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.4) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.4) +
  geom_point(color = "grey", size = 0.8) +
  geom_smooth(method = "lm", color = "indianred", linewidth = 0.6, se = TRUE) +
  annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
           label = stats_label_ad, size = 3.2, fontface = "italic") +
  labs(title = "dlPFC RNA vs Protein log2FC (AD vs Control)",
       x = "RNA log2 fold change", y = "Protein log2 fold change", color = NULL) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(pic_output_dir, "ad_dlpfc_logFC_scatter.tiff"),
       p_scatter_ad, width = 7, height = 6, dpi = 300, compression = "lzw")
