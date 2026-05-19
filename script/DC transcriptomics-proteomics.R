# Set up ----
suppressMessages(library(biomaRt))
suppressMessages(library(DESeq2))
suppressMessages(library(ggvenn))
suppressMessages(library(here))
suppressMessages(library(limma))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))

pic_output_dir <- here("output", "DiverseCohorts", "pics")
table_output_dir <- here("output", "DiverseCohorts", "tables")

dc_rna <- readRDS(here("raw", "DiverseCohorts", "DC-transcriptomics.rds"))
dc_prot <- readRDS(here("raw", "DiverseCohorts", "DC-proteomics.rds"))

# Preprocessing ----
# Filter overlapping cohort for each region
get_overlap <- function(rna_count, rna_meta, prot_expr, prot_meta, region_label) {
  ids_overlap <- intersect(rna_meta$individualID, prot_meta$individualID)
  
  message(region_label,
          " | Transcriptomic: ", nrow(rna_meta),
          " | Proteomic: ", nrow(prot_meta),
          " | Overlap: ", length(ids_overlap))
  
  rna_meta_ol <- rna_meta %>% filter(individualID %in% ids_overlap) %>% arrange(individualID)
  prot_meta_ol <- prot_meta %>% filter(individualID %in% ids_overlap) %>% arrange(individualID)
  rna_count_ol <- rna_count[, rna_meta_ol$specimenID, drop = FALSE]
  prot_expr_ol <- prot_expr[, prot_meta_ol$specimenID, drop = FALSE]
  
  list(rna_count = rna_count_ol,
       rna_meta = rna_meta_ol,
       prot_expr = prot_expr_ol,
       prot_meta = prot_meta_ol)
}

regions <- map(c("dlPFC", "STG"), ~ get_overlap(dc_rna$count[[.x]], dc_rna$meta[[.x]],
                                                dc_prot[[.x]]$expr, dc_prot[[.x]]$meta, .x)) %>% 
  setNames(c("dlPFC", "STG"))

# all(regions$dlPFC$rna_meta$individualID == regions$dlPFC$prot_meta$individualID)
# all(regions$STG$rna_meta$individualID == regions$STG$prot_meta$individualID)
rm(dc_rna, dc_prot)


# Concordance of APOE4-Associated RNA and Protein Change ----
p_val_threshold <- 0.05

## APOE4 DEGs via deseq2 ----
run_deseq2 <- function(rna_count, rna_meta, region_label, contrast) {
  meta <- as.data.frame(rna_meta) %>% droplevels()
  rownames(meta) <- meta$specimenID
  
  count <- round(rna_count[, meta$specimenID])
  mode(count) <- "integer"
  
  site_term <- if (nlevels(meta$dataGenerationSite) > 1) "+ dataGenerationSite" else "" # adjust for data generation site if more than 1 is presented
  design_formula <- as.formula(paste0("~ apoe4 + ADoutcome + sex + age_death_scaled", site_term))
  message(region_label, " design: ", deparse(design_formula))
  
  dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = design_formula)
  
  min_n <- min(table(meta[[contrast[1]]]))
  n_before <- nrow(dds)
  dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ] # filter: keep genes with >= 10 counts in at least as many samples as the smallest contrast group
  message(region_label, " | Before: ", n_before, " | After: ", nrow(dds), " | Removed: ", n_before - nrow(dds))
  
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

rna_res <- imap(regions, ~ run_deseq2(.x$rna_count, .x$rna_meta, toupper(.y),
                                      contrast = c("apoe4", "APOE4pos", "APOE4neg")))

## APOE4 DEPs via limma ----
run_limma <- function(prot_expr, prot_meta, coef) {
  meta <- as.data.frame(prot_meta)
  design_mat <- model.matrix(~ apoe4 + ADoutcome + sex + age_death_num + batch, data = meta)
  fit <- lmFit(prot_expr, design_mat) %>% eBayes()
  
  topTable(fit, coef = coef, number = Inf) %>%
    rownames_to_column("protein_id") %>%
    mutate(negLog10FDR = -log10(adj.P.Val),
           direction = case_when(adj.P.Val < p_val_threshold & logFC > 0 ~ "Upregulated",
                                 adj.P.Val < p_val_threshold & logFC < 0 ~ "Downregulated",
                                 TRUE ~ "Not significant"))
}

prot_res <- map(regions, ~ run_limma(.x$prot_expr, .x$prot_meta, coef = "apoe4APOE4+"))

## Gene ID mapping ----
# RNA
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

id_map <- map(rna_res, "ensembl_gene_id") %>% unlist() %>% unique() %>% rna_to_gene()
rna_res <- map(rna_res, ~ left_join(.x, id_map, by = "ensembl_gene_id"))

# Proteins
protein_to_gene <- function(prot_res_list) {
  map(prot_res_list, ~ mutate(.x,
                              hgnc_symbol = sub("\\|.*", "", protein_id),
                              hgnc_symbol = if_else(hgnc_symbol == "", NA_character_, hgnc_symbol)))
}

prot_res <- protein_to_gene(prot_res)

## Concordance of DEGs and DEPs----
venn_plots <- imap(regions, ~ {
  rna_sig <- filter(rna_res[[.y]], direction != "Not significant", !is.na(hgnc_symbol))$hgnc_symbol
  prot_sig <- filter(prot_res[[.y]], direction != "Not significant", !is.na(hgnc_symbol))$hgnc_symbol
  
  ggvenn(list(RNA = rna_sig, Protein = prot_sig),
         fill_color = c("indianred", "steelblue"),
         stroke_size = 0.5, set_name_size = 4,
         text_size = 0, show_percentage = FALSE) +
    annotate("text", x = -0.8, y = 0,
             label = paste(rna_sig, collapse = "\n"),
             size = 2.2, color = "black") +
    annotate("text", x = 0.8, y = 0,
             label = paste(prot_sig, collapse = "\n"),
             size = 2.2, color = "black") +
    labs(title = .y) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
})

apoe_venn <- wrap_plots(venn_plots, nrow = 1)

# ggsave(file.path(pic_output_dir, "apoe4_DEGs_DEPs_concordance.png"), apoe_venn, width = 8, height = 6)

## Correlation of log FC -----
dedup_by_symbol <- function(df, expr_col) {
  df %>%
    filter(!is.na(hgnc_symbol)) %>%
    group_by(hgnc_symbol) %>%
    slice_max(!!sym(expr_col), n = 1, with_ties = FALSE) %>% # clean duplication in gene symbol during mapping
    ungroup()
}

make_lfc_scatter_apoe4 <- function(rna, prot, region_label) {
  rna_clean <- dedup_by_symbol(rna, "baseMean") %>%
    dplyr::select(hgnc_symbol, rna_lfc = log2FoldChange, rna_dir = direction)
  prot_clean <- dedup_by_symbol(prot, "AveExpr") %>%
    dplyr::select(hgnc_symbol, prot_lfc = logFC, prot_dir = direction)
  
  joined <- inner_join(rna_clean, prot_clean, by = "hgnc_symbol") %>%
    mutate(de_status = case_when(rna_dir != "Not significant" & prot_dir != "Not significant" ~ "DEG + DEP",
                                 rna_dir != "Not significant" ~ "DEG only",
                                 prot_dir != "Not significant" ~ "DEP only",
                                 TRUE ~ "Not significant"),
           label = if_else(de_status != "Not significant", hgnc_symbol, NA_character_))

  test <- cor.test(joined$rna_lfc, joined$prot_lfc, method = "spearman")
  stats_label <- paste0("Spearman r = ", round(test$estimate, 3),
                        ", p = ", round(test$p.value, 4),
                        "\nn = ", nrow(joined), " genes")
  
  ggplot(joined, aes(x = rna_lfc, y = prot_lfc)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.4) +
    geom_point(aes(color = de_status), size = 0.8) +
    geom_smooth(method = "lm", color = "black", linewidth = 0.6, se = TRUE) +
    ggrepel::geom_text_repel(aes(label = label), na.rm = TRUE, size = 3,
                             color = "black", max.overlaps = Inf) +
    scale_color_manual(values = c("DEG" = "steelblue", "DEP" = "indianred",
                                  "Not significant" = "grey")) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
             label = stats_label, size = 3.2, fontface = "italic") +
    labs(title = paste(region_label, "RNA vs Protein log2FC (APOE4+ vs APOE4-)"),
         x = "RNA log2 fold change", y = "Protein log2 fold change", color = NULL) +
    theme_classic() +
    theme(legend.position = "bottom")
}

scatter_plots <- imap(regions, ~ make_lfc_scatter_apoe4(rna_res[[.y]], prot_res[[.y]], .y))

# ggsave(file.path(pic_output_dir, "apoe4_dlpfc_logFC_scatter.png"), scatter_plots[["dlPFC"]], width = 8, height = 6)
# ggsave(file.path(pic_output_dir, "apoe4_stg_logFC_scatter.png"),  scatter_plots[["STG"]],  width = 8, height = 6)


# Concordance of AD-Associated RNA and Protein Change ----
## Some functions are defined in previous analysis

## AD DEGs via deseq2 ----
rna_res_ad <- imap(regions, ~ run_deseq2(.x$rna_count, .x$rna_meta, toupper(.y),
                                         contrast = c("ADoutcome", "AD", "Control")))

## AD DEPs via limma ----
prot_res_ad <- map(regions, ~ run_limma(.x$prot_expr, .x$prot_meta, coef = "ADoutcomeAD"))

## Gene ID mapping ----
id_map_ad <- map(rna_res_ad, "ensembl_gene_id") %>% unlist() %>% unique() %>% rna_to_gene()
rna_res_ad <- map(rna_res_ad, ~ left_join(.x, id_map_ad, by = "ensembl_gene_id"))

prot_res_ad <- protein_to_gene(prot_res_ad)

## Concordance of DEGs and DEPs----
get_concordance <- function(rna, prot, region_label) {
  rna_sig <- filter(rna, direction != "Not significant", !is.na(hgnc_symbol))
  prot_sig <- filter(prot, direction != "Not significant", !is.na(hgnc_symbol))
  
  overlap <- inner_join(
    dplyr::select(rna_sig, hgnc_symbol, rna_direction = direction,
                  rna_lfc = log2FoldChange, rna_padj = padj),
    dplyr::select(prot_sig, hgnc_symbol, prot_direction = direction,
                  prot_lfc = logFC, prot_padj = adj.P.Val),
    by = "hgnc_symbol"
  ) %>%
    mutate(concordant = rna_direction == prot_direction)
  
  message(region_label, " | DEGs: ", nrow(rna_sig), " | DEPs: ", nrow(prot_sig),
          " | Overlap: ", nrow(overlap), " | Concordant: ", sum(overlap$concordant))
  overlap
}

concordance_ad <- imap(regions, ~ get_concordance(rna_res_ad[[.y]], prot_res_ad[[.y]], .y))

## Correlation of log FC ----
make_lfc_scatter_ad <- function(rna, prot, region_label) {
  rna_clean <- dedup_by_symbol(rna, "baseMean") %>%
    dplyr::select(hgnc_symbol, rna_lfc = log2FoldChange, rna_dir = direction)
  prot_clean <- dedup_by_symbol(prot, "AveExpr") %>%
    dplyr::select(hgnc_symbol, prot_lfc = logFC, prot_dir = direction)
  
  joined <- inner_join(rna_clean, prot_clean, by = "hgnc_symbol") %>%
    mutate(de_status = case_when(
      rna_dir != "Not significant" & prot_dir != "Not significant" &
        rna_dir == prot_dir ~ "Concordant (DEG + DEP)",
      rna_dir != "Not significant" & prot_dir != "Not significant" &
        rna_dir != prot_dir ~ "Discordant (DEG + DEP)",
      rna_dir != "Not significant" ~ "DEG only",
      prot_dir != "Not significant" ~ "DEP only",
      TRUE ~ "Not significant"
    ))
  
  test <- cor.test(joined$rna_lfc, joined$prot_lfc, method = "spearman")
  stats_label <- paste0("Spearman r = ", round(test$estimate, 3),
                        ", p = ", round(test$p.value, 4),
                        "\nn = ", nrow(joined), " genes")
  
  ggplot(joined, aes(x = rna_lfc, y = prot_lfc)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.4) +
    geom_point(aes(color = de_status), size = 0.8) +
    geom_smooth(method = "lm", color = "black", linewidth = 0.6, se = TRUE) +
    scale_color_manual(values = c(
      "Concordant (DEG + DEP)" = "#FDAE61",
      "Discordant (DEG + DEP)" = "#ABDDA4",
      "DEG only" = "steelblue",
      "DEP only" = "indianred",
      "Not significant" = "grey"
    )) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
             label = stats_label, size = 3.2, fontface = "italic") +
    labs(title = paste(region_label, "RNA vs Protein log2FC (AD vs Control)"),
         x = "RNA log2 fold change", y = "Protein log2 fold change", color = NULL) +
    theme_classic() +
    theme(legend.position = "bottom")
}

scatter_plots_ad <- imap(regions, ~ make_lfc_scatter_ad(rna_res_ad[[.y]], prot_res_ad[[.y]], .y))

# ggsave(file.path(pic_output_dir, "ad_dlpfc_logFC_scatter.png"), scatter_plots_ad[["dlPFC"]], width = 8, height = 6)
# ggsave(file.path(pic_output_dir, "ad_stg_logFC_scatter.png"),  scatter_plots_ad[["STG"]],  width = 8, height = 6)
