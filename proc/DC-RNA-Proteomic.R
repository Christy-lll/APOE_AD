# Set up ----
suppressMessages(library(tidyverse))
library(here)

pic_output_dir <- here("output", "pics", "DiverseCohorts")
table_output_dir <- here("output", "tables", "DiverseCohorts")

# Transcriptomics data
rna_count_by_region <- readRDS(here("raw", "DiverseCohorts", "rna_count_by_region.rds"))
rna_meta_by_region <- readRDS(here("raw", "DiverseCohorts", "rna_meta_by_region.rds"))

dlpfc_rna_count <- rna_count_by_region[["dorsolateral_prefrontal_cortex"]]
dlpfc_rna_meta <- rna_meta_by_region[["dorsolateral_prefrontal_cortex"]]
stg_rna_count <- rna_count_by_region[["superior_temporal_gyrus"]]
stg_rna_meta <- rna_meta_by_region[["superior_temporal_gyrus"]]

# Proteomics data
dlpfc_prot_expr <- readRDS(here("raw", "DiverseCohorts", "dlpfc_prot_expr.rds"))
dlpfc_prot_meta <- readRDS(here("raw", "DiverseCohorts", "dlpfc_prot_meta.rds"))
stg_prot_expr <- readRDS(here("raw", "DiverseCohorts", "stg_prot_expr.rds"))
stg_prot_meta <- readRDS(here("raw", "DiverseCohorts", "stg_prot_meta.rds"))


# Preprocessing: Filter overlapping cohort for each region ----
get_overlap <- function(rna_count, rna_meta, prot_expr, prot_meta, region_label) {
  ids_overlap <- intersect(rna_meta$individualID, prot_meta$individualID)
  
  message(region_label, 
          " | Transcriptomic: ", nrow(rna_meta),
          " | Proteomic: ", nrow(prot_meta),
          " | Overlap: ", length(ids_overlap))
  
  rna_meta_ol <- filter(rna_meta, individualID %in% ids_overlap)
  prot_meta_ol <- filter(prot_meta, individualID %in% ids_overlap)
  rna_count_ol <- rna_count[, rna_meta_ol$specimenID, drop = FALSE]
  prot_expr_ol <- prot_expr[, prot_meta_ol$specimenID, drop = FALSE]
  
  list(rna_count = rna_count_ol,
       rna_meta  = rna_meta_ol,
       prot_expr = prot_expr_ol,
       prot_meta = prot_meta_ol)
}

dlpfc <- get_overlap(dlpfc_rna_count, dlpfc_rna_meta, dlpfc_prot_expr, dlpfc_prot_meta, "DLPFC")
stg <- get_overlap(stg_rna_count, stg_rna_meta, stg_prot_expr, stg_prot_meta, "STG")

rm(dlpfc_rna_count, dlpfc_rna_meta, dlpfc_prot_expr, dlpfc_prot_meta,
   stg_rna_count, stg_rna_meta, stg_prot_expr, stg_prot_meta,
   rna_count_by_region, rna_meta_by_region)


# Analysis 1: Concordance of DESeq2 & limma results ----
suppressWarnings(library(DESeq2))
library(limma)

p_val_threshold <- 0.05
regions <- list(dlpfc = dlpfc, stg = stg)

## DESeq2 (RNA) ----
run_deseq2 <- function(rna_count, rna_meta, region_label) {
  meta <- as.data.frame(rna_meta) %>% droplevels()
  rownames(meta) <- meta$specimenID
  meta$age_death_scaled <- scale(meta$age_death_num)[, 1] # scale age
  
  count <- round(rna_count[, meta$specimenID]) # round count to nearest integer
  mode(count) <- "integer"
  
  site_term <- if (nlevels(meta$dataGenerationSite) > 1) "+ dataGenerationSite" else ""
  design_formula <- as.formula(paste("~ apoe4 + ADoutcome + sex + age_death_scaled", site_term))
  message(region_label, " design: ", deparse(design_formula)) # STG only contains sample from Mayo
  
  dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = design_formula)
  
  min_n <- min(table(meta$apoe4))
  n_before <- nrow(dds)
  dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ]
  message(region_label, " DESeq2 | Before: ", n_before, " | After: ", nrow(dds), " | Removed: ", n_before - nrow(dds))
  
  dds <- DESeq(dds)
  
  results(dds, contrast = c("apoe4", "APOE4pos", "APOE4neg"), alpha = p_val_threshold) %>%
    as.data.frame() %>%
    rownames_to_column("ensembl_gene_id") %>%
    mutate(ensembl_gene_id = sub("\\..*", "", ensembl_gene_id),
           negLog10FDR = -log10(padj),
           direction = case_when(padj < p_val_threshold & log2FoldChange > 0 ~ "Upregulated",
                                 padj < p_val_threshold & log2FoldChange < 0 ~ "Downregulated",
                                 TRUE ~ "Not significant")) %>%
    arrange(padj)
}

rna_res <- imap(regions, ~ run_deseq2(.x$rna_count, .x$rna_meta, toupper(.y)))

## Limma (Protein) ----
run_limma <- function(prot_expr, prot_meta) {
  meta <- as.data.frame(prot_meta)
  design_mat <- model.matrix(~ apoe4 + ADoutcome + sex + age_death_num + batch, data = meta)
  fit <- lmFit(prot_expr, design_mat) %>% eBayes()
  
  topTable(fit, coef = "apoe4APOE4+", number = Inf) %>%
    rownames_to_column("protein_id") %>%
    mutate(negLog10FDR = -log10(adj.P.Val),
           direction = case_when(adj.P.Val < p_val_threshold & logFC > 0 ~ "Upregulated",
                                 adj.P.Val < p_val_threshold & logFC < 0 ~ "Downregulated",
                                 TRUE ~ "Not significant"))
}

prot_res <- map(regions, ~ run_limma(.x$prot_expr, .x$prot_meta))

## Map Ensembl/Uniprot IDs to HGNC symbols ----
library(biomaRt)
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

all_ensembl_ids <- map(rna_res, "ensembl_gene_id") %>% unlist() %>% unique()
id_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = "ensembl_gene_id",
                values = all_ensembl_ids,
                mart = mart) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  mutate(hgnc_symbol = if_else(hgnc_symbol == "", ensembl_gene_id, hgnc_symbol))

rna_res <- map(rna_res, ~ left_join(.x, id_map, by = "ensembl_gene_id"))

prot_res <- map(prot_res, ~ mutate(.x, 
                                   hgnc_symbol = sub("\\|.*", "", protein_id), 
                                   hgnc_symbol = if_else(hgnc_symbol == "", NA_character_, hgnc_symbol)))

## Save results ----
# iwalk(rna_res,  ~ write.csv(.x, file.path(table_output_dir, paste0("DC_multiomic_", .y, "_deseq2.csv")), row.names = FALSE))
# iwalk(prot_res, ~ write.csv(.x, file.path(table_output_dir, paste0("DC_multiomic_", .y, "_limma.csv")),  row.names = FALSE))

## Concordance ----
get_concordance <- function(rna, prot, region_label) {
  rna_sig <- filter(rna,  direction != "Not significant", !is.na(hgnc_symbol))
  prot_sig <- filter(prot, direction != "Not significant", !is.na(hgnc_symbol))
  
  overlap <- inner_join(
    select(rna_sig, hgnc_symbol, rna_direction= direction, rna_lfc = log2FoldChange, rna_padj = padj),
    select(prot_sig, hgnc_symbol, prot_direction = direction, prot_lfc = logFC, prot_padj = adj.P.Val),
    by = "hgnc_symbol"
  ) %>%
    mutate(concordant = rna_direction == prot_direction)
  
  message(region_label, " | DEGs: ", nrow(rna_sig), " | DEPs: ", nrow(prot_sig),
          " | Overlap: ", nrow(overlap), " | Concordant: ", sum(overlap$concordant))
  overlap
}

concordance <- imap(regions, ~ get_concordance(rna_res[[.y]], prot_res[[.y]], toupper(.y))) 
# No overlap DEP and DEG in either region

## Visualisation ----
library(ggvenn)
venn_plots <- imap(venn_data, ~ {
  rna_sig  <- filter(rna_res[[.y]],  direction != "Not significant", !is.na(hgnc_symbol))$hgnc_symbol
  prot_sig <- filter(prot_res[[.y]], direction != "Not significant", !is.na(hgnc_symbol))$hgnc_symbol
  
  ggvenn(.x, fill_color = c("mistyrose", "steelblue"),
         stroke_size = 0.5, set_name_size = 4,
         text_size = 0, show_percentage = FALSE) +
    annotate("text", x = -0.8, y = 0,
             label = paste(rna_sig,  collapse = "\n"),
             size = 2.3, color = "black") +
    annotate("text", x = 0.8, y = 0,
             label = paste(prot_sig, collapse = "\n"),
             size = 2.3, color = "black") +
    labs(title = toupper(.y)) +
    theme(plot.title = element_text(hjust = 0.5, face = "bold"))
})

patchwork::wrap_plots(venn_plots, nrow = 1)


# Analysis 2: Group-level RNA-protein logFC correlation ----
dedup_by_symbol <- function(df, expr_col) {
  df %>%
    filter(!is.na(hgnc_symbol)) %>%
    group_by(hgnc_symbol) %>%
    slice_max(!!sym(expr_col), n = 1, with_ties = FALSE) %>%
    ungroup()
} # remove duplicates in gene symbols from mapping

## Spearman correlation ----
get_lfc_correlation <- function(rna, prot, region_label) {
  joined <- inner_join(
    dedup_by_symbol(rna, "baseMean") %>% select(hgnc_symbol, rna_lfc = log2FoldChange),
    dedup_by_symbol(prot, "AveExpr")  %>% select(hgnc_symbol, prot_lfc = logFC),
    by = "hgnc_symbol"
  )
  
  test <- cor.test(joined$rna_lfc, joined$prot_lfc, method = "spearman")
  message(region_label, " | n genes: ", nrow(joined),
          " | Spearman r: ", round(test$estimate, 3),
          " | p: ", signif(test$p.value, 3))
  
  list(data = joined, r = test$estimate, p = test$p.value)
}

lfc_cor <- imap(regions, ~ get_lfc_correlation(rna_res[[.y]], prot_res[[.y]], toupper(.y)))

# Sensitivity test for APOE driven correlation 
lfc_cor_noAPOE <- imap(regions, ~ {
  rna <- rna_res[[.y]] %>% filter(hgnc_symbol != "APOE")
  prot <- prot_res[[.y]] %>% filter(hgnc_symbol != "APOE4")
  get_lfc_correlation(rna, prot, paste(toupper(.y), "excl. APOE"))
  }) # correlation is not driven by APOE4

## Scatter plot ----
make_lfc_scatter <- function(rna, prot, region_label) {
  rna_clean <- dedup_by_symbol(rna, "baseMean") %>% 
    select(hgnc_symbol, rna_lfc = log2FoldChange, rna_dir = direction)
  prot_clean <- dedup_by_symbol(prot, "AveExpr")%>% 
    select(hgnc_symbol, prot_lfc = logFC, prot_dir = direction)
  
  joined <- inner_join(rna_clean, prot_clean, by = "hgnc_symbol") %>%
    mutate(de_status = case_when(rna_dir!= "Not significant" ~ "DEG",
                                 prot_dir != "Not significant" ~ "DEP",
                                 TRUE ~ "Not significant"),
           label = if_else(de_status != "Not significant", hgnc_symbol, NA_character_))
  
  test <- cor.test(joined$rna_lfc, joined$prot_lfc, method = "spearman")
  r_val <- round(test$estimate, 3)
  p_val <- round(test$p.value, 4)
  stats_label <- paste0("Spearman r = ", r_val, ", p = ", p_val, "\nn = ", nrow(joined), " genes")
  
  ggplot(joined, aes(x = rna_lfc, y = prot_lfc)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.4) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "darkgrey", linewidth = 0.4) +
    geom_point(aes(color = de_status), alpha = 0.6, size = 1) +
    geom_smooth(method = "lm", color = "black", linewidth = 0.6, se = TRUE) +
    ggrepel::geom_text_repel(aes(label = label), na.rm = TRUE, size = 3,
                             color = "black", max.overlaps = Inf) +
    scale_color_manual(values = c("DEG" = "blue", "DEP" = "red", "Not significant" = "grey")) +
    annotate("text", x = Inf, y = -Inf, hjust = 1.05, vjust = -0.5,
             label = stats_label, size = 3.2, fontface = "italic") +
    labs(title = paste(region_label, "- RNA vs Protein log2FC (APOE4+ vs APOE4-)"),
         x = "RNA log2 fold change", y = "Protein log2 fold change", color = NULL) +
    theme_classic() +
    theme(legend.position = "bottom")
}

scatter_plots <- imap(regions, ~ make_lfc_scatter(rna_res[[.y]], prot_res[[.y]], toupper(.y)))

# ggsave(file.path(pic_output_dir, "DC_multiomic_dlpfc_logFC_scatter.png"), scatter_plots[["dlpfc"]], width = 8, height = 6)
# ggsave(file.path(pic_output_dir, "DC_multiomic_stg_logFC_scatter.png"),  scatter_plots[["stg"]],  width = 8, height = 6)

