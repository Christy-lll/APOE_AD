# Set up ----
suppressMessages(library(BiocParallel))
suppressMessages(library(biomaRt))
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(gtsummary))
suppressMessages(library(here))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(variancePartition))

pic_output_dir <- here("output", "ROSMAP", "pics")
table_output_dir <- here("output", "ROSMAP", "tables")

expr_dat <- read_tsv(here("raw", "ROSMAP", "ROSMAP_Normalized_counts__CQN_.tsv")) %>%
  column_to_rownames(var = "feature")
count_dat <- read_tsv(here("raw", "ROSMAP", "ROSMAP_Filtered_counts__greater_than_1cpm_.tsv")) %>%
  column_to_rownames(var = "feature")

ID_meta <- read.csv(here("raw", "ROSMAP", "ROSMAP_biospecimen_metadata.csv"),
                    header = T, strip.white = T, na.strings = "NA")
rna_meta <- read.csv(here("raw", "ROSMAP", "ROSMAP_assay_rnaSeq_metadata.csv"),
                     header = T, strip.white = T, na.strings = "NA")
sample_meta <- read.csv(here("raw", "ROSMAP", "ROSMAP_clinical_1_.csv"),
                        header = T, strip.white = T, na.strings = "NA")


# Preprocessing ----
## Metadata ----
meta_dat <- rna_meta %>%
  left_join(ID_meta[, c("specimenID", "individualID", "tissue", "assay")], by = "specimenID") %>%
  filter(assay == "rnaSeq") %>%
  left_join(sample_meta, by = "individualID") %>%
  filter(specimenID %in% colnames(count_dat)) %>%
  mutate(
    apoe4 = factor(grepl("4", apoe_genotype), labels = c("APOE4neg", "APOE4pos")),
    apoe_genotype = factor(apoe_genotype, levels = c("22", "23", "33", "24", "34", "44")),
    cogdx = factor(case_when(cogdx == 1 ~ "NCI", cogdx == 2 ~ "MCI", cogdx == 4 ~ "AD"),
                   levels = c("NCI", "MCI", "AD")),
    msex = factor(msex, levels = c("0", "1"), labels = c("Female", "Male")),
    age_death_num = as.numeric(if_else(age_death == "90+", "90", age_death)),
    age_death_scaled = as.numeric(scale(age_death_num)), 
    sequencingBatch = factor(sequencingBatch)
  ) %>%
  filter(!is.na(cogdx),
         sequencingBatch != "0, 6, 7") %>% # drop wired batch with only 1 observation
  mutate(sequencingBatch = droplevels(sequencingBatch)) %>%
  group_by(individualID, tissue) %>%
  slice_max(RIN, n = 1, with_ties = FALSE) %>% # keep highest RIN specimen when individual has multiple specimens per tissue
  ungroup()

## Split data by tissue ----
tissue_keys <- c(
  "dorsolateral prefrontal cortex" = "dlPFC",
  "Head of caudate nucleus" = "hCN",
  "posterior cingulate cortex" = "PCC"
)

meta_by_tissue <- split(meta_dat, meta_dat$tissue) %>%
  setNames(tissue_keys[names(.)])

count_by_tissue <- lapply(meta_by_tissue, function(meta) {
  count_dat[, meta$specimenID, drop = FALSE]
})

expr_by_tissue <- lapply(meta_by_tissue, function(meta) {
  expr_dat[, meta$specimenID, drop = FALSE]
})

rm(count_dat, expr_dat, ID_meta, meta_dat, rna_meta, sample_meta)

## Duplication & alignment checks ----
walk(names(meta_by_tissue), function(tissue) {
  dups <- meta_by_tissue[[tissue]]$individualID[duplicated(meta_by_tissue[[tissue]]$individualID)]
  if (length(dups) > 0) {
    warning(tissue, " has duplicate individualIDs: ", paste(dups, collapse = ", "))
  } else {
    message(tissue, ": no duplicate individualIDs")
  }
  
  if (!identical(colnames(count_by_tissue[[tissue]]), meta_by_tissue[[tissue]]$specimenID)) {
    warning(tissue, ": count columns and metadata specimenIDs do not match")
  } else {
    message(tissue, ": count and metadata aligned (", ncol(count_by_tissue[[tissue]]), " specimens)")
  }
})


# Demographic ----
demographic <- imap(meta_by_tissue, function(meta, tissue_name) {
  meta %>%
    dplyr::select(cogdx, age_death_num, pmi, msex, apoe_genotype) %>%
    tbl_summary(
      by = cogdx,
      statistic = list(
        all_continuous() ~ "{mean} ± {sd}",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = all_continuous() ~ 1,
      label = list(
        age_death_num ~ "Age at Death (top-coded 90)",
        msex ~ "Sex",
        apoe_genotype ~ "APOE Genotype",
        pmi ~ "Post-Mortem Interval"
      ),
      missing = "no"
    ) %>%
    add_overall(last = TRUE) %>%
    bold_labels() %>%
    modify_spanning_header(everything() ~ paste0("**", tissue_name, "**"))
}) %>%
  tbl_merge(tab_spanner = FALSE)

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "brain transcriptomics demographic table.png"))


# Variance Partition Analysis ----
# select region for analysis 
tissue <- "dlPFC"

meta_varpart <- meta_by_tissue[[tissue]] %>%
  as.data.frame() %>%
  column_to_rownames("specimenID")

# variance partition 
param <- SnowParam(workers = parallel::detectCores() - 1, type = "SOCK")
register(param)

factor_labels <- c(apoe4 = "APOE4 Status",
                   cogdx = "Diagnosis",
                   sequencingBatch = "Sequencing Batch",
                   msex = "Sex",
                   Residuals = "Residuals")

vp_form <- ~ (1|cogdx) + (1|apoe4) + (1|sequencingBatch) + (1|msex)
vp_fit <- fitExtractVarPartModel(expr_by_tissue[[tissue]], vp_form, meta_varpart, BPPARAM = param)

vp_medians <- as.data.frame(vp_fit) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
  mutate(Factor = factor_labels[Factor]) %>%
  arrange(desc(Median_VarExplained))

# write.csv(vp_medians, file.path(table_output_dir, "dlpfc_transcriptomics_varpart.csv"), row.names = FALSE)


# DESeq2 ----
## Model fitting ----
design_formula <- ~ cogdx + apoe4 + msex + age_death_scaled + sequencingBatch

dds_by_tissue <- imap(meta_by_tissue, function(meta, tissue_name) {
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$specimenID
  
  count <- count_by_tissue[[tissue_name]][, meta$specimenID]
  stopifnot(all(colnames(count) == meta$specimenID))
  
  min_n <- min(table(meta$apoe4))
  
  dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = design_formula)
  
  n_before <- nrow(dds)
  dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ] # filter: >= 10 counts in at least as many samples as the smallest APOE4 group
  message("Tissue: ", tissue_name,
          " | Before: ", n_before,
          " | After: ", nrow(dds),
          " | Removed: ", n_before - nrow(dds))
  
  DESeq(dds)
})

## Results ----
p_val_threshold <- 0.05

res_by_tissue <- map(dds_by_tissue, function(dds) {
  results(dds,
          contrast = c("apoe4", "APOE4pos", "APOE4neg"),
          alpha = p_val_threshold) %>%
    as.data.frame() %>%
    arrange(padj) %>%
    rownames_to_column("ensembl_gene_id") %>%
    mutate(negLog10FDR = -log10(padj),
           direction = case_when(padj < p_val_threshold & log2FoldChange > 0 ~ "Upregulated",
                                 padj < p_val_threshold & log2FoldChange < 0 ~ "Downregulated",
                                 TRUE ~ "Not significant"))
})

## Gene ID mapping ----
all_ensembl_ids <- map(res_by_tissue, "ensembl_gene_id") %>% unlist() %>% unique()

id_map <- tryCatch({
  mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")
  getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
        filters = "ensembl_gene_id",
        values = all_ensembl_ids,
        mart = mart) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    mutate(hgnc_symbol = if_else(hgnc_symbol == "", ensembl_gene_id, hgnc_symbol))
}, error = function(e) {
  message("Ensembl unavailable, falling back to org.Hs.eg.db: ", conditionMessage(e))
  suppressMessages(library(org.Hs.eg.db))
  AnnotationDbi::select(org.Hs.eg.db,
                        keys = all_ensembl_ids,
                        columns = "SYMBOL",
                        keytype = "ENSEMBL") %>%
    dplyr::rename(ensembl_gene_id = ENSEMBL, hgnc_symbol = SYMBOL) %>%
    distinct(ensembl_gene_id, .keep_all = TRUE) %>%
    mutate(hgnc_symbol = if_else(is.na(hgnc_symbol), ensembl_gene_id, hgnc_symbol))
})

res_by_tissue <- map(res_by_tissue, ~ left_join(.x, id_map, by = "ensembl_gene_id"))

sig_res_by_tissue <- map(res_by_tissue, ~ filter(.x, padj < p_val_threshold))

# iwalk(sig_res_by_tissue, ~ write.csv(.x, file.path(table_output_dir, paste0(tolower(.y), "_transcriptomics_deseq2.csv")), row.names = FALSE))

## Volcano plots ----
make_volcano <- function(res, sig_res, tissue_name) {
  n_up <- sum(sig_res$direction == "Upregulated", na.rm = TRUE)
  n_down <- sum(sig_res$direction == "Downregulated", na.rm = TRUE)
  
  ggplot(res, aes(x = log2FoldChange, y = negLog10FDR)) +
    geom_point(aes(color = direction), size = 1.2) +
    ggrepel::geom_text_repel(data = sig_res, aes(label = hgnc_symbol),
                             size = 3, color = "black") +
    scale_color_manual(values = c("Upregulated" = "indianred",
                                  "Downregulated" = "steelblue",
                                  "Not significant" = "grey")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "grey") +
    labs(title = tissue_name,
         subtitle = sprintf("adj_p < %.2f  |  Up: %d  Down: %d", p_val_threshold, n_up, n_down),
         x = "log2 fold change",
         y = "-log10 adjusted P",
         color = "Direction of change") +
    theme_classic() +
    theme(legend.position = "bottom")
}

volcano_plots <- imap(res_by_tissue, ~ make_volcano(.x, sig_res_by_tissue[[.y]], .y))

combined_volcano <- wrap_plots(volcano_plots, nrow = 1) +
  plot_annotation(title = "DE Genes by APOE4 Status",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

# ggsave(file.path(pic_output_dir, "brain_deseq2_volcano.png"), plot = combined_volcano, width = 8 * length(volcano_plots), height = 9)
