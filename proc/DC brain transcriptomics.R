# Set up ----
suppressMessages(library(biomaRt))
suppressMessages(library(DESeq2))
suppressMessages(library(edgeR))
suppressMessages(library(gtsummary))
suppressMessages(library(here))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))
suppressMessages(library(variancePartition))

pic_output_dir <- here("output", "DiverseCohorts", "pics")
table_output_dir <- here("output", "DiverseCohorts", "tables")

rna_files <- list(
  Mayo_Emory = "Mayo_Emory_counts_filtered.csv",
  Columbia = "Columbia_counts_filtered.csv",
  MSSM = "MSSM_counts_filtered.csv",
  Rush = "Rush_counts_filtered_1_.csv"
)

rna_count <- lapply(rna_files, function(f) {
  dat <- read.csv(here("raw", "DiverseCohorts", f),
                  header = T, strip.white = T, na.strings = "NA", row.names = 1)
  colnames(dat) <- sub("^X", "", colnames(dat)) # strip leading X added by R when column names start with digits
  dat
})

individual_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_individual_metadata.csv"),
                            header = T, strip.white = T, na.strings = c("", "NA"))
specimen_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_biospecimen_metadata.csv"),
                          header = T, strip.white = T, na.strings = c("", "NA"))


# Preprocessing ----
## RNA-seq and metadata ----
# merge RNA-seq data across sites
rna_merged <- rna_count %>%
  lapply(rownames_to_column, "gene") %>%
  purrr::reduce(inner_join, by = "gene") %>% # keep genes present in all 4 sites
  column_to_rownames("gene")

# filter specimen metadata  
rna_specimen_meta <- specimen_meta %>% # keep specimens with both RNA-seq and individual metadata
  filter(assay == "rnaSeq", 
         specimenID %in% colnames(rna_merged),
         individualID %in% individual_meta$individualID)

site_priority <- c("Mayo", "RUSH", "NYGC")
rna_specimen_meta <- rna_specimen_meta %>%
  mutate(site_rank = match(dataGenerationSite, site_priority)) %>%
  group_by(individualID, tissue) %>%
  slice_min(site_rank, n = 1, with_ties = FALSE) %>% # handle duplication (same individual+region appears in multiple sites)
  ungroup() %>%
  dplyr::select(-site_rank)

# split RNA-seq data by brain region 
regions <- rna_specimen_meta %>%
  filter(tissue != "temporal pole") %>% # exclude temporal pole: no observations
  pull(tissue) %>%
  unique() 

region_keys <- c(
  "dorsolateral prefrontal cortex" = "dlPFC",
  "superior temporal gyrus" = "STG",
  "caudate nucleus" = "CN"
)

rna_by_region <- setNames(
  lapply(regions, function(reg) {
    spec_ids <- rna_specimen_meta$specimenID[rna_specimen_meta$tissue == reg]
    m <- as.matrix(rna_merged[, spec_ids, drop = FALSE])
    mode(m) <- "numeric"
    m
  }),
  region_keys[regions]
)

# merge metadata and split by region
meta_base <- individual_meta %>%
  dplyr::select(individualID, sex, race, ageDeath, apoeGenotype, ADoutcome) %>%
  mutate(apoeGenotype = na_if(apoeGenotype, "missing or unknown"),
         ageDeath = na_if(ageDeath, "missing or unknown")) %>%
  filter(!is.na(apoeGenotype), # drop samples with key NAs
         !is.na(ageDeath),
         ADoutcome %in% c("AD", "Control"))

meta_by_region <- setNames(
  lapply(regions, function(reg) {
    rna_specimen_meta %>%
      filter(tissue == reg) %>%
      dplyr::select(individualID, specimenID, tissue, dataGenerationSite) %>%
      inner_join(meta_base, by = "individualID") %>%
      mutate(
        ADoutcome = factor(ADoutcome, levels = c("Control", "AD")),
        apoe4 = factor(apoeGenotype %in% c("24", "34", "44"), labels = c("APOE4neg", "APOE4pos")),
        apoeGenotype = factor(apoeGenotype, levels = c("22", "23", "33", "24", "34", "44")),
        age_death_num = as.numeric(if_else(ageDeath == "90+", "90", ageDeath)),
        age_death_scaled = as.numeric(scale(age_death_num)), 
        sex = factor(sex),
        dataGenerationSite = factor(dataGenerationSite)
      )
  }),
  region_keys[regions]
)

# align RNA-seq data to metadata 
rna_by_region <- map2(rna_by_region, meta_by_region, ~ .x[, .y$specimenID, drop = FALSE])

## Duplication & alignment checks ----
walk(region_keys, function(reg) {
  dups <- meta_by_region[[reg]]$individualID[duplicated(meta_by_region[[reg]]$individualID)]
  if (length(dups) > 0) {
    warning(reg, " has duplicate individualIDs: ", paste(dups, collapse = ", "))
  } else {
    message(reg, ": no duplicate individualIDs")
  }
  
  if (!identical(colnames(rna_by_region[[reg]]), meta_by_region[[reg]]$specimenID)) {
    warning(reg, ": RNA columns and metadata specimenIDs do not match")
  } else {
    message(reg, ": RNA and metadata aligned (", ncol(rna_by_region[[reg]]), " specimens)")
  }
})

## Save cleaned data ----
# saveRDS(list(count = rna_by_region, meta = meta_by_region), here("raw", "DiverseCohorts", "DC-transcriptomics.rds"))
rm(rna_count, rna_files, rna_merged, meta_base, individual_meta, specimen_meta, rna_specimen_meta)


# Demographic ----
demographic <- imap(meta_by_region, function(meta, region_name) {
  meta %>%
    dplyr::select(ADoutcome, age_death_num, sex, apoeGenotype) %>%
    tbl_summary(
      by = ADoutcome,
      statistic = list(
        all_continuous() ~ "{mean} ± {sd}",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = all_continuous() ~ 1,
      label = list(
        age_death_num ~ "Age at Death (top-coded 90)",
        apoeGenotype ~ "APOE Genotype",
        sex ~ "Sex"
      ),
      missing = "no"
    ) %>%
    add_overall(last = TRUE) %>%
    bold_labels() %>%
    modify_spanning_header(everything() ~ paste0("**", region_name, "**"))
}) %>%
  tbl_merge(tab_spanner = FALSE)

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "brain transcriptomics demographic table.png"))


# Variance Partition Analysis ----
# select region for analysis 
reg <- "dlPFC"

meta_varpart <- meta_by_region[[reg]] %>%
  as.data.frame() %>%
  column_to_rownames("specimenID")

# voom normalization on count data
dge_varpart <- DGEList(counts = rna_by_region[[reg]])
dge_varpart <- calcNormFactors(dge_varpart, method = "TMM") # normalise 
design_varpart <- model.matrix(~ ADoutcome + apoe4 + sex + age_death_scaled + dataGenerationSite,
                               data = meta_varpart)
voom_varpart <- voom(dge_varpart, design = design_varpart)

# variance partition 
param <- BiocParallel::SnowParam(workers = parallel::detectCores() - 1, type = "SOCK")
BiocParallel::register(param)

factor_labels <- c(apoe4 = "APOE4 Status",
                   ADoutcome = "Diagnosis",
                   dataGenerationSite = "Data Generation Site",
                   sex = "Sex",
                   Residuals = "Residuals")

vp_form <- ~ (1|ADoutcome) + (1|apoe4) + (1|dataGenerationSite) + (1|sex)
vp_fit <- fitExtractVarPartModel(voom_varpart, vp_form, meta_varpart, BPPARAM = param)

vp_medians <- as.data.frame(vp_fit) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
  mutate(Factor = factor_labels[Factor]) %>%
  arrange(desc(Median_VarExplained))

# write.csv(vp_medians, file.path(table_output_dir, "dlPFC_transcriptomics_varpart.csv"), row.names = FALSE)


# DESeq2 ----
## Model fitting ----
design_formula <- ~ ADoutcome + apoe4 + sex + age_death_scaled + dataGenerationSite

dds_by_region <- map2(meta_by_region, rna_by_region, function(meta, count) {
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$specimenID
  
  count <- round(count[, meta$specimenID])
  mode(count) <- "integer"
  stopifnot(all(colnames(count) == meta$specimenID))
  
  min_n <- min(table(meta$apoe4)) 
  
  dds <- DESeqDataSetFromMatrix(countData = count, colData = meta, design = design_formula)
  
  n_before <- nrow(dds)
  dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ] # filter: genes with >= 10 counts in at least as many samples as the smallest APOE4 group
  message("Region: ", unique(meta$tissue),
          " | Before: ", n_before,
          " | After: ", nrow(dds),
          " | Removed: ", n_before - nrow(dds))
  
  DESeq(dds)
})

## Results ----
p_val_threshold <- 0.05

res_by_region <- map(dds_by_region, function(dds) {
  results(dds,
          contrast = c("apoe4", "APOE4pos", "APOE4neg"),
          alpha = p_val_threshold) %>%
    as.data.frame() %>%
    arrange(padj) %>%
    rownames_to_column("ensembl_gene_id") %>%
    mutate(
      ensembl_gene_id = sub("\\..*", "", ensembl_gene_id), # strip Ensembl version suffix for ID mapping
      negLog10FDR = -log10(padj),
      direction = case_when(
        padj < p_val_threshold & log2FoldChange > 0 ~ "Upregulated",
        padj < p_val_threshold & log2FoldChange < 0 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
})

## Gene ID mapping ----
all_ensembl_ids <- map(res_by_region, "ensembl_gene_id") %>% unlist() %>% unique()

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

res_by_region <- map(res_by_region, ~ left_join(.x, id_map, by = "ensembl_gene_id"))

sig_res_by_region <- map(res_by_region, ~ filter(.x, padj < p_val_threshold))

# iwalk(sig_res_by_region, ~ write.csv(.x, file.path(table_output_dir, paste0(tolower(.y), "_transcriptomics_deseq2.csv")), row.names = FALSE))

## Volcano plots ----
make_volcano <- function(res, sig_res, region_name) {
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
    labs(title = region_name,
         subtitle = sprintf("adj_p < %.2f  |  Up: %d  Down: %d", p_val_threshold, n_up, n_down),
         x = "log2 fold change",
         y = "-log10 adjusted P",
         color = "Direction of change") +
    theme_classic() +
    theme(legend.position = "bottom")
}

volcano_plots <- imap(res_by_region, ~ make_volcano(.x, sig_res_by_region[[.y]], .y))

combined_volcano <- wrap_plots(volcano_plots, nrow = 1) +
  plot_annotation(title = "DE Genes by APOE4 Status",
                  theme = theme(plot.title = element_text(hjust = 0.5)))

# ggsave(file.path(pic_output_dir, "brain_deseq2_volcano.png"), plot = combined_volcano, width = 8 * length(volcano_plots), height = 9)
