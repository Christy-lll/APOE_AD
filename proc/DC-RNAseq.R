# Set up ----
suppressMessages(library(tidyverse))
library(here)

pic_output_dir <- here("output", "pics", "DiverseCohorts")
table_output_dir <- here("output", "tables", "DiverseCohorts")

Columbia_dat <- read.csv(here("raw", "DiverseCohorts", "Columbia_counts_filtered.csv"),
                         header = T, strip.white = T, na.strings = "NA", row.names = 1)
Mayo_Emory_dat <- read.csv(here("raw", "DiverseCohorts", "Mayo_Emory_counts_filtered.csv"),
                      header = T, strip.white = T, na.strings = "NA", row.names = 1)
colnames(Mayo_Emory_dat) <- sub("^X", "", colnames(Mayo_Emory_dat))
MSSM_dat <- read.csv(here("raw", "DiverseCohorts", "MSSM_counts_filtered.csv"),
                      header = T, strip.white = T, na.strings = "NA", row.names = 1)
colnames(MSSM_dat) <- sub("^X", "", colnames(MSSM_dat))
Rush_dat <- read.csv(here("raw", "DiverseCohorts", "Rush_counts_filtered_1_.csv"),
                     header = T, strip.white = T, na.strings = "NA", row.names = 1)

individual_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_individual_metadata.csv"), 
                            header = T, strip.white = T, na.strings = c("", "NA"))
specimen_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_biospecimen_metadata.csv"), 
                          header = T, strip.white = T, na.strings = c("", "NA"))

# Preprocessing ----
## RNA-seq data ----
# Merge RNA-seq data from all sites 
rna_merged <- list(Mayo_Emory_dat, Columbia_dat, MSSM_dat, Rush_dat) %>%
  lapply(rownames_to_column, "gene") %>%
  purrr::reduce(inner_join, by = "gene") %>%
  column_to_rownames("gene") # 46527 overlap genes in all sites

# Filter specimen_meta to RNA-seq only and present in actual data
rna_specimen_meta <- specimen_meta %>%
  filter(assay == "rnaSeq",
         specimenID %in% colnames(rna_merged),
         individualID %in% individual_meta$individualID)

# Handle duplicates (same individual, same region, multiple sites)
dup_check <- rna_specimen_meta %>%
  group_by(individualID, tissue) %>%
  filter(n() > 1) %>%
  ungroup()

message(n_distinct(dup_check$individualID), " individuals with duplicate region-site entries, ",
        nrow(dup_check), " affected rows total.")

site_priority <- c("Mayo", "RUSH", "NYGC") # double check this preference

rna_specimen_meta <- rna_specimen_meta %>%
  mutate(site_rank = match(dataGenerationSite, site_priority)) %>%
  group_by(individualID, tissue) %>%
  slice_min(site_rank, n = 1, with_ties = FALSE) %>%
  ungroup() %>%
  dplyr::select(-site_rank)

# Split RNA-seq by region and convert to matrix
regions <- rna_specimen_meta %>%
  filter(tissue != "temporal pole") %>% # no observation for this region
  pull(tissue) %>%
  unique()

rna_by_region <- setNames(
  lapply(regions, function(reg) {
    spec_ids <- rna_specimen_meta$specimenID[rna_specimen_meta$tissue == reg]
    m <- as.matrix(rna_merged[, spec_ids, drop = FALSE])
    mode(m) <- "numeric"
    m
  }),
  gsub(" ", "_", regions)
)

## Metadata ----
# Build one metadata table per region
meta_base <- individual_meta %>%
  dplyr::select(individualID, sex, race, ageDeath, apoeGenotype, ADoutcome) %>%
  mutate(apoeGenotype = na_if(apoeGenotype, "missing or unknown"),
         ageDeath = na_if(ageDeath, "missing or unknown")) %>%
  filter(!is.na(apoeGenotype),
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
        sex = factor(sex),                          
        dataGenerationSite = factor(dataGenerationSite)
        )
  }),
  gsub(" ", "_", regions)
)

# Drops RNA specimens that failed metadata filters
rna_by_region <- setNames(
  lapply(regions, function(reg) {
    reg_name <- gsub(" ", "_", reg)
    spec_ids <- meta_by_region[[reg_name]]$specimenID
    rna_by_region[[reg_name]][, spec_ids, drop = FALSE]
  }),
  gsub(" ", "_", regions)
)

## Sanity check & save data----
# No duplicate individualIDs within each region
invisible(lapply(names(meta_by_region), function(reg) {
  dups <- meta_by_region[[reg]]$individualID[duplicated(meta_by_region[[reg]]$individualID)]
  if (length(dups) > 0) {
    warning(reg, " has duplicate individualIDs: ", paste(dups, collapse = ", "))
  } else {
    message(reg, ": no duplicate individualIDs")
  }
}))

# RNA columns match metadata rows exactly for each region
invisible(lapply(names(rna_by_region), function(reg) {
  rna_specs  <- colnames(rna_by_region[[reg]])
  meta_specs <- meta_by_region[[reg]]$specimenID
  if (!identical(rna_specs, meta_specs)) {
    warning(reg, ": RNA columns and metadata specimenIDs do not match")
  } else {
    message(reg, ": RNA and metadata aligned (", length(rna_specs), " specimens)")
  }
}))

# Check dataGenerationSite
meta_by_region[["dorsolateral_prefrontal_cortex"]] %>%
  dplyr::count(.data$dataGenerationSite, .data$apoe4) %>%
  dplyr::group_by(.data$dataGenerationSite) %>%
  dplyr::mutate(prop = n / sum(n)) %>%
  tidyr::pivot_wider(names_from = apoe4, values_from = prop, values_fill = 0)

# Remove unused dataset
rm(meta_base, individual_meta, specimen_meta, dup_check,
   Columbia_dat, Mayo_Emory_dat, MSSM_dat, Rush_dat)

# Save cleaned data ----
# saveRDS(rna_by_region,  here("raw", "DiverseCohorts", "rna_count_by_region.rds"))
# saveRDS(meta_by_region, here("raw", "DiverseCohorts", "rna_meta_by_region.rds"))

# Demographic ----
library(gtsummary)

demographic <- imap(meta_by_region, function(meta, region_name) {
  meta %>%
    select(apoe4, ADoutcome, age_death_num, sex, apoeGenotype) %>%
    tbl_summary(
      by = apoe4,
      statistic = list(
        all_continuous() ~ "{mean} ± {sd}",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits= all_continuous() ~ 1,
      label = list(
        ADoutcome ~ "Diagnosis",
        age_death_num ~ "Age at Death",
        apoeGenotype ~ "APOE Genotype",
        sex ~ "Sex"
      ),
      missing = "no"
    ) %>%
    add_overall(last = TRUE) %>%
    bold_labels() %>%
    modify_spanning_header(everything() ~ paste("**", gsub("_", " ", region_name), "**"))
}) %>%
  tbl_merge(tab_spanner = FALSE)

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "DC_RNAseq_demographic table.png"))


# Variance Partition Analysis (DLPFC only) ----
## Analysis ----
library(variancePartition)
library(BiocParallel)

param <- SnowParam(workers = parallel::detectCores() - 1, type = "SOCK")
register(param)

vp_form <- ~ (1|ADoutcome) + (1|apoe4) + (1|dataGenerationSite) + (1|sex)

reg <- "dorsolateral_prefrontal_cortex"

meta_dlpfc <- as.data.frame(meta_by_region[[reg]])
rownames(meta_dlpfc) <- meta_dlpfc$specimenID
expr_dlpfc <- rna_by_region[[reg]][, meta_dlpfc$specimenID]

vp_dlpfc <- fitExtractVarPartModel(expr_dlpfc, vp_form, meta_dlpfc, BPPARAM = param)

## Results ----
factor_labels <- c(ADoutcome = "Diagnosis",
                   apoe4 = "APOE4 Status",
                   dataGenerationSite = "Data Generation Site",
                   sex = "Sex",
                   Residuals = "Residuals")

vp_medians_dlpfc <- as.data.frame(vp_dlpfc) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
  mutate(Region = reg,
         Factor = factor_labels[Factor],
         Median_VarExplained = round(Median_VarExplained, 4)) %>%
  arrange(desc(Median_VarExplained))

# write.csv(vp_medians_dlpfc, file.path(table_output_dir, "DC_RNAseq_dlpfc_varpart.csv"), row.names = FALSE)

# DESeq2 ----
suppressMessages(library(DESeq2))
suppressMessages(library(biomaRt))

## Analysis ----
design_formula <- ~ ADoutcome + apoe4 + sex + age_death_scaled + dataGenerationSite # adjust for age and sex

min_group_size <- lapply(meta_by_region, function(meta)
  min(table(meta$apoe4)))

dds_by_region <- mapply(function(meta, count, min_n) {
  meta <- as.data.frame(meta)
  rownames(meta) <- meta$specimenID
  meta$age_death_scaled <- scale(meta$age_death_num)[, 1]  # centre and scale age
  count <- round(count[, meta$specimenID]) # round to nearest integer
  mode(count) <- "integer"
  stopifnot(all(colnames(count) == meta$specimenID))
  
  dds <- DESeqDataSetFromMatrix(countData = count,
                                colData = meta,
                                design = design_formula)
  
  n_before <- nrow(dds)
  dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ]
  n_after <- nrow(dds)
  message("Region: ", unique(meta$tissue),
          " | Before: ", n_before,
          " | After: ", n_after,
          " | Removed: ", n_before - n_after)
  
  DESeq(dds)
}, meta_by_region, rna_by_region, min_group_size, SIMPLIFY = FALSE)

## Results ----
p_val_threshold <- 0.05

res_by_region <- lapply(dds_by_region, function(dds) {
  results(dds,
          contrast = c("apoe4", "APOE4pos", "APOE4neg"),
          alpha = p_val_threshold) %>%
    as.data.frame() %>%
    arrange(padj) %>%
    rownames_to_column("ensembl_gene_id") %>%
    mutate(
      ensembl_gene_id = sub("\\..*", "", ensembl_gene_id), # strip version suffix for ID mapping
      negLog10FDR = -log10(padj),
      direction = case_when(
        padj < p_val_threshold & log2FoldChange > 0 ~ "Upregulated",
        padj < p_val_threshold & log2FoldChange < 0 ~ "Downregulated",
        TRUE ~ "Not significant"
      )
    )
})

## Gene ID mapping ----
all_ensembl_ids <- lapply(res_by_region, `[[`, "ensembl_gene_id") %>%
  unlist() %>%
  unique() 

mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

id_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = "ensembl_gene_id",
                values = all_ensembl_ids,
                mart = mart) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  mutate(hgnc_symbol = if_else(hgnc_symbol == "", ensembl_gene_id, hgnc_symbol))

res_by_region <- lapply(res_by_region, function(res)
  left_join(res, id_map, by = "ensembl_gene_id"))

sig_res_by_region <- lapply(res_by_region, function(res)
  filter(res, direction != "Not significant"))

# for (region in names(res_by_region)) {
#   region_file <- gsub(" ", "_", region)
#   write.csv(res_by_region[[region]],
#             file.path(table_output_dir, paste0(region_file, "_DC_RNAseq_deseq2.csv")),
#             row.names = FALSE)
# }

## Volcano plots ----
make_volcano <- function(res, sig_res, region_name) {
  n_up   <- sum(res$direction == "Upregulated",   na.rm = TRUE)
  n_down <- sum(res$direction == "Downregulated", na.rm = TRUE)
  
  ggplot(res, aes(x = log2FoldChange, y = negLog10FDR)) +
    geom_point(aes(color = direction), alpha = 0.6, size = 1.2) +
    ggrepel::geom_text_repel(data = sig_res, aes(label = hgnc_symbol),
                             size = 3, color = "black") +
    scale_color_manual(values = c("Upregulated" = "red",
                                  "Downregulated" = "blue",
                                  "Not significant" = "grey")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "grey") +
    labs(title = gsub("_", " ", region_name),
         subtitle = paste0("adj_p < 0.05  |  Up: ", n_up, "  Down: ", n_down),
         x = "log2 fold change",
         y = "-log10 adjusted P",
         colour = "Direction of change") +
    theme_classic() +
    theme(legend.position = "bottom")
}

volcano_plots <- mapply(make_volcano,
                        res = res_by_region,
                        sig_res = sig_res_by_region,
                        region_name = names(res_by_region),
                        SIMPLIFY = FALSE)

library(patchwork)
combined_volcano <- patchwork::wrap_plots(volcano_plots, nrow = 1) +
  plot_annotation(title = "DE Genes by APOE4 Status",
                  theme = theme(plot.title = element_text(hjust = 0.5))) 

# ggsave(file.path(pic_output_dir, "DC_RNAseq_deseq2_volcano.png"),
#        plot = combined_volcano,
#        width = 8 * length(volcano_plots),
#        height = 9)
