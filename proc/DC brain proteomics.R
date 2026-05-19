# Set up ----
suppressMessages(library(gtsummary))
suppressMessages(library(here))
suppressMessages(library(limma))
suppressMessages(library(tidyverse))

pic_output_dir <- here("output", "DiverseCohorts", "pics")
table_output_dir <- here("output", "DiverseCohorts", "tables")

individual_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_individual_metadata.csv"),
                            header = T, strip.white = T, na.strings = c("", "NA"))
specimen_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_biospecimen_metadata.csv"),
                          header = T, strip.white = T, na.strings = c("", "NA"))
assay_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_assay_TMTproteomics_metadata.csv"),
                       header = T, strip.white = T, na.strings = c("", "NA"))
dlpfc_prot_raw <- read.csv(here("raw", "DiverseCohorts", "normAbundances_post-FP_frontal.csv"),
                           header = T, strip.white = T, na.strings = "NA", row.names = 1)
stg_prot_raw <- read.csv(here("raw", "DiverseCohorts", "normAbundances_post-FP_temporal.csv"),
                         header = T, strip.white = T, na.strings = "NA", row.names = 1)

# Preprocessing ----
## Metadata ----
meta_base <- individual_meta %>%
  select(individualID, sex, race, ageDeath, apoeGenotype, ADoutcome) %>%
  mutate(apoeGenotype = na_if(apoeGenotype, "missing or unknown"),
         ageDeath = na_if(ageDeath, "missing or unknown")) %>%
  filter(!is.na(apoeGenotype), !is.na(ageDeath), ADoutcome %in% c("AD", "Control")) %>%
  mutate(apoe4 = factor(apoeGenotype %in% c("24", "34", "44"), labels = c("APOE4-", "APOE4+")),
         apoeGenotype = factor(apoeGenotype, levels = c("22", "23", "33", "24", "34", "44")),
         age_death_num = as.numeric(if_else(ageDeath == "90+", "90", ageDeath)),
         ADoutcome = factor(ADoutcome, levels = c("Control", "AD")))

non_control <- assay_meta %>% 
  filter(isAssayControl == FALSE) %>%
  pull(specimenID)

# dlPFC metadata
# prefer Emory/emdp specimens when duplicates exist per individual
dlpfc_spec <- specimen_meta %>%
  filter(tissue == "dorsolateral prefrontal cortex",
         specimenID %in% non_control,
         dataGenerationSite == "Emory") %>%
  group_by(individualID) %>%
  mutate(is_dup = n() > 1) %>%
  filter(!is_dup | grepl("^emdp", specimenID)) %>%
  ungroup() %>%
  select(individualID, specimenID)

dlpfc_meta <- meta_base %>%
  inner_join(dlpfc_spec, by = "individualID") %>%
  inner_join(assay_meta %>% select(specimenID, batch), by = "specimenID")
# any(duplicated(dlpfc_meta$individualID))

# STG metadata
stg_spec <- specimen_meta %>%
  filter(tissue == "superior temporal gyrus",
         specimenID %in% non_control) %>%
  filter(n() == 1, .by = individualID) %>%
  select(individualID, specimenID)

stg_meta <- meta_base %>%
  inner_join(stg_spec, by = "individualID") %>%
  inner_join(assay_meta %>% select(specimenID, batch), by = "specimenID")
# any(duplicated(stg_meta$individualID))


## Proteomics data ----
process_expr_mat <- function(proteomic_dat, meta_dat) {
  median_impute <- function(x) {
    x[is.na(x)] <- median(x, na.rm = TRUE)
    x
  }
  
  proteomic_dat %>%
    select(all_of(meta_dat$specimenID)) %>% 
    mutate(across(where(is.numeric), ~ na_if(.x, 0))) %>% # code 0 as NA
    .[rowMeans(is.na(.)) < 0.3, ] %>% # remove proteins with >30% missing
    apply(1, median_impute) %>%  # median imputation 
    t() %>%
    log2() # log transformation
}

# dlPFC proteomics
dlpfc_expr <- process_expr_mat(dlpfc_prot_raw, dlpfc_meta)
# all(colnames(dlpfc_expr) == dlpfc_meta$specimenID)

# STG proteomics
stg_expr <- process_expr_mat(stg_prot_raw, stg_meta)
# all(colnames(stg_expr) == stg_meta$specimenID)


## Merge and save cleaned data ----
regions <- list(
  dlPFC = list(expr = dlpfc_expr, meta = dlpfc_meta),
  STG = list(expr = stg_expr, meta = stg_meta)
)

# saveRDS(regions, here("raw", "DiverseCohorts", "DC-proteomics.rds"))

rm(individual_meta, assay_meta, specimen_meta, meta_base, non_control,
   dlpfc_spec, stg_spec, dlpfc_prot_raw, stg_prot_raw, 
   dlpfc_meta, dlpfc_expr, stg_meta, stg_expr)

# Demographic ----
demographic <- imap(regions, ~ {
  .x$meta %>%
    select(ADoutcome, age_death_num, sex, apoeGenotype) %>%
    tbl_summary(
      by = ADoutcome,
      statistic = list(
        all_continuous() ~ "{mean} ± {sd}",
        all_categorical() ~ "{n} ({p}%)"
      ),
      digits = all_continuous() ~ 1,
      label = list(
        age_death_num ~ "Age at Death (top-coded 90)",
        sex ~ "Sex",
        apoeGenotype ~ "APOE genotype"
      ),
      missing = "no"
    ) %>%
    add_overall(last = TRUE) %>%
    bold_labels() %>%
    modify_spanning_header(everything() ~ paste0("**", .y, "**"))
}) %>%
  tbl_merge(tab_spanner = FALSE)

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "brain proteomics demographic table.png"))


# Variance Partition Analysis ----
factor_labels <- c(apoe4 = "APOE4 Status",
                   ADoutcome = "Diagnosis",
                   batch = "Batch",
                   sex = "Sex",
                   Residuals = "Residuals")

vp_form <- ~ (1 | batch) + (1 | apoe4) + (1 | ADoutcome) + (1 | sex)

vp_results <- imap(regions, ~ {
  vp_meta <- .x$meta %>% column_to_rownames("specimenID")
  vp_fit <- variancePartition::fitExtractVarPartModel(.x$expr, vp_form, vp_meta)
  
  vp_medians <- as.data.frame(vp_fit) %>%
    summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
    mutate(Factor = factor_labels[Factor]) %>%
    arrange(desc(Median_VarExplained))
  
  list(fit = vp_fit, medians = vp_medians)
})

# iwalk(vp_results, ~ write.csv(.x$medians, file.path(table_output_dir, paste0(tolower(.y), "_proteomics_varpart.csv")), row.names = FALSE))

# Batch confounding check
iwalk(regions, ~ {
  message(.y)
  .x$meta %>%
    group_by(batch, apoe4) %>%
    summarise(n = n(), .groups = "drop") %>%
    pivot_wider(names_from = apoe4, values_from = n, values_fill = 0) %>%
    print(n = Inf)
})


# Limma ----
## Model fitting & results ----
p_val_threshold <- 0.05

limma_results <- imap(regions, ~ {
  design_mat <- model.matrix(~ apoe4 + age_death_num + sex + ADoutcome + batch, # adjust for age, sex, diagnosis and batch
                             data = .x$meta)
  
  fit <- lmFit(.x$expr, design_mat) %>% eBayes()
  
  topTable(fit, coef = "apoe4APOE4+", number = Inf) %>%
    rownames_to_column(var = "Protein") %>%
    mutate(negLog10FDR = -log10(adj.P.Val),
           `Direction of change` = case_when(
             adj.P.Val < p_val_threshold & logFC > 0 ~ "Upregulated",
             adj.P.Val < p_val_threshold & logFC < 0 ~ "Downregulated",
             TRUE ~ "Not significant"))
})

sig_limma_results <- map(limma_results, ~ filter(.x, adj.P.Val < p_val_threshold))

# iwalk(sig_limma_results, ~ write.csv(.x, file.path(table_output_dir, paste0(tolower(.y), "_proteomics_limma.csv")), row.names = FALSE))

## Volcano plots ----
make_volcano <- function(res, sig_res, region_label, n_samples) {
  ggplot(res, aes(x = logFC, y = negLog10FDR)) +
    geom_point(aes(color = `Direction of change`), size = 1.2) +
    ggrepel::geom_text_repel(data = sig_res,
                             aes(label = Protein),
                             size = 3,
                             color = "black") +
    scale_color_manual(values = c("Upregulated" = "indianred",
                                  "Downregulated" = "steelblue",
                                  "Not significant" = "grey")) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "grey") +
    labs(title = sprintf("%d Differentially Expressed %s Proteins by APOE4 Status (n = %d)",
                         nrow(sig_res), region_label, n_samples),
         x = "log2 fold change",
         y = "-log10 adjusted P",
         color = "Direction of change") +
    theme_classic() +
    theme(legend.position = "bottom")
}

# DLPFC: exclude APOE4 from plot (top DEP distorts scale) 
p_volcano_dlpfc <- make_volcano(limma_results$dlPFC %>% filter(Protein != "APOE4|APOE4"),
                                sig_limma_results$dlPFC %>% filter(Protein != "APOE4|APOE4"),
                                "dlPFC", nrow(regions$dlPFC$meta)) +
  labs(subtitle = "Threshold: adj_p < 0.05 | APOE4 not shown on plot")

# ggsave(file.path(pic_output_dir, "dlpfc_limma_volcano.png"), p_volcano_dlpfc, width = 8, height = 6)

# STG
p_volcano_stg <- make_volcano(limma_results$STG, sig_limma_results$STG,
                              "STG", nrow(regions$STG$meta)) +
  labs(subtitle = "Threshold: adj_p < 0.05")

# ggsave(file.path(pic_output_dir, "stg_limma_volcano.png"), p_volcano_stg, width = 8, height = 6)

## PCA ----
dlpfc_deps <- sig_limma_results$dlPFC$Protein

pca_deps <- prcomp(t(regions$dlPFC$expr[dlpfc_deps, ]), center = TRUE, scale. = TRUE)

p_pca <- as.data.frame(pca_deps$x[, 1:2]) %>%
  bind_cols(regions$dlPFC$meta %>% select(apoe4, apoeGenotype)) %>%
  ggplot(aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA on %d dlPFC DEPs", length(dlpfc_deps)),
       color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "dlpfc_proteomics_pca.png"), p_pca, width = 7, height = 6)
