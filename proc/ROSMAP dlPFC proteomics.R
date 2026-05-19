# Set up ----
suppressMessages(library(gtsummary))
suppressMessages(library(here))
suppressMessages(library(limma))
suppressMessages(library(tidyverse))

pic_output_dir <- here("output", "ROSMAP", "pics")
table_output_dir <- here("output", "ROSMAP", "tables")

proteomic_dat <- read.csv(here("raw", "ROSMAP", "A1.PD-RAW_normalized_abundance-11672x500TMTchannels-no_batch_correction_1_.csv"),
                          header = T, strip.white = T, na.strings = c("", "NA"), row.names = 1)
protein_meta <- read.csv(here("raw", "ROSMAP", "ROSMAP_assay_proteomics_TMTquantitation_metadata.csv"),
                         header = T, strip.white = T, na.strings = c("", "NA"))
sample_meta <- read.csv(here("raw", "ROSMAP", "ROSMAP_clinical_1_.csv"),
                        header = T, strip.white = T, na.strings = c("", "NA"))


# Preprocessing ----
## Metadata ----
meta_dat <- protein_meta %>%
  filter(!isAssayControl) %>%
  mutate(individualID = str_extract(specimenID, "[^.]+$")) %>%
  select(specimenID, individualID, batch, batchChannel) %>%
  left_join(sample_meta, by = "individualID") %>%
  mutate(
    apoe4 = factor(grepl("4", apoe_genotype), labels = c("APOE4-", "APOE4+")),
    apoe_genotype = factor(apoe_genotype, levels = c("22", "23", "33", "24", "34", "44")),
    cogdx = case_when(cogdx == 1 ~ "NCI", cogdx == 2 ~ "MCI", cogdx == 4 ~ "AD"),
    cogdx = factor(cogdx, levels = c("NCI", "MCI", "AD")),
    msex = factor(msex, levels = c("0", "1"), labels = c("Female", "Male")),
    age_death_num = as.numeric(if_else(age_death == "90+", "90", age_death))
  ) %>%
  filter(!is.na(cogdx))

## Proteomic data ----
median_impute <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

expr_mat <- as.matrix(proteomic_dat)[, meta_dat$batchChannel] %>% # match metadata
  .[rowMeans(is.na(.)) < 0.3, ] %>% # drop proteins with >30% NA
  apply(1, median_impute) %>% # median imputation
  t() %>%
  log2() # log transformation 

rm(proteomic_dat, protein_meta, sample_meta)


# Demographic ----
demographic <- meta_dat %>%
  dplyr::select(cogdx, age_death_num, msex, educ, pmi, apoe_genotype) %>%
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
      educ ~ "Education (years)",
      pmi ~ "Post-Mortem Interval",
      apoe_genotype ~ "APOE genotype"
    ),
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels()

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "dlpfc proteomics demographic table.png"))


# Variance Partition Analysis ----
factor_labels <- c(apoe4 = "APOE4 Status",
                   batch = "Batch",
                   cogdx = "Diagnosis",
                   msex = "Sex",
                   Residuals = "Residuals")

vp_form <- ~ (1 | batch) + (1 | apoe4) + (1 | cogdx) + (1 | msex)

meta_vp <- meta_dat %>% column_to_rownames("batchChannel")
vp_fit <- variancePartition::fitExtractVarPartModel(expr_mat, vp_form, meta_vp)

vp_medians <- as.data.frame(vp_fit) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
  mutate(Factor = factor_labels[Factor]) %>%
  arrange(desc(Median_VarExplained))

# write.csv(vp_medians, file.path(table_output_dir, "dlpfc_proteomics_varpart.csv"), row.names = FALSE)

# Batch confounding check
meta_dat %>%
  group_by(batch, apoe4) %>%
  summarise(n = n(), .groups = "drop") %>%
  pivot_wider(names_from = apoe4, values_from = n, values_fill = 0) %>%
  print(n = Inf)


# Limma ----
## Model fitting & results ----
p_val_threshold <- 0.05

design_mat <- model.matrix(~ apoe4 + age_death_num + msex + cogdx + batch,
                           data = meta_dat) # adjust for sex, age at death, diagnosis and batch
fit <- lmFit(expr_mat, design_mat) %>% eBayes()

res_limma <- topTable(fit, coef = "apoe4APOE4+", number = Inf) %>%
  rownames_to_column(var = "Protein") %>%
  mutate(negLog10FDR = -log10(adj.P.Val),
         `Direction of change` = case_when(
           adj.P.Val < p_val_threshold & logFC > 0 ~ "Upregulated",
           adj.P.Val < p_val_threshold & logFC < 0 ~ "Downregulated",
           TRUE ~ "Not significant"))

sig_res_limma <- res_limma %>%
  filter(adj.P.Val < p_val_threshold)

# write.csv(sig_res_limma, file.path(table_output_dir, "dlpfc_proteomics_limma.csv"), row.names = FALSE)

## Volcano plot ----
p_limma <- ggplot(res_limma, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = `Direction of change`), size = 1.2) +
  ggrepel::geom_text_repel(data = sig_res_limma,
                           aes(label = Protein),
                           size = 3,
                           color = "black") +
  scale_color_manual(values = c("Upregulated" = "indianred",
                                "Downregulated" = "steelblue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "grey") +
  labs(title = sprintf("%d Differentially Expressed dlPFC Proteins by APOE4 Status (n = %d)",
                       nrow(sig_res_limma), ncol(expr_mat)),
       subtitle = "Threshold: adj_p < 0.05",
       x = "log2 fold change",
       y = "-log10 adjusted P",
       color = "Direction of change") +
  theme_classic() +
  theme(legend.position = "bottom")

# ggsave(file.path(pic_output_dir, "dlpfc_limma_volcano.png"), p_limma, width = 8, height = 6)

## PCA ----
pca_deps <- prcomp(t(expr_mat[sig_res_limma$Protein, ]), center = TRUE, scale. = TRUE)

p_pca <- as.data.frame(pca_deps$x[, 1:2]) %>%
  bind_cols(meta_dat %>% select(apoe4, apoe_genotype)) %>%
  ggplot(aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA on %d dlPFC DEPs (n = %d)", nrow(sig_res_limma), ncol(expr_mat)),
       color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "dlpfc_pca.png"), p_pca, width = 7, height = 6)

