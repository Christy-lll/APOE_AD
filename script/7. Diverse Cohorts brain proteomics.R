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

# Preprocessing ----
## Metadata ----
meta_base <- individual_meta %>%
  select(individualID, sex, race, ageDeath, apoeGenotype, ADoutcome, PMI) %>%
  mutate(apoeGenotype = na_if(apoeGenotype, "missing or unknown"),
         ageDeath = na_if(ageDeath, "missing or unknown")) %>%
  filter(!is.na(apoeGenotype), !is.na(ageDeath), ADoutcome %in% c("AD", "Control")) %>%
  mutate(apoe4 = factor(apoeGenotype %in% c("24", "34", "44"), labels = c("APOE ε4-", "APOE ε4+")),
         apoeGenotype = factor(apoeGenotype, levels = c("22", "23", "33", "24", "34", "44")),
         age_death_num = as.numeric(if_else(ageDeath == "90+", "90", ageDeath)),
         ADoutcome = factor(ADoutcome, levels = c("Control", "AD")),
         PMI = as.numeric(PMI))

non_control <- assay_meta %>%
  filter(isAssayControl == FALSE) %>%
  pull(specimenID)

dlpfc_spec <- specimen_meta %>%
  filter(tissue == "dorsolateral prefrontal cortex",
         specimenID %in% non_control,
         dataGenerationSite == "Emory") %>%
  group_by(individualID) %>%
  mutate(is_dup = n() > 1) %>%
  filter(!is_dup | grepl("^emdp", specimenID)) %>%
  ungroup() %>%
  select(individualID, specimenID)

meta_dat <- meta_base %>%
  inner_join(dlpfc_spec, by = "individualID") %>%
  inner_join(assay_meta %>% select(specimenID, batch), by = "specimenID")
# any(duplicated(meta_dat$individualID))

## Proteomic data ----
median_impute <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

expr_mat <- dlpfc_prot_raw %>%
  filter(!rownames(.) %in% c("APOE4|APOE4")) %>%
  select(all_of(meta_dat$specimenID)) %>%
  mutate(across(where(is.numeric), ~ na_if(.x, 0))) %>% # code 0 as NA
  .[rowMeans(is.na(.)) < 0.3, ] %>% # remove proteins with >30% missing
  apply(1, median_impute) %>% # median imputation
  t() %>%
  log2() # log transformation

## Save data ----
# saveRDS(list(expr = expr_mat, meta = meta_dat), here("raw", "DiverseCohorts", "DC-proteomics.rds"))

rm(individual_meta, assay_meta, specimen_meta, meta_base, non_control, dlpfc_spec, dlpfc_prot_raw)


# Demographic ----
demographic <- meta_dat %>%
  select(ADoutcome, age_death_num, sex, apoeGenotype, PMI, race) %>%
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
      apoeGenotype ~ "APOE genotype",
      PMI ~ "Post mortem interval (hours)",
      race ~ "Race"
    ),
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels()

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "dlpfc proteomics demographic table.png"))


# Variance Partition Analysis ----
factor_labels <- c(apoe4 = "APOE4 Status",
                   ADoutcome = "Diagnosis",
                   batch = "Batch",
                   sex = "Sex",
                   Residuals = "Residuals")

vp_form <- ~ (1 | batch) + (1 | apoe4) + (1 | ADoutcome) + (1 | sex)

meta_vp <- meta_dat %>% column_to_rownames("specimenID")
vp_fit <- variancePartition::fitExtractVarPartModel(expr_mat, vp_form, meta_vp)

vp_medians <- as.data.frame(vp_fit) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
  mutate(Factor = factor_labels[Factor]) %>%
  arrange(desc(Median_VarExplained))


# Limma ----
## Model fitting & results ----
p_val_threshold <- 0.05

design_mat <- model.matrix(~ apoe4 + age_death_num + sex + ADoutcome + batch,
                           data = meta_dat) # adjust for age, sex, diagnosis and batch
fit <- lmFit(expr_mat, design_mat) %>% eBayes()

res_limma <- topTable(fit, coef = "apoe4APOE ε4+", number = Inf) %>%
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
n_up <- sum(sig_res_limma$`Direction of change` == "Upregulated")
n_down <- sum(sig_res_limma$`Direction of change` == "Downregulated")

p_limma_volcano <- ggplot(res_limma, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = `Direction of change`), size = 1.2) +
  ggrepel::geom_text_repel(data = sig_res_limma,
                           aes(label = Protein),
                           size = 3, color = "black") +
  scale_color_manual(values = c("Upregulated" = "indianred",
                                "Downregulated" = "steelblue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "grey") +
  labs(title = sprintf("%d dlPFC DAPs in AMP-AD Diverse Cohorts Study",
                       nrow(sig_res_limma)),
       subtitle = sprintf("Up: %d  Down: %d", n_up, n_down),
       x = "log2 fold change",
       y = "-log10 FDR",
       color = NULL) +
  theme_classic() +
  theme(legend.position = "bottom")

# ggsave(file.path(pic_output_dir, "dlpfc_proteomics_limma_volcano.tiff"), p_limma_volcano, width = 7, height = 6, dpi = 300, compression = "lzw")


# Mutual Information Analysis ----
## Filter Controls ----
ctrl_idx <- meta_dat$ADoutcome == "Control"
meta_ctrl <- meta_dat[ctrl_idx, ]
expr_ctrl <- expr_mat[, ctrl_idx]

## Batch correction ----
design_bc <- model.matrix(~ apoe4 + sex + age_death_num, data = meta_ctrl)
expr_ctrl_bc <- removeBatchEffect(x = expr_ctrl, batch = meta_ctrl$batch, design = design_bc)

## Compute MI ----
mi_dat <- expr_ctrl_bc %>% t() %>% as.data.frame() %>% mutate(apoe4 = meta_ctrl$apoe4)

mi_scores_all <- FSelectorRcpp::information_gain(apoe4 ~ ., data = mi_dat, type = "infogain") %>%
  arrange(desc(importance)) %>%
  dplyr::rename(Protein = attributes)

mi_scores <- mi_scores_all %>%
  filter(importance > 0)

## Barplot ----
mi_top30 <- mi_scores_all %>% slice_max(importance, n = 30)

p_mi_bar <- ggplot(mi_top30, aes(x = importance, y = fct_reorder(Protein, importance))) +
  geom_col(fill = "steelblue", alpha = 0.9) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.6) +
  labs(title = "AMP-AD Diverse Cohorts Study dlPFC top 30 MI-selected proteins",
       x = "MI score", y = NULL) +
  theme_classic()

# ggsave(file.path(pic_output_dir, "dlpfc_proteomics_MI_barplot.tiff"), p_mi_bar, width = 7, height = 6, dpi = 300, compression = "lzw")

