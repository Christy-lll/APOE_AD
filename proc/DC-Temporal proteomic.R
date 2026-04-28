# Set up ----
suppressMessages(library(tidyverse))
library(here)

pic_output_dir <- here("output", "pics", "DiverseCohorts")
table_output_dir <- here("output", "tables", "DiverseCohorts")

proteomic_dat <- read.csv(here("raw", "DiverseCohorts", "normAbundances_post-FP_temporal.csv"), 
                          header = T, strip.white = T, na.strings =  "NA", row.names = 1)
individual_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_individual_metadata.csv"), 
                            header = T, strip.white = T, na.strings = c("", "NA"))
specimen_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_biospecimen_metadata.csv"), 
                          header = T, strip.white = T, na.strings = c("", "NA"))
assay_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_assay_TMTproteomics_metadata.csv"), 
                       header = T, strip.white = T, na.strings = c("", "NA"))

# Preprocessing ----
## Metadata ----
specimen_meta <- specimen_meta %>%
  filter(specimenID %in% colnames(proteomic_dat),
         specimenID %in% assay_meta$specimenID,
         individualID %in% individual_meta$individualID) %>%
  select(individualID, specimenID, tissue, assay) %>%
  left_join(assay_meta %>% select(specimenID, batch), by = "specimenID") %>%
  filter(n() == 1, .by = individualID) # keep only 1 sample per individual

meta_dat <- assay_meta %>%
  filter(specimenID %in% specimen_meta$specimenID,
         isAssayControl == FALSE) %>%
  select(specimenID, platform, TMTType) %>%
  inner_join(specimen_meta, by = "specimenID") %>%
  inner_join(individual_meta, by = "individualID") %>%
  filter(ADoutcome %in% c("AD", "Control"),
         ageDeath != "missing or unknown",
         apoeGenotype != "missing or unknown") %>%
  mutate(apoe4 = factor(grepl("4", apoeGenotype), labels = c("APOE4-", "APOE4+")),
         apoeGenotype = factor(apoeGenotype, levels = c("22", "23", "33", "24", "34", "44")),
         age_death_90plus = as.integer(ageDeath == "90+"),
         age_death_num = as.numeric(if_else(ageDeath == "90+", "90", ageDeath)),
         diag_apoe = factor(paste0(apoe4, " ", ADoutcome),
                            levels = c("APOE4- Control", "APOE4- AD",
                                       "APOE4+ Control", "APOE4+ AD"))) 

## Proteomic data ----
# filter sample & code 0 as NA
expr_mat <- proteomic_dat %>% 
  select(all_of(meta_dat$specimenID)) %>%
  mutate(across(where(is.numeric), ~na_if(.x, 0))) 

# median imputation 
median_impute <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

expr_mat <- expr_mat %>%
  .[rowMeans(is.na(.)) < 0.3, ] %>% # remove proteins with > 30% missing
  apply(1, median_impute) %>% 
  t() 

# log-transform
expr_mat <- log2(expr_mat) 

# boxplot(expr_mat, las = 2, cex.axis = 0.5, outline = FALSE)
# final dim: 8628 proteins x 208 individuals
# all(colnames(expr_mat) == meta_dat$specimenID)


# Demographic table ----
library(gtsummary)
demographic <- meta_dat %>%
  select(diag_apoe, age_death_num, age_death_90plus, sex, apoeGenotype) %>%
  tbl_summary(
    by = diag_apoe,
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age_death_num ~ "Age at Death (max 90)",
      age_death_90plus ~ "Age over 90",
      sex ~ "Sex",
      apoeGenotype ~ "APOE genotype"
    ),
    missing = "no"
  ) %>%
  add_p() %>%         
  add_overall(last = TRUE) %>%   
  bold_labels()

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "DC_temporal demographic table.png"))


# Variance Partition Analysis ----
meta_vp <- meta_dat %>% 
  column_to_rownames("specimenID")

vp_form <- ~ (1 | batch) + (1 | apoe4) + (1 | ADoutcome) + (1 | sex)
vp_fit <- variancePartition::fitExtractVarPartModel(expr_mat, vp_form, meta_vp)

factor_labels <- c(apoe4 = "APOE4 status",
                   ADoutcome = "Diagnosis",
                   batch = "Batch",
                   sex = "Sex",
                   Residuals = "Residuals")

vp_medians <- as.data.frame(vp_fit) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
  mutate(Factor = factor_labels[Factor]) %>%
  arrange(desc(Median_VarExplained)) 

p_varpart <- variancePartition::plotVarPart(vp_fit, main = "Variance Partition: DC-STG Proteomics")
p_varpart$data <- p_varpart$data %>%
  mutate(variable = factor_labels[as.character(variable)])


# Batch-diagnosis/apoe4 confounding check
meta_dat %>%
  group_by(batch, diag_apoe) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = diag_apoe, values_from = n, values_fill = 0) %>%
  print(n = 19)

# write.csv(vp_medians, file.path(table_output_dir, "DC_temporal_varpart_medians.csv"), row.names = FALSE)
# ggsave(file.path(pic_output_dir, "DC_temporal_varpart_violin.png"), p_varpart, width = 7, height = 5)


# Limma analysis ----
## Model: APOE4 effect adjusted for diagnosis, sex, age and batch ----
library(limma)
design_mat <- model.matrix(~ apoe4 + ADoutcome + sex + age_death_num + age_death_90plus + batch , 
                           data = meta_dat) # adjust for diagnosis, sex, age, and batch

fit <- lmFit(expr_mat, design_mat)
fit <- eBayes(fit)

p_val_threshold <- 0.05
logfc_threshold <- 0.13

res_limma <- topTable(fit, coef = "apoe4APOE4+", number = Inf) %>%
  rownames_to_column(var = "Protein") %>%
  mutate(negLog10FDR = -log10(adj.P.Val),
         `Direction of change` = case_when(adj.P.Val < p_val_threshold & logFC > logfc_threshold ~ "Upregulated",
                                           adj.P.Val < p_val_threshold & logFC < -logfc_threshold ~ "Downregulated",
                                           TRUE ~ "Not significant"))

# No significant results
# write.csv(res_limma, file.path(table_output_dir, "DC_temporal_all_limma.csv"), row.names = FALSE)
