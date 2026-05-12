# Set up ----
suppressMessages(library(tidyverse))
library(here)

pic_output_dir <- here("output", "pics", "DiverseCohorts")
table_output_dir <- here("output", "tables", "DiverseCohorts")

proteomic_dat <- read.csv(here("raw", "DiverseCohorts", "normAbundances_post-FP_frontal.csv"), 
                          header = T, strip.white = T, na.strings =  "NA", row.names = 1)
individual_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_individual_metadata.csv"), 
                            header = T, strip.white = T, na.strings = c("", "NA"))
specimen_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_biospecimen_metadata.csv"), 
                            header = T, strip.white = T, na.strings = c("", "NA"))
assay_meta <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_assay_TMTproteomics_metadata.csv"), 
                            header = T, strip.white = T, na.strings = c("", "NA"))

# Preprocessing ----
## Metadata ----
temp_meta <- specimen_meta %>%
  filter(dataGenerationSite == "Emory", tissue == "dorsolateral prefrontal cortex") %>%
  group_by(individualID) %>%
  mutate(status = n() > 1) %>% # Emory is a priority if there's more than 1 entry per ID
  filter(!status | (status & grepl("^emdp", specimenID))) %>%
  ungroup() %>%
  select(individualID, specimenID)

meta_dat <- individual_meta %>%
  select(individualID, sex, race, ageDeath, apoeGenotype, ADoutcome) %>%
  mutate(apoeGenotype = na_if(apoeGenotype, "missing or unknown"),
         ageDeath = na_if(ageDeath, "missing or unknown")) %>%
  filter(!is.na(apoeGenotype), !is.na(ageDeath), ADoutcome %in% c("AD", "Control")) %>%
  inner_join(temp_meta, by = "individualID") %>%
  inner_join(assay_meta[, c("specimenID", "batch")], by = "specimenID") %>%
  mutate(apoe4 = factor(apoeGenotype %in% c("24", "34", "44"), labels = c("APOE4-", "APOE4+")),
         apoeGenotype = factor(apoeGenotype, levels = c("22", "23", "33", "24", "34", "44")),
         age_death_num = as.numeric(if_else(ageDeath == "90+", "90", ageDeath)),
         diag_apoe = factor(paste0(apoe4, " ", ADoutcome),
                            levels = c("APOE4- Control", "APOE4- AD",
                                       "APOE4+ Control", "APOE4+ AD")))

# any(duplicated(meta_dat$individualID))
rm(temp_meta, individual_meta, assay_meta, specimen_meta)


## Proteomic data ----
expr_mat <- proteomic_dat %>% 
  select(all_of(meta_dat$specimenID))
# any(expr_mat == 0, na.rm = TRUE) --> no 0

# median imputation & log transformation
median_impute <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

expr_mat <- expr_mat %>%
  .[rowMeans(is.na(.)) < 0.3, ] %>% # remove proteins with > 30% missing
  apply(1, median_impute) %>% 
  t() 

expr_mat <- log2(expr_mat) 

# boxplot(expr_mat, las = 2, cex.axis = 0.5, outline = FALSE)
# final dim: 8426 proteins x 779 individuals
# all(colnames(expr_mat) == meta_dat$specimenID)

## Save cleaned data----
# saveRDS(expr_mat,  here("raw", "DiverseCohorts", "dlpfc_prot_expr.rds"))
# saveRDS(meta_dat, here("raw", "DiverseCohorts", "dlpfc_prot_meta.rds"))

# Demographic table ----
library(gtsummary)
demographic <- meta_dat %>%
  select(diag_apoe, age_death_num, sex, apoeGenotype) %>%
  tbl_summary(
    by = diag_apoe,
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
  add_p() %>%         
  add_overall(last = TRUE) %>%   
  bold_labels()

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "DC_frontal demographic table.png"))


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

p_varpart <- variancePartition::plotVarPart(vp_fit, main = "Variance Partition: DC-DLPFC Proteomics")
p_varpart$data <- p_varpart$data %>%
  mutate(variable = factor_labels[as.character(variable)])

meta_dat %>%
  group_by(batch, apoe4) %>%
  summarise(n = n()) %>%
  pivot_wider(names_from = apoe4, values_from = n, values_fill = 0) %>%
  print(n = 72)

# write.csv(vp_medians, file.path(table_output_dir, "DC_frontal_varpart_medians.csv"), row.names = FALSE)
# ggsave(file.path(pic_output_dir, "DC_frontal_varpart_violin.png"), p_varpart, width = 7, height = 5)


# Limma analysis ----
## Model: APOE4 effect adjusted for diagnosis, sex, age and batch ----
library(limma)
design_mat <- model.matrix(~ apoe4 + ADoutcome + sex + age_death_num + batch, 
                           data = meta_dat) # adjust for diagnosis, sex, age, and batch

fit <- lmFit(expr_mat, design_mat)
fit <- eBayes(fit)

p_val_threshold <- 0.05

res_limma <- topTable(fit, coef = "apoe4APOE4+", number = Inf) %>%
  rownames_to_column(var = "Protein") %>%
  mutate(negLog10FDR = -log10(adj.P.Val),
         `Direction of change` = case_when(adj.P.Val < p_val_threshold & logFC > 0 ~ "Upregulated",
                                           adj.P.Val < p_val_threshold & logFC < 0 ~ "Downregulated",
                                           TRUE ~ "Not significant"))

sig_res_limma <- res_limma %>%
  filter(adj.P.Val < p_val_threshold) 

# write.csv(res_limma, file.path(table_output_dir, "DC_frontal_all_limma.csv"), row.names = FALSE)
# write.csv(sig_res_limma, file.path(table_output_dir, "DC_frontal_limma.csv"), row.names = FALSE)

## Volcano plot ----
res_limma <- res_limma %>% filter(Protein != "APOE4|APOE4") # APOE4 is top DEP and distort volcano plot
sig_res_limma <- sig_res_limma %>% filter(Protein != "APOE4|APOE4")

p_limma <- ggplot(res_limma, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = `Direction of change`), alpha = 0.6, size = 1.2) +
  ggrepel::geom_text_repel(data = sig_res_limma,
                           aes(label = Protein),
                           size = 3,
                           color = "black") +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = 0, 
             linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(p_val_threshold), 
             linetype = "dashed", color = "grey") +
  labs(title = sprintf("%d Differentially Expressed dlPCF Proteins by APOE4 Status (n = %d)", 
                       nrow(sig_res_limma) + 1, nrow(meta_dat)),
       subtitle = "Thresholds: adj_p < 0.05 | APOE4 not shown on plot",
       x = "log2 fold change",
       y = "-log10 adjusted P",
       color = "Direction of change"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

# ggsave(file.path(pic_output_dir, "DC_frontal_limma_volcano.png"), p_limma, width = 8, height = 6)
