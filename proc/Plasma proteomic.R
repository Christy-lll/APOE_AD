# Set up ----
suppressMessages(library(tidyverse))
library(here)

pic_output_dir <- here("output", "pics", "ROSMAP")
table_output_dir <- here("output", "tables", "ROSMAP")

proteomic_dat <- read.csv(here("raw", "ROSMAP", "OhNM2025_ROSMAP_plasma_Soma7k_protein_level_ANML_log10.csv"), 
                          header = T, strip.white = T, na.strings = "")
sample_meta <- read.csv(here("raw", "ROSMAP", "OhNM2025_ROSMAP_plasma_Soma7k_sample_metadata_1_.csv"), 
                        header = T, strip.white = T, na.strings = "") 
protein_meta <- read.csv(here("raw", "ROSMAP", "OhNM2025_ROSMAP_plasma_Soma7k_protein_metadata_1_.csv"), 
                         header = T, strip.white = T, na.strings = "")

# Preprocessing ----
## Map seqIDs to Uniprot IDs ----
colnames(proteomic_dat) <- sub("^X(\\d)", "\\1", colnames(proteomic_dat))
colnames(proteomic_dat) <- gsub("\\.", "-", colnames(proteomic_dat))
protein_vars <- !colnames(proteomic_dat) %in% "projid_visit"
colnames(proteomic_dat)[protein_vars] <- protein_meta$UniProt[match(colnames(proteomic_dat)[protein_vars], protein_meta$SeqId)]
proteomic_dat <- proteomic_dat[, !is.na(colnames(proteomic_dat)) &
                                 !duplicated(colnames(proteomic_dat)) &
                                 colnames(proteomic_dat) != ""] # from 7289 to 6402 proteins

protein_meta_clean <- protein_meta %>%
  dplyr::select(UniProt, EntrezGeneSymbol, TargetFullName) %>%
  distinct(UniProt, .keep_all = TRUE)

## Merge data ----
dat <- merge(proteomic_dat, sample_meta, by = "projid_visit")

## Handle duplication in projid ----
# table(table(dat$projid)) --> 47 with duplication (2-3)
dat <- dat %>%
  filter(!is.na(projid)) %>% # drop 3 samples with no projid
  group_by(projid) %>%
  slice_min(order_by = Visit, n = 1, with_ties = FALSE) %>% # keep 1st visit if having more than 1
  ungroup() 

## Recode diagnosis, APOE4 status & sex ----
dat <- dat %>%
  filter(!is.na(apoe_genotype)) %>% # drop 65 samples with no apoe_genotype
  filter(Diagnosis %in% c("AD", "MCI", "NCI")) %>% # drop 30 samples 
  mutate(Diagnosis = factor(Diagnosis, 
                            levels = c("NCI", "MCI", "AD")), 
         apoe4 = factor(grepl("4", apoe_genotype),
                        levels = c(FALSE, TRUE),
                        labels = c("APOE4-", "APOE4+")),
         apoe_genotype = factor(apoe_genotype, 
                                levels = c("22", "23", "33", "24", "34", "44")),
         diag_apoe = factor(paste0(apoe4, " ", Diagnosis),
                            levels = c("APOE4- NCI", "APOE4- MCI", "APOE4- AD", 
                                       "APOE4+ NCI", "APOE4+ MCI", "APOE4+ AD")),
         msex = factor(msex,
                       levels = c("0", "1"),
                       labels = c("Female", "Male")))

## Define variable groups ----
protein_vars <- colnames(proteomic_dat)[colnames(proteomic_dat) != "projid_visit"]
meta_vars <- setdiff(colnames(dat), protein_vars)

# final n = 793, proteins = 6402

# Demographic & PCA ----
## Demographic table ----
library(gtsummary)
demographic <- dat %>%
  select(diag_apoe, age_at_visit, msex, educ, apoe_genotype) %>%
  tbl_summary(
    by = diag_apoe,
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age_at_visit ~ "Age at Visit",
      msex ~ "Sex",
      educ ~ "Education (years)",
      apoe_genotype ~ "APOE genotype"
    ),
    missing = "no"
  ) %>%
  add_p() %>%         
  add_overall(last = TRUE) %>%   
  bold_labels()

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "plasma demographic table.png"))

## PCA ----
pca_all <- prcomp(dat[, protein_vars], center = TRUE, scale. = TRUE)

pca_all_df <- as.data.frame(pca_all$x[, 1:2]) %>%
  bind_cols(dat %>% select(projid_visit, diag_apoe, apoe4, apoe_genotype)) 

p_pca_all <- ggplot(pca_all_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = "PCA with 6402 plasma proteins",
       x = "PC1",
       y = "PC2",
       color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "plasma_allproteins_pca.png"), p_pca_all, width = 7, height = 6)

# Limma (adjusted for diagnosis) ----
## Model fitting ----
expr_mat <- as.matrix(dat[, protein_vars]) # no missing values anyNA(expr_mat)
rownames(expr_mat) <- dat$projid_visit
expr_mat <- t(expr_mat)
design_mat <- model.matrix(~ apoe4 + Diagnosis + age_at_visit + msex + educ, data = dat) # adjust educ?

library(limma)
fit <- lmFit(expr_mat, design_mat)
fit <- eBayes(fit)

## Result ----
p_val_threshold <- 0.05
logfc_threshold <- 0.1

res_limma <- topTable(fit, coef = "apoe4APOE4+", number = Inf) %>%
  rownames_to_column(var = "UniprotID") %>%
  mutate(negLog10FDR = -log10(adj.P.Val),
         `Direction of change` = case_when(adj.P.Val < p_val_threshold & logFC > logfc_threshold ~ "Upregulated",
                                           adj.P.Val < p_val_threshold & logFC < -logfc_threshold ~ "Downregulated",
                                           TRUE ~ "Not significant"))

sig_res_limma <- res_limma %>%
  filter(adj.P.Val < p_val_threshold, abs(logFC) > logfc_threshold) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  select(UniprotID, EntrezGeneSymbol, TargetFullName, everything())

# write.csv(sig_res_limma, file.path(table_output_dir, "plasma_limma.csv"), row.names = FALSE)

## Volcano plot ----
p_limma <- ggplot(res_limma, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = `Direction of change`), alpha = 0.6, size = 1.2) +
  ggrepel::geom_text_repel(data = sig_res_limma,
                  aes(label = EntrezGeneSymbol),
                  size = 3,
                  color = "black") +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
             linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(p_val_threshold), 
             linetype = "dashed", color = "grey") +
  labs(title = "13 Differentially Expressed Plasma Proteins by APOE4 Status (limma)",
    subtitle = "Thresholds: adj_p < 0.05 & |log2 FC| > 0.1",
    x = "log2 fold change",
    y = "-log10 adjusted P",
    color = "Direction of change"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

# ggsave(file.path(pic_output_dir, "plasma_limma_volcano.png"), p_limma, width = 8, height = 6)

## PCA check ----
limma_sig_proteins <- sig_res_limma$UniprotID
pca_limma <- prcomp(dat[, limma_sig_proteins], center = TRUE, scale. = TRUE)
pca_limma_df <- as.data.frame(pca_limma$x[, 1:2]) %>%
  bind_cols(dat %>% select(projid_visit, diag_apoe, apoe4, apoe_genotype)) 

p_pca_limma <- ggplot(pca_limma_df, aes(x = PC1, y = PC2, color = apoe_genotype)) +
  geom_point() +
  labs(title = sprintf("PCA with %d differentially expressed plasma proteins (limma)", length(limma_sig_proteins)),
       x = "PC1",
       y = "PC2",
       color = "APOE genotype") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "plasma_limma_pca.png"), p_pca_limma, width = 7, height = 6)


# Machine Learning ----
## Feature selection (on NCI) ----
mi_dat <- dat %>% 
  filter(Diagnosis == "NCI") %>%
  select(apoe4, all_of(protein_vars))

mi_scores_all <- FSelectorRcpp::information_gain(apoe4 ~ ., data = mi_dat, type = "infogain") %>%
  arrange(desc(importance)) %>%
  filter(importance > 0) %>%
  left_join(protein_meta_clean, by = c("attributes" = "UniProt")) %>%
  mutate(selected = importance > 0.08)

mi_scores <- mi_scores_all %>%
  filter(importance > 0.08) %>% # threshold selected by sensitivity analysis --> 14 proteins
  select(attributes, EntrezGeneSymbol, TargetFullName, importance)

mi_sig_proteins <- mi_scores$attributes
mi_sig_proteins_name <- mi_scores$EntrezGeneSymbol
# intersect(mi_sig_proteins, sig_res_limma$UniprotID)
# write.csv(mi_scores, file.path(table_output_dir, "plasma_MI.csv"), row.names = FALSE)

# Barplot 
p_mi_bar <- ggplot(mi_scores_all, aes(x = importance, 
                                      y = forcats::fct_reorder(EntrezGeneSymbol, importance),
                                      fill = selected)) +
  geom_col(alpha = 0.9) +
  geom_vline(xintercept = 0.08, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey75"),
                    name = "Selected for ML") +
  labs(title = sprintf("%d Proteins with MI > 0", nrow(mi_scores_all)),
       x = "Mutual Information (MI)",
       y = NULL) +
  theme_classic() 

# ggsave(file.path(pic_output_dir, "plasma_MI_barplot.png"), p_mi_bar, width = 10, height = 8)

# PCA check 
pca_MI <- prcomp(dat[, mi_sig_proteins], center = TRUE, scale. = TRUE)
pca_MI_df <- as.data.frame(pca_MI$x[, 1:2]) %>%
  bind_cols(dat %>% select(projid_visit, diag_apoe, apoe4, apoe_genotype)) 

p_pca_MI <- ggplot(pca_MI_df, aes(x = PC1, y = PC2, color = apoe_genotype)) +
  geom_point() +
  labs(title = sprintf("PCA with %d MI-selected plasma proteins", length(mi_sig_proteins)),
       x = "PC1",
       y = "PC2",
       color = "APOE genotype") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "plasma_MI_pca.png"), p_pca_MI, width = 7, height = 6)

## Model building (on NCI) ----
library(caret)
train_dat <- mi_dat %>%
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

cv_ctrl <- trainControl(
  method = "repeatedcv", 
  number = 5, 
  repeats = 5,
  classProbs = TRUE, 
  summaryFunction = twoClassSummary,
  sampling = "smote",
  savePredictions = "final"
)

set.seed(16)
model_fit <- train(
  class ~ .,
  data = train_dat,
  method = "rf",
  metric = "ROC",
  trControl = cv_ctrl,
  tuneLength = 10,
  importance = TRUE
)

## Test on MCI ----
test_MCI <- dat %>%
  filter(Diagnosis == "MCI") %>%
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

mci_class <- predict(model_fit, newdata = test_MCI)
mci_probs <- predict(model_fit, newdata = test_MCI, type = "prob")[, "pos"]
mci_cm <- confusionMatrix(mci_class, test_MCI$class, positive = "pos")
mci_roc <- pROC::roc(response = test_MCI$class, predictor = mci_probs, levels = c("neg", "pos"))

mci_metrics <- tibble(
  Model = "Random Forest (APOE4- MCI vs. APOE4+ MCI)",
  N_train = nrow(train_dat),
  N_test = nrow(test_MCI),
  N_features = length(mi_sig_proteins),
  Sensitivity = mci_cm$byClass["Sensitivity"],
  Specificity = mci_cm$byClass["Specificity"],
  PPV = mci_cm$byClass["Pos Pred Value"],
  NPV = mci_cm$byClass["Neg Pred Value"],
  Train_AUC = model_fit$results %>% 
    filter(mtry == model_fit$bestTune$mtry) %>% pull(ROC),
  Test_AUC = as.numeric(mci_roc$auc),
  AUC_gap = Train_AUC - Test_AUC
)

## Test on AD ----
test_AD <- dat %>%
  filter(Diagnosis == "AD") %>%
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

ad_class <- predict(model_fit, newdata = test_AD)
ad_probs <- predict(model_fit, newdata = test_AD, type = "prob")[, "pos"]
ad_cm <- confusionMatrix(ad_class, test_AD$class, positive = "pos")
ad_roc <- pROC::roc(response = test_AD$class, predictor = ad_probs, levels = c("neg", "pos"))

ad_metrics <- tibble(
  Model = "Random Forest (APOE4- AD vs. APOE4+ AD)",
  N_train = nrow(train_dat),
  N_test = nrow(test_AD),
  N_features = length(mi_sig_proteins),
  Sensitivity = ad_cm$byClass["Sensitivity"],
  Specificity = ad_cm$byClass["Specificity"],
  PPV = ad_cm$byClass["Pos Pred Value"],
  NPV = ad_cm$byClass["Neg Pred Value"],
  Train_AUC = model_fit$results %>% 
    filter(mtry == model_fit$bestTune$mtry) %>% pull(ROC),
  Test_AUC = as.numeric(ad_roc$auc),
  AUC_gap = Train_AUC - Test_AUC
)

## Extension (train on female; test on male) ----
female_dat <- dat %>%
  filter(msex == "Female") %>% 
  select(apoe4, all_of(protein_vars)) %>%
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

set.seed(16)
sex_model_fit <- train(
  class ~ .,
  data = female_dat,
  method = "rf",
  metric = "ROC",
  trControl = cv_ctrl,
  tuneLength = 10,
  importance = TRUE
)

male_dat <- dat %>%
  filter(msex == "Male") %>% 
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

male_class <- predict(sex_model_fit, newdata = male_dat)
male_probs <- predict(sex_model_fit, newdata = male_dat, type = "prob")[, "pos"]
male_cm <- confusionMatrix(male_class, male_dat$class, positive = "pos")
male_roc <- pROC::roc(response = male_dat$class, predictor = male_probs, levels = c("neg", "pos"))

male_metrics <- tibble(
  Model = "Random Forest (APOE4- Male vs. APOE4+ Male)",
  N_train = nrow(female_dat),
  N_test = nrow(male_dat),
  N_features = length(mi_sig_proteins),
  Sensitivity = male_cm$byClass["Sensitivity"],
  Specificity = male_cm$byClass["Specificity"],
  PPV = male_cm$byClass["Pos Pred Value"],
  NPV = male_cm$byClass["Neg Pred Value"],
  Train_AUC = sex_model_fit$results %>%
    filter(mtry == sex_model_fit$bestTune$mtry) %>% pull(ROC),
  Test_AUC = as.numeric(male_roc$auc),
  AUC_gap = Train_AUC - Test_AUC
)

## Extension (train on ROSMAP; test on ADNI) ---- 
ROSMAP_dat <- dat %>%
  select(apoe4, all_of(protein_vars)) %>%
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins))) %>%
  rename(any_of(setNames(mi_sig_proteins, mi_sig_proteins_name))) %>%
  select(-ZW10) # ADNI data doesn't have ZW10

set.seed(16)
ROSMAP_model_fit <- train(
  class ~ .,
  data = ROSMAP_dat,
  method = "rf",
  metric = "ROC",
  trControl = cv_ctrl,
  tuneLength = 10,
  importance = TRUE
)

ADNI_dat <- read.csv(here("raw", "ADNI_log_norm_dat.csv"), header = T)

ADNI_dat <- ADNI_dat %>%
  transmute(class = factor(ifelse(apoe == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(any_of(mi_sig_proteins_name)))

ADNI_class <- predict(ROSMAP_model_fit, newdata = ADNI_dat)
ADNI_probs <- predict(ROSMAP_model_fit, newdata = ADNI_dat, type = "prob")[, "pos"]
ADNI_cm <- confusionMatrix(ADNI_class, ADNI_dat$class, positive = "pos")
ADNI_roc <- pROC::roc(response = ADNI_dat$class, predictor = ADNI_probs, levels = c("neg", "pos"))

ADNI_metrics <- tibble(
  Model = "Random Forest (External: ADNI)",
  N_train = nrow(ROSMAP_dat),
  N_test = nrow(ADNI_dat),
  N_features = ncol(ADNI_dat) - 1,
  Sensitivity = ADNI_cm$byClass["Sensitivity"],
  Specificity = ADNI_cm$byClass["Specificity"],
  PPV = ADNI_cm$byClass["Pos Pred Value"],
  NPV = ADNI_cm$byClass["Neg Pred Value"],
  Train_AUC = ROSMAP_model_fit$results %>% 
    filter(mtry == ROSMAP_model_fit$bestTune$mtry) %>% pull(ROC),
  Test_AUC = as.numeric(ADNI_roc$auc),
  AUC_gap = Train_AUC - Test_AUC
  )

## ML Metrics ----
all_metrics <- bind_rows(mci_metrics, ad_metrics, male_metrics, ADNI_metrics)

# write.csv(all_metrics, file.path(table_output_dir, "plasma_ML_metrics.csv"), row.names = FALSE)


# Mediation analysis ----
# Exposure: APOE4 status (binary); Mediator: Proteins (continuous); Outcome: AD Diagnosis (binary)

## Feature screening with limma (APOE4 --> proteins) ----
expr_mat <- as.matrix(dat[, protein_vars]) 
rownames(expr_mat) <- dat$projid_visit
expr_mat <- t(expr_mat)
design_mat_med <- model.matrix(~ apoe4 + age_at_visit + msex + educ, data = dat) # no adjustment for diagnosis

library(limma)
fit_med <- lmFit(expr_mat, design_mat_med)
fit_med <- eBayes(fit_med)

res_limma_med <- topTable(fit_med, coef = "apoe4APOE4+", number = Inf) %>%
  rownames_to_column(var = "UniprotID") %>%
  mutate(negLog10FDR = -log10(adj.P.Val)) %>%
  filter(adj.P.Val < 0.05) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, TargetFullName, everything()) # 53 significant proteins

## DEPs & diag_apoe heatmap ----
heatmap_mat <- dat %>%
  dplyr::select(diag_apoe, all_of(med_proteins)) %>%
  group_by(diag_apoe) %>%
  summarise(across(everything(), mean, na.rm = TRUE)) %>%
  column_to_rownames("diag_apoe") %>%
  t()

rownames(heatmap_mat) <- protein_meta_clean$EntrezGeneSymbol[
  match(rownames(heatmap_mat), protein_meta_clean$UniProt)
]

p_heatmap <- pheatmap::pheatmap(heatmap_mat,
                                scale = "row",         
                                cluster_cols = FALSE,   
                                cluster_rows = TRUE,
                                show_colnames = TRUE,
                                fontsize_row = 9,
                                main = "Z-scored mean expression of DEPs across diagnosis and APOE4 subgroups")

# ggsave(file.path(pic_output_dir, "plasma_limma_heatmap.png"), plot = p_heatmap, width = 8, height = 8)


## Mediation model: Does the protein mediate APOE4's effect on AD ----
med_dat <- dat %>% 
  filter(Diagnosis %in% c("NCI", "AD")) %>%
  mutate(diagnosis_bin = as.integer(Diagnosis == "AD")) # 0 = NCI, 1 = AD

med_proteins <- res_limma_med$UniprotID

library(mediation)
run_mediation <- function(protein) {
  df <- med_dat %>% dplyr::select(diagnosis_bin, apoe4, age_at_visit, msex, educ, mediator = all_of(protein))
  
  med.fit <- lm(mediator ~ apoe4 + age_at_visit + msex + educ, data = df)
  out.fit <- glm(diagnosis_bin ~ mediator + apoe4 + age_at_visit + msex + educ, data = df, family = binomial(link = "logit"))
  
  med.out <- mediate(med.fit, out.fit, treat = "apoe4", mediator = "mediator",
                     control.value = "APOE4-", treat.value = "APOE4+", 
                     boot = TRUE, sims = 1500)
  
  tibble(UniprotID = protein,
         ACME_est = med.out$d.avg, # Average casual mediation effect: APOE4's effect on AD through this protein
         ACME_ci_low = med.out$d.avg.ci[1],
         ACME_ci_up = med.out$d.avg.ci[2],
         ACME_p = med.out$d.avg.p,
         ADE_est = med.out$z.avg, # Average direct effect: The effect of APOE4 on AD not explained by this protein
         ADE_p = med.out$z.avg.p,
         prop_med = med.out$n.avg, # Proportion of total effect mediated by this protein
         prop_ci_low = med.out$n.avg.ci[1],
         prop_ci_up = med.out$n.avg.ci[2],
         total_effect = med.out$tau.coef, # ACME + ADE
         total_effect_p = med.out$tau.p)
}

library(furrr)
plan(multisession, workers = parallel::detectCores() - 1)

set.seed(16)
res_mediation <- map_dfr(med_proteins, run_mediation) %>%
  mutate(ACME_fdr = p.adjust(ACME_p, method = "BH")) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, TargetFullName, 
                ACME_est, ACME_ci_low, ACME_ci_up, ACME_p, ACME_fdr, ADE_est, ADE_p, 
                prop_med, prop_ci_low, prop_ci_up, total_effect, total_effect_p) %>%
  arrange(ACME_fdr)

plan(sequential)

## Result ----
sig_res_mediation <- res_mediation %>%
  filter(ACME_fdr < 0.05) # 12 proteins with significant ACME FDR

# write.csv(sig_res_mediation, file.path(table_output_dir, "plasma_mediation_APOE_prot_AD.csv"), row.names = FALSE)

# Forest plot
library(patchwork)
p_ACME <- ggplot(sig_res_mediation, aes(x = ACME_est, y = forcats::fct_reorder(EntrezGeneSymbol, ACME_est))) +
  geom_point(size = 2) + 
  geom_errorbarh(aes(xmin = ACME_ci_low, xmax = ACME_ci_up), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "Average causal mediation effect (ACME)", y = NULL) +
  theme_classic()

p_PROP <- ggplot(sig_res_mediation, aes(x = prop_med * 100, y = forcats::fct_reorder(EntrezGeneSymbol, ACME_est))) +
  geom_point(size = 2) + 
  geom_errorbarh(aes(xmin = prop_ci_low * 100, xmax = prop_ci_up * 100), height = 0.2) +
  labs(x = "Proportion mediated (%)", y = NULL) +
  theme_classic() 

p_forest <- p_ACME + p_PROP +
  plot_annotation(title = "Mediation Analysis: APOE4 → Protein → AD Diagnosis",
                  subtitle = "12 proteins with ACME FDR < 0.05 | 95% bootstrap CI")

# ggsave(file.path(pic_output_dir, "plasma_mediation_APOE_prot_AD.png"), p_forest, width = 9, height = 6)

# Volcano plot 
res_limma_med_volcano <- topTable(fit_med, coef = "apoe4APOE4+", number = Inf) %>%
  rownames_to_column(var = "UniprotID") %>%
  mutate(negLog10FDR = -log10(adj.P.Val),
         `Direction of change` = case_when(adj.P.Val < 0.05 & logFC > 0 ~ "Upregulated",
                                           adj.P.Val < 0.05 & logFC < 0 ~ "Downregulated",
                                           TRUE ~ "Not significant"), 
         mediation_sig = UniprotID %in% sig_res_mediation$UniprotID) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt"))

p_limma_med <- ggplot(res_limma_med_volcano, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = `Direction of change`), alpha = 0.6, size = 1.2) +
  ggrepel::geom_text_repel(data = filter(res_limma_med_volcano, mediation_sig),
                           aes(label = EntrezGeneSymbol),
                           size = 3, color = "black") +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not significant" = "grey")) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(title = "APOE4-associated plasma proteins (unadjusted for diagnosis)",
       subtitle = "Labelled: proteins with significant mediation effect (ACME FDR < 0.05)",
       x = "log2 fold change",
       y = "-log10 adjusted P",
       color = " ") +
  theme_classic() +
  theme(legend.position = "bottom")

# ggsave(file.path(pic_output_dir, "plasma_mediation_volcano.png"), p_limma_med, width = 8, height = 7)

## Sensitivity check ----
# df_check <- med_dat %>% dplyr::select(diagnosis_bin, apoe4, age_at_visit, msex, educ, [UNIPROT])
# med.fit <- lm([UNIPROT] ~ apoe4 + age_at_visit + msex + educ, data = df_check)
# out.fit <- glm(diagnosis_bin ~ [UNIPROT] + apoe4 + age_at_visit + msex + educ,
#                data = df_check, family = binomial(link = "probit"))
# med.out <- mediate(med.fit, out.fit, treat = "apoe4", mediator = "[UNIPROT]",
#                    control.value = "APOE4-", treat.value = "APOE4+",
#                    boot = TRUE, sims = 1000)
# sens.out <- medsens(med.out, rho.by = 0.05, effect.type = "indirect")
# summary(sens.out)


