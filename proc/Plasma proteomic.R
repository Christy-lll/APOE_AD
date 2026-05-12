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
design_mat <- model.matrix(~ apoe4 + Diagnosis + age_at_visit + msex, data = dat) # adjust educ?

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


# AD/APOE4 effect decomposition ----
## Build interaction limma model ----
library(limma) 

med_dat <- dat %>%
  filter(Diagnosis %in% c("NCI", "AD")) %>%
  mutate(Diagnosis = droplevels(Diagnosis),
         diagnosis_bin = as.integer(Diagnosis == "AD")) # 0 = NCI, 1 = AD

expr_mat <- med_dat %>%
  dplyr::select(all_of(protein_vars)) %>%
  as.matrix() %>% t()
rownames(expr_mat) <- protein_vars
colnames(expr_mat) <- med_dat$projid_visit

design_mat <- model.matrix(~ 0 + Diagnosis + apoe4 + Diagnosis:apoe4 + age_at_visit + msex,
                           data = med_dat)

colnames(design_mat) <- gsub("Diagnosis|apoe4", "", colnames(design_mat)) %>%
  gsub(":", "_x_", .) %>%
  gsub("\\+", "pos", .)

contrasts_full <- makeContrasts(
  AD_vs_NCI_APOE4neg = AD - NCI, # AD effect in APOE4-
  AD_vs_NCI_APOE4pos = AD + AD_x_APOE4pos - NCI, # AD effect in APOE4+
  APOE4_effect_in_NCI = APOE4pos, # APOE4 effect in NCI
  APOE4_effect_in_AD = APOE4pos + AD_x_APOE4pos, # APOE4 effect in AD
  Interaction = AD_x_APOE4pos, # Differential APOE4 effect: AD vs NCI
  levels = design_mat
)

fit <- lmFit(expr_mat, design_mat) %>%
  contrasts.fit(contrasts_full) %>%
  eBayes()

## Extract results ----
p_val_threshold <- 0.05

get_res <- function(fit, coef) {
  topTable(fit, coef = coef, number = Inf) %>%
    rownames_to_column(var = "UniprotID") %>%
    mutate(negLog10FDR = -log10(adj.P.Val),
           Status = case_when(
             adj.P.Val < p_val_threshold & logFC >  0 ~ "Upregulated",
             adj.P.Val < p_val_threshold & logFC <  0 ~ "Downregulated",
             TRUE ~ "Not significant")) %>%
    left_join(protein_meta_clean, by = c("UniprotID" = "UniProt"))
}

res_AD_vs_NCI_APOE4neg <- get_res(fit, "AD_vs_NCI_APOE4neg")
res_AD_vs_NCI_APOE4pos <- get_res(fit, "AD_vs_NCI_APOE4pos")
res_APOE4_effect_in_NCI <- get_res(fit, "APOE4_effect_in_NCI")
res_APOE4_effect_in_AD <- get_res(fit, "APOE4_effect_in_AD")
res_interaction <- get_res(fit, "Interaction")

## Volcano plots ----
library(patchwork)

plot_volcano <- function(res, title) {
  sig <- res %>% filter(adj.P.Val < p_val_threshold)
  n_up <- sum(sig$logFC > 0)
  n_down <- sum(sig$logFC < 0)
  
  ggplot(res, aes(x = logFC, y = negLog10FDR)) +
    geom_point(aes(color = Status), alpha = 0.6, size = 1.2) +
    ggrepel::geom_text_repel(data = sig, aes(label = EntrezGeneSymbol),
                             size = 3, color = "black") +
    scale_color_manual(values = c("Upregulated" = "red", 
                                  "Downregulated" = "blue",
                                  "Not significant" = "grey")) +
    geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "grey") +
    labs(title = title,
         subtitle = sprintf("adj_p < 0.05  |  Up: %d   Down: %d", n_up, n_down),
         x = "log2 fold change", y = "-log10 adjusted P", color = "") +
    theme_classic() +
    theme(legend.position = "right")
}

p_AD_APOE4neg <- plot_volcano(res_AD_vs_NCI_APOE4neg, "AD vs NCI | APOE4-")
p_AD_APOE4pos <- plot_volcano(res_AD_vs_NCI_APOE4pos, "AD vs NCI | APOE4+")
p_APOE4_in_NCI <- plot_volcano(res_APOE4_effect_in_NCI, "APOE4+ vs APOE4- | NCI")
p_APOE4_in_AD <- plot_volcano(res_APOE4_effect_in_AD, "APOE4+ vs APOE4- | AD")
p_interaction <- plot_volcano(res_interaction, "Interaction: AD × APOE4")

p_volcano_grid <- (p_AD_APOE4neg  | p_AD_APOE4pos) / (p_APOE4_in_NCI | p_APOE4_in_AD) +
  plot_annotation(title = "AD and APOE4 Effect Decomposition")

# ggsave(file.path(pic_output_dir, "plasma_decomposition_volcano.png"), plot = p_volcano_grid, width = 15, height = 12)

## PCA on APOE4 effect DEPs, across strata ----
APOE4_DEPs <- res_APOE4_effect_in_AD %>%
  filter(adj.P.Val < p_val_threshold) %>%
  pull(UniprotID)

plot_pca_stratum <- function(dat, proteins, stratum, title) {
  dat_s <- dat %>% filter(Diagnosis == stratum)
  pca <- prcomp(dat_s[, proteins], center = TRUE, scale. = TRUE)
  var_exp <- summary(pca)$importance[2, 1:2] * 100
  
  as.data.frame(pca$x[, 1:2]) %>%
    mutate(apoe4 = dat_s$apoe4) %>%
    ggplot(aes(x = PC1, y = PC2, color = apoe4)) +
    geom_point(alpha = 0.7, size = 1.5) +
    labs(title = title,
         subtitle = sprintf("%d APOE4- vs APOE4+ (AD) DEPs", length(proteins)),
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2]),
         color = "") +
    theme_classic() +
    theme(legend.position = "bottom")
}

p_pca_AD  <- plot_pca_stratum(med_dat, APOE4_DEPs, "AD",  "PCA in AD cohort")
p_pca_NCI <- plot_pca_stratum(med_dat, APOE4_DEPs, "NCI", "PCA in NCI cohort")

p_pca_grid <- p_pca_NCI | p_pca_AD

# ggsave(file.path(pic_output_dir, "plasma_pca_APOE.png"), plot = p_pca_grid, width = 15, height = 8)

## PCA on AD effect DEPs, across strata ----
AD_DEPs <- res_AD_vs_NCI_APOE4neg %>%
  filter(adj.P.Val < p_val_threshold,
         abs(logFC) > 0.1) %>%
  pull(UniprotID)

plot_pca_apoe4_stratum <- function(dat, proteins, apoe4_stratum, title) {
  dat_s <- dat %>% filter(apoe4 == apoe4_stratum)
  pca <- prcomp(dat_s[, proteins], center = TRUE, scale. = TRUE)
  var_exp <- summary(pca)$importance[2, 1:2] * 100
  
  as.data.frame(pca$x[, 1:2]) %>%
    mutate(Diagnosis = dat_s$Diagnosis) %>%
    ggplot(aes(x = PC1, y = PC2, color = Diagnosis)) +
    geom_point(alpha = 0.7, size = 1.5) +
    labs(title = title,
         subtitle = sprintf("%d AD vs NCI (APOE4-) DEPs with |log FC| > 0.1", length(proteins)),
         x = sprintf("PC1 (%.1f%%)", var_exp[1]),
         y = sprintf("PC2 (%.1f%%)", var_exp[2]),
         color = "") +
    theme_classic() +
    theme(legend.position = "bottom")
}

p_pca_apoe4neg <- plot_pca_apoe4_stratum(med_dat, AD_DEPs, "APOE4-", "PCA in APOE4− cohort")
p_pca_apoe4pos <- plot_pca_apoe4_stratum(med_dat, AD_DEPs, "APOE4+", "PCA in APOE4+ cohort")

p_pca_AD_grid <- p_pca_apoe4neg | p_pca_apoe4pos

# ggsave(file.path(pic_output_dir, "plasma_pca_AD.png"), plot = p_pca_AD_grid, width = 15, height = 8)


# Mediation analysis 1: Different candidate proteins ----
library(mediation)
library(furrr)

## Pathway 1 (Upstream): APOE4 --> Proteins --> AD ----
# Filter candidate protein (with pre-AD APOE4 effect)
candidates_P1 <- res_APOE4_effect_in_NCI %>%
  filter(adj.P.Val < 0.05) %>%
  pull(UniprotID)

# Mediation model 
run_mediation_APOE <- function(protein) {
  df <- med_dat %>%
    dplyr::select(diagnosis_bin, apoe4, mediator = all_of(protein), age_at_visit, msex)
  
  med.fit <- lm(mediator ~ apoe4 + age_at_visit + msex, data = df)
  out.fit <- glm(diagnosis_bin ~ mediator + apoe4 + age_at_visit + msex, 
                 data = df, family = binomial(link = "logit"))
  
  med.out <- mediate(med.fit, out.fit, treat = "apoe4", mediator = "mediator",
                     control.value = "APOE4-", treat.value = "APOE4+",
                     boot = TRUE, sims = 1000)
  
  tibble(UniprotID = protein,
         ACME_est = med.out$d.avg,
         ACME_ci_low = med.out$d.avg.ci[1],
         ACME_ci_up = med.out$d.avg.ci[2],
         ACME_p = med.out$d.avg.p,
         ADE_est = med.out$z.avg,
         ADE_p = med.out$z.avg.p,
         prop_med = med.out$n.avg,
         prop_ci_low = med.out$n.avg.ci[1],
         prop_ci_up = med.out$n.avg.ci[2],
         total_effect = med.out$tau.coef,
         total_effect_p = med.out$tau.p)
}

plan(multisession, workers = parallel::detectCores() - 1)
res_mediation_APOE <- future_map_dfr(candidates_P1, run_mediation_APOE,
                                     .options = furrr_options(seed = 16)) %>%
  mutate(ACME_fdr = p.adjust(ACME_p, method = "BH")) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, TargetFullName,
                ACME_est, ACME_ci_low, ACME_ci_up, ACME_p, ACME_fdr, ADE_est, ADE_p,
                prop_med, prop_ci_low, prop_ci_up, total_effect, total_effect_p) %>%
  arrange(ACME_fdr)
plan(sequential)

sig_res_mediation_APOE <- res_mediation_APOE %>%
  filter(ACME_fdr < 0.05)
# write.csv(sig_res_mediation_APOE, file.path(table_output_dir, "plasma_mediation_APOE_prot_AD.csv"), row.names = FALSE)

# Forest plot
p_ACME_APOE <- ggplot(sig_res_mediation_APOE, aes(x = ACME_est, y = forcats::fct_reorder(EntrezGeneSymbol, ACME_est))) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = ACME_ci_low, xmax = ACME_ci_up), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "Average causal mediation effect (ACME)", y = NULL, 
       title = "Mediation Analysis: APOE4 → Protein → AD Diagnosis",
       subtitle = paste0(nrow(sig_res_mediation_APOE), " proteins with ACME FDR < 0.05 | 95% bootstrap CI")) +
  theme_classic()

# ggsave(file.path(pic_output_dir, "plasma_mediation_APOE_prot_AD.png"), p_ACME_APOE, width = 6, height = 6)


## Pathway 2 (Downstream): APOE4 --> AD --> Proteins ---- 
# Filter candidate protein (APOE4 mediated in AD)
candidates_P2 <- res_APOE4_effect_in_AD %>%
  filter(adj.P.Val < 0.05) %>%
  pull(UniprotID) # how much of APOE4 effect is mediated through AD diagnosis

# Mediation model 
run_mediation_AD <- function(protein) {
  df <- med_dat %>%
    dplyr::select(diagnosis_bin, apoe4, outcome_prot = all_of(protein), age_at_visit, msex)
  
  med.fit <- glm(diagnosis_bin ~ apoe4 + age_at_visit + msex, data = df, family = binomial(link = "logit"))
  out.fit <- lm(outcome_prot ~ diagnosis_bin + apoe4 + age_at_visit + msex, data = df)
  
  med.out <- mediate(med.fit, out.fit, treat = "apoe4", mediator = "diagnosis_bin",
                     control.value = "APOE4-", treat.value = "APOE4+",
                     boot = TRUE, sims = 1000)
  
  tibble(UniprotID = protein,
         ACME_est = med.out$d.avg,
         ACME_ci_low = med.out$d.avg.ci[1],
         ACME_ci_up = med.out$d.avg.ci[2],
         ACME_p = med.out$d.avg.p,
         ADE_est = med.out$z.avg,
         ADE_p = med.out$z.avg.p,
         prop_med = med.out$n.avg,
         prop_ci_low = med.out$n.avg.ci[1],
         prop_ci_up = med.out$n.avg.ci[2],
         total_effect = med.out$tau.coef,
         total_effect_p = med.out$tau.p)
}

plan(multisession, workers = parallel::detectCores() - 1)
res_mediation_AD <- future_map_dfr(candidates_P2, run_mediation_AD,
                                   .options = furrr_options(seed = 16)) %>%
  mutate(ACME_fdr = p.adjust(ACME_p, method = "BH")) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, TargetFullName,
                ACME_est, ACME_ci_low, ACME_ci_up, ACME_p, ACME_fdr, ADE_est, ADE_p,
                prop_med, prop_ci_low, prop_ci_up, total_effect, total_effect_p) %>%
  arrange(ACME_fdr)
plan(sequential)

sig_res_mediation_AD <- res_mediation_AD %>%
  filter(ACME_fdr < 0.05)
# write.csv(sig_res_mediation_AD, file.path(table_output_dir, "plasma_mediation_APOE_AD_prot.csv"), row.names = FALSE)

# Forest plot
p_ACME_AD <- ggplot(sig_res_mediation_AD, aes(x = ACME_est, y = forcats::fct_reorder(EntrezGeneSymbol, ACME_est))) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = ACME_ci_low, xmax = ACME_ci_up), height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "Average causal mediation effect (ACME)", y = NULL,
       title = "Mediation Analysis: APOE4 → AD Diagnosis → Protein",
       subtitle = paste0(nrow(sig_res_mediation_AD), " proteins with ACME FDR < 0.05 | 95% bootstrap CI")) +
  theme_classic()

# ggsave(file.path(pic_output_dir, "plasma_mediation_APOE_AD_prot.png"), p_ACME_AD, width = 6, height = 6)

## Pathways Overlap ----
mediation_classification <- bind_rows(
  res_mediation_APOE %>% mutate(pathway = "P1_upstream"),
  res_mediation_AD %>% mutate(pathway = "P2_downstream")
) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, pathway,
                ACME_est, ACME_ci_low, ACME_ci_up, ACME_fdr, prop_med) %>%
  pivot_wider(names_from = pathway,
              values_from = c(ACME_est, ACME_ci_low, ACME_ci_up, ACME_fdr, prop_med)) %>%
  mutate(tested_P1 = UniprotID %in% candidates_P1,
         tested_P2 = UniprotID %in% candidates_P2,
         sig_P1 = replace_na(ACME_fdr_P1_upstream < 0.05, FALSE),
         sig_P2 = replace_na(ACME_fdr_P2_downstream < 0.05, FALSE),
         classification = case_when(
           sig_P1 & sig_P2  ~ "Both pathways",
           sig_P1 & !sig_P2 &  tested_P2 ~ "Upstream only (tested, ns in P2)",
           sig_P1 & !sig_P2 & !tested_P2 ~ "Upstream only (P2 not tested)",
           !sig_P1 & sig_P2 &  tested_P1 ~ "Downstream only (tested, ns in P1)",
           !sig_P1 & sig_P2 & !tested_P1 ~ "Downstream only (P1 not tested)",
           TRUE ~ "Neither"
         ))

# Upstream only
dat_P1 <- mediation_classification %>%
  filter(classification == "Upstream only (P2 not tested)") %>%
  left_join(sig_res_mediation_APOE %>% dplyr::select(UniprotID),
            by = "UniprotID") %>%
  mutate(label = EntrezGeneSymbol,
         y = forcats::fct_reorder(label, ACME_est_P1_upstream))

p_A <- ggplot(dat_P1, aes(x = ACME_est_P1_upstream, y = y)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = ACME_ci_low_P1_upstream, xmax = ACME_ci_up_P1_upstream),
                 height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "ACME", y = NULL,
       title = "Upstream only: APOE4 → Protein → AD",
       subtitle = "proteins are not tested in other pathway") +
  theme_classic()

# Downstream only
dat_P2 <- mediation_classification %>%
  filter(classification == "Downstream only (tested, ns in P1)") %>%
  left_join(sig_res_mediation_AD %>% dplyr::select(UniprotID),
            by = "UniprotID") %>%
  mutate(label = EntrezGeneSymbol,
         y = forcats::fct_reorder(label, ACME_est_P2_downstream))

p_B <- ggplot(dat_P2, aes(x = ACME_est_P2_downstream, y = y)) +
  geom_point(size = 2) +
  geom_errorbarh(aes(xmin = ACME_ci_low_P2_downstream, xmax = ACME_ci_up_P2_downstream),
                 height = 0.2) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "ACME", y = NULL,
       title = "Downstream only: APOE4 → AD → Protein",
       subtitle = "proteins are tested but ns in other pathway") +
  theme_classic()

# Both
dat_both_wide <- mediation_classification %>%
  filter(classification == "Both pathways")

p_C <- ggplot(dat_both_wide, aes(x = ACME_est_P1_upstream, y = ACME_est_P2_downstream)) +
  geom_point(size = 2.5) +
  ggrepel::geom_text_repel(aes(label = EntrezGeneSymbol), size = 3) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  labs(x = "ACME - Upstream (APOE4 → Protein → AD)",
       y = "ACME - Downstream (APOE4 → AD → Protein)",
       title = "Proteins appear in both pathways",
       subtitle = "") +
  theme_classic()

# Save result
p_mediation_final <- p_A + p_B + p_C +
  plot_annotation(title = "Mediation Analysis")

# ggsave(file.path(pic_output_dir, "plasma_mediation_final.png"), p_mediation_final, width = 15, height = 6)


# Mediation analysis 2: All APOE4-associated proteins in both pathways ----
library(mediation)
library(furrr)

candidates_all <- res_limma %>%
  filter(adj.P.Val < 0.05) %>%
  pull(UniprotID)

## Pathway 1: APOE4 --> Protein --> AD ----
run_mediation_P1 <- function(protein) {
  df <- med_dat %>%
    dplyr::select(diagnosis_bin, apoe4, mediator = all_of(protein), age_at_visit, msex)
  
  med.fit <- lm(mediator ~ apoe4 + age_at_visit + msex, data = df)
  out.fit <- glm(diagnosis_bin ~ mediator + apoe4 + age_at_visit + msex,
                 data = df, family = binomial(link = "logit"))
  
  med.out <- mediate(med.fit, out.fit, treat = "apoe4", mediator = "mediator",
                     control.value = "APOE4-", treat.value = "APOE4+",
                     boot = TRUE, sims = 1000)
  
  tibble(UniprotID = protein,
         ACME_est = med.out$d.avg,
         ACME_ci_low = med.out$d.avg.ci[1],
         ACME_ci_up = med.out$d.avg.ci[2],
         ACME_p = med.out$d.avg.p,
         ADE_est = med.out$z.avg,
         ADE_p = med.out$z.avg.p,
         prop_med = med.out$n.avg,
         prop_ci_low = med.out$n.avg.ci[1],
         prop_ci_up = med.out$n.avg.ci[2],
         total_effect = med.out$tau.coef,
         total_effect_p = med.out$tau.p)
}

## Pathway 2: APOE4 --> AD --> Protein ----
run_mediation_P2 <- function(protein) {
  df <- med_dat %>%
    dplyr::select(diagnosis_bin, apoe4, outcome_prot = all_of(protein), age_at_visit, msex)
  
  med.fit <- glm(diagnosis_bin ~ apoe4 + age_at_visit + msex, data = df, family = binomial(link = "logit"))
  out.fit <- lm(outcome_prot ~ diagnosis_bin + apoe4 + age_at_visit + msex, data = df)
  
  med.out <- mediate(med.fit, out.fit, treat = "apoe4", mediator = "diagnosis_bin",
                     control.value = "APOE4-", treat.value = "APOE4+",
                     boot = TRUE, sims = 1000)
  
  tibble(UniprotID = protein,
         ACME_est = med.out$d.avg,
         ACME_ci_low = med.out$d.avg.ci[1],
         ACME_ci_up = med.out$d.avg.ci[2],
         ACME_p = med.out$d.avg.p,
         ADE_est = med.out$z.avg,
         ADE_p = med.out$z.avg.p,
         prop_med = med.out$n.avg,
         prop_ci_low = med.out$n.avg.ci[1],
         prop_ci_up = med.out$n.avg.ci[2],
         total_effect = med.out$tau.coef,
         total_effect_p = med.out$tau.p)
}

plan(multisession, workers = parallel::detectCores() - 1)

res_P1_all <- future_map_dfr(candidates_all, run_mediation_P1,
                             .options = furrr_options(seed = 16)) %>%
  mutate(ACME_fdr = p.adjust(ACME_p, method = "BH")) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt"))

res_P2_all <- future_map_dfr(candidates_all, run_mediation_P2,
                             .options = furrr_options(seed = 16)) %>%
  mutate(ACME_fdr = p.adjust(ACME_p, method = "BH")) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt"))

plan(sequential)

# write.csv(res_P1_all, file.path(table_output_dir, "plasma_mediation_all_P1.csv"), row.names = FALSE)
# write.csv(res_P2_all, file.path(table_output_dir, "plasma_mediation_all_P2.csv"), row.names = FALSE)

## Heatmap ----
star_label <- function(p) case_when(p < 0.001 ~ "***", p < 0.01 ~ "**", p < 0.05 ~ "*", TRUE ~ "")

path_labels <- c("Path1" = "APOE4 -> Protein -> AD", 
                 "Path2" = "APOE4 -> AD -> Protein")

sig_either <- union(res_P1_all %>% filter(ACME_fdr < 0.05) %>% pull(UniprotID),
                    res_P2_all %>% filter(ACME_fdr < 0.05) %>% pull(UniprotID))

heatmap_dat <- bind_rows(
  res_P1_all %>% filter(UniprotID %in% sig_either) %>% mutate(pathway = "Path1"),
  res_P2_all %>% filter(UniprotID %in% sig_either) %>% mutate(pathway = "Path2")
) %>%
  mutate(pct_med = prop_med * 100,
         sig_star = star_label(ACME_fdr))

protein_wide <- heatmap_dat %>%
  select(UniprotID, EntrezGeneSymbol, pathway, pct_med, ACME_fdr) %>%
  pivot_wider(names_from = pathway, values_from = c(pct_med, ACME_fdr)) %>%
  mutate(EntrezGeneSymbol = as.character(EntrezGeneSymbol),
         sig_P1 = replace_na(ACME_fdr_Path1 < 0.05, FALSE),
         sig_P2 = replace_na(ACME_fdr_Path2 < 0.05, FALSE),
         sort_group = case_when(sig_P1 & sig_P2 ~ 1, sig_P1 ~ 2, TRUE ~ 3),
         sort_val = replace_na(pct_med_Path1, 0),
         axis_color = case_when(sig_P1 & sig_P2 ~ "black", sig_P1 ~ "red", TRUE ~ "blue")) %>%
  arrange(sort_group, desc(sort_val))

protein_order <- protein_wide$EntrezGeneSymbol
axis_colors <- protein_wide$axis_color

heatmap_dat <- heatmap_dat %>%
  mutate(EntrezGeneSymbol = factor(EntrezGeneSymbol, levels = protein_order))

p_heatmap <- ggplot(heatmap_dat, aes(x = EntrezGeneSymbol, y = pathway, fill = pct_med)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sig_star), color = "black", size = 3, vjust = 0.75) +
  scale_fill_gradient2(
    low = "#9e9ac8", mid = "white", high = "#fb6a4a",
    midpoint = 0, name = "Percentage\nmediated(%)"
  ) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits = c("Path2", "Path1"), labels = path_labels) +
  labs(title = "APOE4-associated proteins that are involved in mediation paths",
       x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8, color = axis_colors),
    axis.text.y = element_text(size = 10),
    axis.line = element_blank(),
    axis.ticks  = element_blank(),
    legend.position = "right"
  )

# ggsave(file.path(pic_output_dir, "plasma_mediation_heatmap.png"), p_heatmap, width = 8, height = 6)
