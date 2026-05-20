# Set up ----
suppressMessages(library(caret))
suppressMessages(library(furrr))
suppressMessages(library(gtsummary))
suppressMessages(library(here))
suppressMessages(library(limma))
suppressMessages(library(mediation))
suppressMessages(library(patchwork))
suppressMessages(library(tidyverse))

pic_output_dir <- here("output", "ROSMAP", "pics")
table_output_dir <- here("output", "ROSMAP", "tables")

proteomic_dat <- read.csv(here("raw", "ROSMAP", "OhNM2025_ROSMAP_plasma_Soma7k_protein_level_ANML_log10.csv"),
                          header = T, strip.white = T, na.strings = "")
sample_meta <- read.csv(here("raw", "ROSMAP", "OhNM2025_ROSMAP_plasma_Soma7k_sample_metadata_1_.csv"),
                        header = T, strip.white = T, na.strings = "")
protein_meta <- read.csv(here("raw", "ROSMAP", "OhNM2025_ROSMAP_plasma_Soma7k_protein_metadata_1_.csv"),
                         header = T, strip.white = T, na.strings = "")

# Preprocessing ----
## Map seqIDs to UniProtIDs ----
protein_col_mask <- !colnames(proteomic_dat) %in% "projid_visit"
colnames(proteomic_dat) <- sub("^X(\\d)", "\\1", colnames(proteomic_dat))
colnames(proteomic_dat) <- gsub("\\.", "-", colnames(proteomic_dat))
colnames(proteomic_dat)[protein_col_mask] <- protein_meta$UniProt[
  match(colnames(proteomic_dat)[protein_col_mask], protein_meta$SeqId)
]

# remove missing or duplicated uniprot ID
proteomic_dat <- proteomic_dat[, !is.na(colnames(proteomic_dat)) & 
                                 !duplicated(colnames(proteomic_dat)) & 
                                 colnames(proteomic_dat) != ""] 

protein_meta_clean <- protein_meta %>%
  dplyr::select(UniProt, EntrezGeneSymbol, TargetFullName) %>%
  distinct(UniProt, .keep_all = TRUE)

## Merge data & filter samples ----
dat <- merge(proteomic_dat, sample_meta, by = "projid_visit")

# check duplication in sample --> table(table(dat$projid)) 
dat <- dat %>%
  filter(!is.na(projid)) %>% 
  group_by(projid) %>%
  slice_min(order_by = Visit, n = 1, with_ties = FALSE) %>% # keep 1st visit if having more than 1
  ungroup()

# filter NAs
dat <- dat %>%
  filter(!is.na(apoe_genotype)) %>% # drop samples with no apoe_genotype
  filter(Diagnosis %in% c("AD", "MCI", "NCI")) %>% # drop samples with other diagnosis
  mutate(Diagnosis = factor(Diagnosis,
                            levels = c("NCI", "MCI", "AD")),
         apoe4 = factor(grepl("4", apoe_genotype),
                        levels = c(FALSE, TRUE),
                        labels = c("APOE4-", "APOE4+")),
         apoe_genotype = factor(apoe_genotype,
                                levels = c("22", "23", "33", "24", "34", "44")),
         msex = factor(msex,
                       levels = c("0", "1"),
                       labels = c("Female", "Male")))

# define variable groups
protein_vars <- colnames(proteomic_dat)[colnames(proteomic_dat) != "projid_visit"]
meta_vars <- setdiff(colnames(dat), protein_vars)
# final n = 793, proteins = 6402
rm(protein_meta, sample_meta, proteomic_dat, protein_col_mask)


# Demographic ----
demographic <- dat %>%
  dplyr::select(Diagnosis, age_at_visit, msex, educ, apoe_genotype) %>%
  tbl_summary(
    by = Diagnosis,
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
  add_overall(last = TRUE) %>%
  bold_labels()

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "plasma proteomics demographic table.png"))


# Limma ----
## Model fitting & results ----
expr_mat <- as.matrix(dat[, protein_vars]) 
rownames(expr_mat) <- dat$projid_visit
expr_mat <- t(expr_mat)
design_mat <- model.matrix(~ apoe4 + Diagnosis + age_at_visit + msex, data = dat) # adjust for diagnosis, age and sex 

fit <- lmFit(expr_mat, design_mat)
fit <- eBayes(fit)

p_val_threshold <- 0.05
logfc_threshold <- 0

res_limma <- topTable(fit, coef = "apoe4APOE4+", number = Inf) %>%
  rownames_to_column(var = "UniprotID") %>%
  mutate(negLog10FDR = -log10(adj.P.Val),
         `Direction of change` = case_when(
           adj.P.Val < p_val_threshold & logFC > logfc_threshold ~ "Upregulated",
           adj.P.Val < p_val_threshold & logFC < -logfc_threshold ~ "Downregulated",
           TRUE ~ "Not significant"))

sig_res_limma <- res_limma %>%
  filter(adj.P.Val < p_val_threshold, abs(logFC) > logfc_threshold) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, TargetFullName, everything())

# write.csv(sig_res_limma, file.path(table_output_dir, "plasma_limma.csv"), row.names = FALSE)

## Volcano plot ----
n_up <- sum(sig_res_limma$`Direction of change` == "Upregulated", na.rm = TRUE)
n_down <- sum(sig_res_limma$`Direction of change` == "Downregulated", na.rm = TRUE)

p_limma <- ggplot(res_limma, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = `Direction of change`), size = 1.2) +
  ggrepel::geom_text_repel(data = sig_res_limma,
                           aes(label = EntrezGeneSymbol),
                           size = 3,
                           color = "black") +
  scale_color_manual(values = c("Upregulated" = "indianred",
                                "Downregulated" = "steelblue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold),
             linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(p_val_threshold),
             linetype = "dashed", color = "grey") +
  labs(title = sprintf("%d Differentially Expressed Plasma Proteins by APOE4 Status",
                       nrow(sig_res_limma)),
       subtitle = sprintf("Threshold: adj_p < %.2f & |log2 FC| > %i  |  Up: %i  Down: %i", 
                          p_val_threshold, logfc_threshold, n_up, n_down),
       x = "log2 fold change",
       y = "-log10 adjusted P",
       color = "Direction of change") +
  theme_classic() +
  theme(legend.position = "bottom")

# ggsave(file.path(pic_output_dir, "plasma_limma_volcano.png"), p_limma, width = 10, height = 8)

## PCA ----
# With all proteins
pca_all <- prcomp(dat[, protein_vars], center = TRUE, scale. = TRUE)

pca_all_df <- as.data.frame(pca_all$x[, 1:2]) %>%
  bind_cols(dat %>% dplyr::select(projid_visit, apoe4, apoe_genotype))

p_pca_all <- ggplot(pca_all_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA with %d plasma proteins", length(protein_vars)),
       color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "plasma_all_proteins_pca.png"), p_pca_all, width = 7, height = 6)

# With DEPs only
limma_sig_proteins <- sig_res_limma$UniprotID
pca_limma <- prcomp(dat[, limma_sig_proteins], center = TRUE, scale. = TRUE)
pca_limma_df <- as.data.frame(pca_limma$x[, 1:2]) %>%
  bind_cols(dat %>% dplyr::select(projid_visit, apoe4, apoe_genotype))

p_pca_limma <- ggplot(pca_limma_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA with %d differentially expressed plasma proteins",
                       length(limma_sig_proteins)),
       color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "plasma_limma_pca.png"), p_pca_limma, width = 7, height = 6)


# Machine Learning ----
## Feature selection using mutual information (on NCI) ----
mi_dat <- dat %>%
  filter(Diagnosis == "NCI") %>% 
  dplyr::select(apoe4, all_of(protein_vars))

mi_scores_all <- FSelectorRcpp::information_gain(apoe4 ~ ., data = mi_dat, type = "infogain") %>%
  arrange(desc(importance)) %>%
  filter(importance > 0) %>%
  left_join(protein_meta_clean, by = c("attributes" = "UniProt")) %>%
  mutate(selected = importance > 0.1)

mi_scores <- mi_scores_all %>%
  filter(importance > 0.1) %>% # Use proteins with MI > 0.1 for machine learning
  dplyr::select(attributes, EntrezGeneSymbol, TargetFullName, importance)

mi_sig_proteins <- mi_scores$attributes
mi_sig_proteins_name <- mi_scores$EntrezGeneSymbol

# write.csv(mi_scores, file.path(table_output_dir, "plasma_MI.csv"), row.names = FALSE)

## Barplot & PCA ----
p_mi_bar <- ggplot(mi_scores_all, aes(x = importance,
                                      y = forcats::fct_reorder(EntrezGeneSymbol, importance),
                                      fill = selected)) +
  geom_col(alpha = 0.9) +
  geom_vline(xintercept = 0.1, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey"),
                    name = "Selected for ML") +
  labs(title = sprintf("%d Proteins with MI > 0", nrow(mi_scores_all)),
       x = "Mutual Information (MI)",
       y = NULL) +
  theme_classic()

# ggsave(file.path(pic_output_dir, "plasma_MI_barplot.png"), p_mi_bar, width = 10, height = 8)

pca_MI <- prcomp(dat[, mi_sig_proteins], center = TRUE, scale. = TRUE)
pca_MI_df <- as.data.frame(pca_MI$x[, 1:2]) %>%
  bind_cols(dat %>% dplyr::select(projid_visit, apoe4, apoe_genotype))

p_pca_MI <- ggplot(pca_MI_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA with %d MI-selected plasma proteins", length(mi_sig_proteins)),
       color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "plasma_MI_pca.png"), p_pca_MI, width = 7, height = 6)

## Build model on NCI ----
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

## Train on female, test on male ----
female_dat <- dat %>%
  filter(msex == "Female") %>%
  dplyr::select(apoe4, all_of(protein_vars)) %>%
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

## Train on ROSMAP, test on ADNI ----
ROSMAP_dat <- dat %>%
  dplyr::select(apoe4, all_of(protein_vars)) %>%
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins))) %>%
  rename(any_of(setNames(mi_sig_proteins, mi_sig_proteins_name))) %>%
  dplyr::select(-ZW10) # ADNI data does not have ZW10

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
# Load log10 transformed, normalised ADNI CSF proteomic data
ADNI_dat <- read.csv(here("raw", "ADNI_log_norm_dat.csv"), header = T) %>%
  transmute(class = factor(ifelse(apoe == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(any_of(mi_sig_proteins_name)))

ADNI_class <- predict(ROSMAP_model_fit, newdata = ADNI_dat)
ADNI_probs <- predict(ROSMAP_model_fit, newdata = ADNI_dat, type = "prob")[, "pos"]
ADNI_cm <- confusionMatrix(ADNI_class, ADNI_dat$class, positive = "pos")
ADNI_roc <- pROC::roc(response = ADNI_dat$class, predictor = ADNI_probs, levels = c("neg", "pos"))

ADNI_metrics <- tibble(
  Model = "Random Forest (APOE4- ADNI vs. APOE4+ ADNI)",
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

## ML metrics summary ----
all_metrics <- bind_rows(mci_metrics, ad_metrics, male_metrics, ADNI_metrics)

# write.csv(all_metrics, file.path(table_output_dir, "plasma_ML_metrics.csv"), row.names = FALSE)


# APOE4 / AD Effect Decomposition ----
NoMCI_dat <- dat %>%
  filter(Diagnosis %in% c("NCI", "AD")) %>%
  mutate(Diagnosis = droplevels(Diagnosis),
         diagnosis_bin = as.integer(Diagnosis == "AD")) # 0 = NCI, 1 = AD

noMCI_expr_mat <- NoMCI_dat %>%
  dplyr::select(all_of(protein_vars)) %>%
  as.matrix() %>% t()

rownames(noMCI_expr_mat) <- protein_vars
colnames(noMCI_expr_mat) <- NoMCI_dat$projid_visit

## Analysis 1: Limma on full cohorts with interaction ----
# fit model
design_mat_decomp <- model.matrix(~ 0 + Diagnosis + apoe4 + Diagnosis:apoe4 + age_at_visit + msex, # adjust for age and sex
                           data = NoMCI_dat)

colnames(design_mat_decomp) <- gsub("Diagnosis|apoe4", "", colnames(design_mat_decomp)) %>%
  gsub(":", "_x_", .) %>%
  gsub("\\+", "pos", .)

contrasts_full <- makeContrasts(
  AD_vs_NCI_APOE4neg = AD - NCI,                 # AD effect in APOE4-
  AD_vs_NCI_APOE4pos = AD + AD_x_APOE4pos - NCI, # AD effect in APOE4+
  APOE4_effect_in_NCI = APOE4pos,                # APOE4 effect in NCI
  APOE4_effect_in_AD = APOE4pos + AD_x_APOE4pos, # APOE4 effect in AD
  Interaction = AD_x_APOE4pos,                   # Differential APOE4 effect: AD vs NCI
  levels = design_mat_decomp
)

fit <- lmFit(noMCI_expr_mat, design_mat_decomp) %>%
  contrasts.fit(contrasts_full) %>%
  eBayes()

# extract results
p_val_threshold <- 0.05
logfc_threshold <- 0

get_res <- function(fit, coef) {
  topTable(fit, coef = coef, number = Inf) %>%
    rownames_to_column(var = "UniprotID") %>%
    mutate(negLog10FDR = -log10(adj.P.Val),
           Status = case_when(
             adj.P.Val < p_val_threshold & logFC > logfc_threshold ~ "Upregulated",
             adj.P.Val < p_val_threshold & logFC < -logfc_threshold ~ "Downregulated",
             TRUE ~ "Not significant")) %>%
    left_join(protein_meta_clean, by = c("UniprotID" = "UniProt"))
}

res_AD_vs_NCI_APOE4neg <- get_res(fit, "AD_vs_NCI_APOE4neg")
res_AD_vs_NCI_APOE4pos <- get_res(fit, "AD_vs_NCI_APOE4pos")
res_APOE4_effect_in_NCI <- get_res(fit, "APOE4_effect_in_NCI")
res_APOE4_effect_in_AD <- get_res(fit, "APOE4_effect_in_AD")
res_interaction <- get_res(fit, "Interaction")

# volcano plot
plot_volcano <- function(res, title) {
  sig <- res %>% filter(adj.P.Val < p_val_threshold)
  n_up <- sum(sig$logFC > 0)
  n_down <- sum(sig$logFC < 0)
  
  ggplot(res, aes(x = logFC, y = negLog10FDR)) +
    geom_point(aes(color = Status), size = 1.2) +
    ggrepel::geom_text_repel(data = sig, aes(label = EntrezGeneSymbol),
                             size = 3, color = "black") +
    scale_color_manual(values = c("Upregulated" = "indianred",
                                  "Downregulated" = "steelblue",
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

p_volcano_grid <- (p_APOE4_in_NCI | p_APOE4_in_AD) / (p_AD_APOE4neg | p_AD_APOE4pos) + 
  plot_annotation(title = "APOE4 and AD Effect Decomposition")

# ggsave(file.path(pic_output_dir, "plasma_decomposition_volcano.png"), plot = p_volcano_grid, width = 15, height = 12)

# PCA on APOE4 effect DEPs
APOE4_DEPs <- res_APOE4_effect_in_AD %>%
  filter(adj.P.Val < p_val_threshold) %>%
  pull(UniprotID)

plot_pca_stratum <- function(dat, proteins, stratum, title) {
  dat_s <- dat %>% filter(Diagnosis == stratum)
  pca <- prcomp(dat_s[, proteins], center = TRUE, scale. = TRUE)
  
  as.data.frame(pca$x[, 1:2]) %>%
    mutate(apoe4 = dat_s$apoe4) %>%
    ggplot(aes(x = PC1, y = PC2, color = apoe4)) +
    geom_point() +
    labs(title = title,
         subtitle = sprintf("%d APOE4- AD vs APOE4+ AD DEPs", length(proteins)),
         color = "") +
    theme_classic() +
    theme(legend.position = "bottom")
}

p_pca_AD <- plot_pca_stratum(NoMCI_dat, APOE4_DEPs, "AD", "PCA in AD cohort")
p_pca_NCI <- plot_pca_stratum(NoMCI_dat, APOE4_DEPs, "NCI", "PCA in NCI cohort")
p_pca_grid <- p_pca_NCI | p_pca_AD

# ggsave(file.path(pic_output_dir, "plasma_pca_APOE+-.png"), plot = p_pca_grid, width = 15, height = 8)

# PCA on AD effect DEPs
AD_DEPs <- res_AD_vs_NCI_APOE4neg %>%
  filter(adj.P.Val < p_val_threshold) %>%
  pull(UniprotID)

plot_pca_apoe4_stratum <- function(dat, proteins, apoe4_stratum, title) {
  dat_s <- dat %>% filter(apoe4 == apoe4_stratum)
  pca <- prcomp(dat_s[, proteins], center = TRUE, scale. = TRUE)
  
  as.data.frame(pca$x[, 1:2]) %>%
    mutate(Diagnosis = dat_s$Diagnosis) %>%
    ggplot(aes(x = PC1, y = PC2, color = Diagnosis)) +
    geom_point() +
    labs(title = title,
         subtitle = sprintf("%d APOE4- AD vs APOE4- NCI DEPs", length(proteins)),
         color = "") +
    theme_classic() +
    theme(legend.position = "bottom")
}

p_pca_apoe4neg <- plot_pca_apoe4_stratum(NoMCI_dat, AD_DEPs, "APOE4-", "PCA in APOE4- cohort")
p_pca_apoe4pos <- plot_pca_apoe4_stratum(NoMCI_dat, AD_DEPs, "APOE4+", "PCA in APOE4+ cohort")
p_pca_AD_grid <- p_pca_apoe4neg | p_pca_apoe4pos

# ggsave(file.path(pic_output_dir, "plasma_pca_ADNCI.png"), plot = p_pca_AD_grid, width = 15, height = 8)

# Analysis 2: Stratified Limma ----
# fit model
run_stratified_limma <- function(dat, protein_vars, stratum_var, stratum_val, contrast_coef) {
  dat_s <- dat %>%
    filter(.data[[stratum_var]] == stratum_val) %>%
    mutate(across(where(is.factor), droplevels))
  
  expr <- dat_s %>% select(all_of(protein_vars)) %>% as.matrix() %>% t()
  
  design_formula <- if (stratum_var == "apoe4") {
    ~ Diagnosis + age_at_visit + msex
  } else {
    ~ apoe4 + age_at_visit + msex
  }
  
  lmFit(expr, model.matrix(design_formula, data = dat_s)) %>%
    eBayes() %>%
    topTable(coef = contrast_coef, number = Inf) %>%
    rownames_to_column("UniprotID") %>%
    left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
    mutate(stratum = paste0(stratum_var, ": ", stratum_val),
           negLog10FDR = -log10(adj.P.Val),
           direction = case_when(adj.P.Val < p_val_threshold & logFC > 0 ~ "Upregulated",
                                 adj.P.Val < p_val_threshold & logFC < 0 ~ "Downregulated",
                                 TRUE ~ "Not significant"))
}

ad_by_apoe4 <- map(c("APOE4-", "APOE4+"),
                   ~ run_stratified_limma(NoMCI_dat, protein_vars, "apoe4", .x, "DiagnosisAD")) %>%
  setNames(c("APOE4neg", "APOE4pos"))

apoe4_by_dx <- map(c("NCI", "AD"),
                   ~ run_stratified_limma(NoMCI_dat, protein_vars, "Diagnosis", .x, "apoe4APOE4+")) %>%
  setNames(c("NCI", "AD"))

# extract results
list(ad_by_apoe4 = ad_by_apoe4, apoe4_by_dx = apoe4_by_dx) %>%
  map(~ map(.x, ~ count(.x, direction)))

dep_ad_noncarrier <- filter(ad_by_apoe4$APOE4neg, direction != "Not significant") %>% 
  pull(UniprotID)

dep_apoe4_nci <- filter(apoe4_by_dx$NCI, direction != "Not significant") %>% 
  pull(UniprotID)

# logFC correlation in DEPs
make_strat_scatter <- function(res_a, res_b, label_a, label_b, title) {
  joined <- inner_join(
    res_a %>% select(UniprotID, logFC_a = logFC),
    res_b %>% select(UniprotID, logFC_b = logFC, p_b = adj.P.Val),
    by = "UniprotID"
  ) %>% mutate(sig = p_b < 0.05)
  
  cor_test <- cor.test(joined$logFC_a, joined$logFC_b, method = "spearman")
  stats_label <- paste0("Spearman r = ", round(cor_test$estimate, 2),
                        ", p = ", formatC(cor_test$p.value, digits = 1))
  
  ggplot(joined, aes(x = logFC_a, y = logFC_b)) +
    geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "grey") +
    geom_point(aes(color = sig), size = 1.2) +
    geom_smooth(method = "lm", color = "black", linewidth = 0.6, se = TRUE) +
    scale_color_manual(values = c("TRUE" = "indianred", "FALSE" = "steelblue"),
                       labels = c("TRUE" = paste0("Sig in ", label_b), "FALSE" = "Not significant"),
                       guide = guide_legend(override.aes = list(size = 2))) +
    annotate("text", x = -Inf, y = Inf, hjust = -0.1, vjust = 1.5,
             label = stats_label, size = 3, fontface = "italic") +
    labs(title = title, x = paste0("logFC: ", label_a), y = paste0("logFC: ", label_b), color = NULL) +
    theme_classic() +
    theme(legend.position = "bottom")
}

p_ad_strat <- make_strat_scatter(
  ad_by_apoe4$APOE4neg %>% filter(UniprotID %in% dep_ad_noncarrier),
  ad_by_apoe4$APOE4pos %>% filter(UniprotID %in% dep_ad_noncarrier),
  "APOE4-", "APOE4+", "AD vs NCI logFC (among APOE4- DEPs)"
)

p_apoe4_strat <- make_strat_scatter(
  apoe4_by_dx$NCI %>% filter(UniprotID %in% dep_apoe4_nci),
  apoe4_by_dx$AD %>% filter(UniprotID %in% dep_apoe4_nci),
  "NCI", "AD", "APOE4+ vs APOE4- logFC (among NCI DEPs)"
)

# ggsave(file.path(pic_output_dir, "plasma_stratified_logFC_scatter.png"), p_ad_strat | p_apoe4_strat, width = 14, height = 6)


# Mediation Analysis ----
## APOE4-associated proteins ----
design_mat_med <- model.matrix(~ apoe4 + age_at_visit + msex, data = NoMCI_dat)

fit_med <- lmFit(noMCI_expr_mat, design_mat_med)
fit_med <- eBayes(fit_med)

candidate_prot <- topTable(fit_med, coef = "apoe4APOE4+", number = Inf) %>%
  rownames_to_column(var = "UniprotID") %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  filter(adj.P.Val < 0.05) %>%
  pull(UniprotID)

## Pathway 1 (upstream): APOE4 --> Protein --> AD ----
run_mediation_P1 <- function(protein) {
  df <- NoMCI_dat %>%
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

res_P1_all <- future_map_dfr(candidate_prot, run_mediation_P1,
                             .options = furrr_options(seed = 16)) %>%
  mutate(ACME_fdr = p.adjust(ACME_p, method = "BH")) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  arrange(ACME_fdr) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, TargetFullName, ACME_fdr, everything())

plan(sequential)

# write.csv(res_P1_all, file.path(table_output_dir, "plasma_mediation_prot_AD.csv"), row.names = FALSE)

## Pathway 2 (downstream): APOE4 --> AD --> Protein ----
run_mediation_P2 <- function(protein) {
  df <- NoMCI_dat %>%
    dplyr::select(diagnosis_bin, apoe4, outcome_prot = all_of(protein), age_at_visit, msex)
  
  med.fit <- glm(diagnosis_bin ~ apoe4 + age_at_visit + msex, data = df,
                 family = binomial(link = "logit"))
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

res_P2_all <- future_map_dfr(candidate_prot, run_mediation_P2,
                             .options = furrr_options(seed = 16)) %>%
  mutate(ACME_fdr = p.adjust(ACME_p, method = "BH")) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  arrange(ACME_fdr) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, TargetFullName, ACME_fdr, everything())

plan(sequential)

# write.csv(res_P2_all, file.path(table_output_dir, "plasma_mediation_AD_prot.csv"), row.names = FALSE)

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
  dplyr::select(UniprotID, EntrezGeneSymbol, pathway, pct_med, ACME_fdr) %>%
  pivot_wider(names_from = pathway, values_from = c(pct_med, ACME_fdr)) %>%
  mutate(EntrezGeneSymbol = as.character(EntrezGeneSymbol),
         sig_P1 = replace_na(ACME_fdr_Path1 < 0.05, FALSE),
         sig_P2 = replace_na(ACME_fdr_Path2 < 0.05, FALSE),
         sort_group = case_when(sig_P1 & sig_P2 ~ 1, sig_P1 ~ 2, sig_P2 ~ 3),
         sort_val = replace_na(pct_med_Path1, 0)) %>%
  arrange(sort_group, sort_val)

protein_order <- protein_wide$EntrezGeneSymbol

heatmap_dat <- heatmap_dat %>%
  mutate(EntrezGeneSymbol = factor(EntrezGeneSymbol, levels = protein_order))

p_heatmap <- ggplot(heatmap_dat, aes(x = EntrezGeneSymbol, y = pathway, fill = pct_med)) +
  geom_tile(color = "white", linewidth = 0.3) +
  geom_text(aes(label = sig_star), color = "black", size = 3, vjust = 0.75) +
  scale_fill_gradient2(
    low = "steelblue", mid = "white", high = "indianred",
    midpoint = 0, name = "Percentage\nmediated(%)"
  ) +
  scale_x_discrete(position = "bottom") +
  scale_y_discrete(limits = c("Path2", "Path1"), labels = path_labels) +
  labs(title = "APOE4-associated proteins involved in mediation pathways",
       x = NULL, y = NULL) +
  theme_classic() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 8),
    axis.text.y = element_text(size = 10),
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    legend.position = "right"
  )

# ggsave(file.path(pic_output_dir, "plasma_mediation_heatmap.png"), p_heatmap, width = 8, height = 6)
