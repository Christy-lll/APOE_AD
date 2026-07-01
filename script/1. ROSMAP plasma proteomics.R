# Set up ----
suppressMessages(library(caret))
suppressMessages(library(furrr))
suppressMessages(library(gtsummary))
suppressMessages(library(here))
suppressMessages(library(limma))
suppressMessages(library(mediation))
suppressMessages(library(patchwork))
suppressMessages(library(pROC))
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

proteomic_dat <- proteomic_dat[, !is.na(colnames(proteomic_dat)) &
                                 !duplicated(colnames(proteomic_dat)) &
                                 colnames(proteomic_dat) != ""]

protein_meta_clean <- protein_meta %>%
  dplyr::select(UniProt, EntrezGeneSymbol, TargetFullName) %>%
  distinct(UniProt, .keep_all = TRUE)

## Merge data & filter samples ----
dat <- merge(proteomic_dat, sample_meta, by = "projid_visit")

# keep first visit per individual
dat <- dat %>%
  filter(!is.na(projid)) %>%
  group_by(projid) %>%
  slice_min(order_by = Visit, n = 1, with_ties = FALSE) %>%
  ungroup()

dat <- dat %>%
  filter(!is.na(apoe_genotype),
         Diagnosis %in% c("AD", "MCI", "NCI")) %>%
  mutate(Diagnosis = factor(Diagnosis, levels = c("NCI", "MCI", "AD")),
         apoe4 = factor(grepl("4", apoe_genotype),
                        levels = c(FALSE, TRUE), labels = c("APOE ε4-", "APOE ε4+")),
         apoe_genotype = factor(apoe_genotype, levels = c("22", "23", "33", "24", "34", "44")),
         msex = factor(msex, levels = c("0", "1"), labels = c("Female", "Male")))

# define protein columns
protein_vars <- colnames(proteomic_dat)[colnames(proteomic_dat) != "projid_visit"]

rm(protein_meta, sample_meta, proteomic_dat, protein_col_mask)


# Demographic ----
demographic <- dat %>%
  dplyr::select(Diagnosis, age_at_visit, msex, apoe_genotype) %>%
  tbl_summary(
    by = Diagnosis,
    statistic = list(all_continuous() ~ "{mean} ± {sd}",
                     all_categorical() ~ "{n} ({p}%)"),
    digits = all_continuous() ~ 1,
    label = list(age_at_visit ~ "Age at Visit",
                 msex ~ "Sex",
                 apoe_genotype ~ "APOE genotype"),
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels()


# Limma ----
## Model fitting & results ----
p_val_threshold <- 0.05

expr_mat <- t(as.matrix(dat[, protein_vars]))
rownames(expr_mat) <- protein_vars
colnames(expr_mat) <- dat$projid_visit

design_mat <- model.matrix(~ apoe4 + Diagnosis + age_at_visit + msex, data = dat)

fit <- lmFit(expr_mat, design_mat) %>% eBayes()

res_limma <- topTable(fit, coef = "apoe4APOE ε4+", number = Inf) %>%
  rownames_to_column(var = "UniprotID") %>%
  mutate(negLog10FDR = -log10(adj.P.Val),
         `Direction of change` = case_when(
           adj.P.Val < p_val_threshold & logFC > 0 ~ "Upregulated",
           adj.P.Val < p_val_threshold & logFC < 0 ~ "Downregulated",
           TRUE ~ "Not significant"))

sig_res_limma <- res_limma %>%
  filter(adj.P.Val < p_val_threshold) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, TargetFullName, everything())

write.csv(sig_res_limma, file.path(table_output_dir, "plasma_limma.csv"), row.names = FALSE)

## Volcano plot ----
n_up <- sum(sig_res_limma$`Direction of change` == "Upregulated", na.rm = TRUE)
n_down <- sum(sig_res_limma$`Direction of change` == "Downregulated", na.rm = TRUE)

p_limma_volcano <- ggplot(res_limma, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = `Direction of change`), size = 1.2) +
  ggrepel::geom_text_repel(data = sig_res_limma,
                           aes(label = EntrezGeneSymbol),
                           size = 3, color = "black") +
  scale_color_manual(values = c("Upregulated" = "indianred",
                                "Downregulated" = "steelblue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "grey") +
  labs(title = sprintf("%d differentially abundant plasma proteins by APOE ε4 carriage",
                       nrow(sig_res_limma)),
       subtitle = sprintf("Up: %d  Down: %d", n_up, n_down),
       x = "log2 fold change",
       y = "-log10 FDR",
       color = NULL) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(pic_output_dir, "plasma_limma_all_volcano.tiff"), p_limma_volcano, width = 7, height = 6, dpi = 300, compression = "lzw")

## PCA ----
# all proteins
pca_all <- prcomp(dat[, protein_vars], center = TRUE, scale. = TRUE)
pca_all_df <- as.data.frame(pca_all$x[, 1:2]) %>%
  bind_cols(dat %>% dplyr::select(apoe4, apoe_genotype))

p_pca_all <- ggplot(pca_all_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA of all %d plasma proteins by APOE ε4 carriage", length(protein_vars)),
       color = "") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(pic_output_dir, "plasma_all_proteins_pca.tiff"), p_pca_all, width = 7, height = 6, dpi = 300, compression = "lzw")

# DAPs only
limma_sig_proteins <- sig_res_limma$UniprotID
pca_limma <- prcomp(dat[, limma_sig_proteins], center = TRUE, scale. = TRUE)
pca_limma_df <- as.data.frame(pca_limma$x[, 1:2]) %>%
  bind_cols(dat %>% dplyr::select(apoe4, apoe_genotype))

p_pca_limma <- ggplot(pca_limma_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA of %d APOE ε4-associated plasma DAPs by APOE ε4 carriage",
                       length(limma_sig_proteins)),
       color = "") +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(pic_output_dir, "plasma_limma_pca.tiff"), p_pca_limma, width = 7, height = 6, dpi = 300, compression = "lzw")


# Machine Learning (random forest classifier) ----
## MI-based feature selection (on NCI) ----
mi_dat <- dat %>%
  filter(Diagnosis == "NCI") %>%
  dplyr::select(apoe4, all_of(protein_vars))

mi_scores_all <- FSelectorRcpp::information_gain(apoe4 ~ ., data = mi_dat, type = "infogain") %>%
  arrange(desc(importance)) %>%
  filter(importance > 0) %>%
  left_join(protein_meta_clean, by = c("attributes" = "UniProt")) %>%
  mutate(selected = importance > 0.1)

mi_scores <- mi_scores_all %>%
  filter(importance > 0.1) %>%
  dplyr::select(attributes, EntrezGeneSymbol, TargetFullName, importance)

mi_sig_proteins <- mi_scores$attributes
mi_sig_proteins_name <- mi_scores$EntrezGeneSymbol

write.csv(mi_scores, file.path(table_output_dir, "plasma_MI.csv"), row.names = FALSE)

## Train on NCI ----
train_dat <- mi_dat %>%
  transmute(class = factor(ifelse(apoe4 == "APOE ε4+", "pos", "neg"), levels = c("pos", "neg")),
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
model_fit <- train(class ~ ., data = train_dat, method = "rf",
                   metric = "ROC", trControl = cv_ctrl, tuneLength = 10, importance = TRUE)

## Test on MCI ----
test_MCI <- dat %>%
  filter(Diagnosis == "MCI") %>%
  transmute(class = factor(ifelse(apoe4 == "APOE ε4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

mci_probs <- predict(model_fit, newdata = test_MCI, type = "prob")[, "pos"]
mci_class <- predict(model_fit, newdata = test_MCI)
mci_cm <- confusionMatrix(mci_class, test_MCI$class, positive = "pos")
mci_roc <- roc(response = test_MCI$class, predictor = mci_probs, levels = c("neg", "pos"))

mci_metrics <- tibble(
  Model = "APOE ε4+ MCI vs APOE ε4- MCI",
  N_train = nrow(train_dat),
  N_test = nrow(test_MCI),
  N_features = length(mi_sig_proteins),
  Sensitivity = mci_cm$byClass["Sensitivity"],
  Specificity = mci_cm$byClass["Specificity"],
  PPV = mci_cm$byClass["Pos Pred Value"],
  NPV = mci_cm$byClass["Neg Pred Value"],
  Train_AUC = model_fit$results %>% filter(mtry == model_fit$bestTune$mtry) %>% pull(ROC),
  Test_AUC = as.numeric(mci_roc$auc),
  AUC_gap = Train_AUC - Test_AUC
)

## Test on AD ----
test_AD <- dat %>%
  filter(Diagnosis == "AD") %>%
  transmute(class = factor(ifelse(apoe4 == "APOE ε4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

ad_probs <- predict(model_fit, newdata = test_AD, type = "prob")[, "pos"]
ad_class <- predict(model_fit, newdata = test_AD)
ad_cm <- confusionMatrix(ad_class, test_AD$class, positive = "pos")
ad_roc <- roc(response = test_AD$class, predictor = ad_probs, levels = c("neg", "pos"))

ad_metrics <- tibble(
  Model = "APOE ε4+ AD vs APOE ε4- AD",
  N_train = nrow(train_dat),
  N_test = nrow(test_AD),
  N_features = length(mi_sig_proteins),
  Sensitivity = ad_cm$byClass["Sensitivity"],
  Specificity = ad_cm$byClass["Specificity"],
  PPV = ad_cm$byClass["Pos Pred Value"],
  NPV = ad_cm$byClass["Neg Pred Value"],
  Train_AUC = model_fit$results %>% filter(mtry == model_fit$bestTune$mtry) %>% pull(ROC),
  Test_AUC = as.numeric(ad_roc$auc),
  AUC_gap = Train_AUC - Test_AUC
)

## Train on female, test on male ----
female_dat <- dat %>%
  filter(msex == "Female") %>%
  transmute(class = factor(ifelse(apoe4 == "APOE ε4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

set.seed(16)
sex_model_fit <- train(class ~ ., data = female_dat, method = "rf",
                       metric = "ROC", trControl = cv_ctrl, tuneLength = 10, importance = TRUE)

male_dat <- dat %>%
  filter(msex == "Male") %>%
  transmute(class = factor(ifelse(apoe4 == "APOE ε4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

male_probs <- predict(sex_model_fit, newdata = male_dat, type = "prob")[, "pos"]
male_class <- predict(sex_model_fit, newdata = male_dat)
male_cm <- confusionMatrix(male_class, male_dat$class, positive = "pos")
male_roc <- roc(response = male_dat$class, predictor = male_probs, levels = c("neg", "pos"))

male_metrics <- tibble(
  Model = "APOE ε4+ Male vs APOE ε4- Male",
  N_train = nrow(female_dat),
  N_test = nrow(male_dat),
  N_features = length(mi_sig_proteins),
  Sensitivity = male_cm$byClass["Sensitivity"],
  Specificity = male_cm$byClass["Specificity"],
  PPV = male_cm$byClass["Pos Pred Value"],
  NPV = male_cm$byClass["Neg Pred Value"],
  Train_AUC = sex_model_fit$results %>% filter(mtry == sex_model_fit$bestTune$mtry) %>% pull(ROC),
  Test_AUC = as.numeric(male_roc$auc),
  AUC_gap = Train_AUC - Test_AUC
)

## Train on ROSMAP, test on ADNI (CSF) ----
ROSMAP_dat <- dat %>%
  transmute(class = factor(ifelse(apoe4 == "APOE ε4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins))) %>%
  rename(any_of(setNames(mi_sig_proteins, mi_sig_proteins_name))) %>%
  dplyr::select(-ZW10) # ZW10 not available in ADNI

set.seed(16)
ROSMAP_model_fit <- train(class ~ ., data = ROSMAP_dat, method = "rf",
                          metric = "ROC", trControl = cv_ctrl, tuneLength = 10, importance = TRUE)

# log10-transformed, normalised ADNI CSF proteomic data
ADNI_dat <- read.csv(here("raw", "ADNI_log_norm_dat.csv"), header = T) %>%
  transmute(class = factor(ifelse(apoe == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(any_of(mi_sig_proteins_name)))

ADNI_probs <- predict(ROSMAP_model_fit, newdata = ADNI_dat, type = "prob")[, "pos"]
ADNI_class <- predict(ROSMAP_model_fit, newdata = ADNI_dat)
ADNI_cm <- confusionMatrix(ADNI_class, ADNI_dat$class, positive = "pos")
ADNI_roc <- roc(response = ADNI_dat$class, predictor = ADNI_probs, levels = c("neg", "pos"))

ADNI_metrics <- tibble(
  Model = "APOE ε4+ CSF vs APOE ε4- CSF",
  N_train = nrow(ROSMAP_dat),
  N_test = nrow(ADNI_dat),
  N_features = ncol(ADNI_dat) - 1,
  Sensitivity = ADNI_cm$byClass["Sensitivity"],
  Specificity = ADNI_cm$byClass["Specificity"],
  PPV = ADNI_cm$byClass["Pos Pred Value"],
  NPV = ADNI_cm$byClass["Neg Pred Value"],
  Train_AUC = ROSMAP_model_fit$results %>% filter(mtry == ROSMAP_model_fit$bestTune$mtry) %>% pull(ROC),
  Test_AUC = as.numeric(ADNI_roc$auc),
  AUC_gap = Train_AUC - Test_AUC
)

## ML metrics summary ----
all_metrics <- bind_rows(mci_metrics, ad_metrics, male_metrics, ADNI_metrics)

auc_plot_dat <- all_metrics %>%
  dplyr::select(Model, Train_AUC, Test_AUC) %>%
  pivot_longer(c(Train_AUC, Test_AUC), names_to = "Set", values_to = "AUC") %>%
  mutate(Set = factor(recode(Set, Train_AUC = "Train", Test_AUC = "Test"),
                      levels = c("Train", "Test")),
         Model = factor(Model, levels = rev(unique(all_metrics$Model))))

p_auc <- ggplot(auc_plot_dat, aes(x = AUC, y = Model, color = Set, shape = Set)) +
  geom_point(size = 3.5) +
  scale_x_continuous(limits = c(0.8, 1.01), breaks = seq(0.8, 1.0, by = 0.05)) +
  scale_color_manual(values = c("Train" = "steelblue", "Test" = "indianred")) +
  labs(title = "Random forest classifier performance",
       x = "AUC",
       y = NULL,
       color = NULL,
       shape = NULL) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(pic_output_dir, "plasma_ML_auc_dotplot.tiff"), p_auc, width = 7, height = 6, dpi = 300, compression = "lzw")


# Mediation Analysis ----
## Filter to NCI and AD only ----
NoMCI_dat <- dat %>%
  filter(Diagnosis %in% c("NCI", "AD")) %>%
  mutate(Diagnosis = droplevels(Diagnosis),
         diagnosis_bin = as.integer(Diagnosis == "AD"))

noMCI_expr_mat <- NoMCI_dat %>%
  dplyr::select(all_of(protein_vars)) %>%
  as.matrix() %>% t()

rownames(noMCI_expr_mat) <- protein_vars
colnames(noMCI_expr_mat) <- NoMCI_dat$projid_visit

## Identify APOE4-associated proteins in NCI + AD ----
p_val_threshold <- 0.05
design_mat_med <- model.matrix(~ apoe4 + age_at_visit + msex, data = NoMCI_dat)
fit_med <- lmFit(noMCI_expr_mat, design_mat_med) %>% eBayes()

all_prot_med <- topTable(fit_med, coef = "apoe4APOE ε4+", number = Inf) %>%
  rownames_to_column(var = "UniprotID") %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  mutate(negLog10FDR = -log10(adj.P.Val),
         `Direction of change` = case_when(
           adj.P.Val < p_val_threshold & logFC > 0 ~ "Upregulated",
           adj.P.Val < p_val_threshold & logFC < 0 ~ "Downregulated",
           TRUE ~ "Not significant"))

candidate_prot <- all_prot_med %>% filter(adj.P.Val < p_val_threshold)

# Volcano plot
n_up <- sum(candidate_prot$`Direction of change` == "Upregulated", na.rm = TRUE)
n_down <- sum(candidate_prot$`Direction of change` == "Downregulated", na.rm = TRUE)

p_limma_med <- ggplot(all_prot_med, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = `Direction of change`), size = 1.5) +
  ggrepel::geom_text_repel(data = candidate_prot,
                           aes(label = EntrezGeneSymbol),
                           size = 3, color = "black") +
  scale_color_manual(values = c("Upregulated" = "indianred",
                                "Downregulated" = "steelblue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "grey") +
  labs(title = sprintf("%d differentially abundant plasma proteins by APOE ε4 carriage",
                       nrow(candidate_prot)),
       subtitle = sprintf("Up: %d  Down: %d", n_up, n_down),
       x = "log2 fold change",
       y = "-log10 FDR",
       color = NULL) +
  theme_classic() +
  theme(legend.position = "bottom")

ggsave(file.path(pic_output_dir, "plasma_limma_mediation_volcano.tiff"), p_limma_med, width = 7, height = 6, dpi = 300, compression = "lzw")

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

res_P1_all <- future_map_dfr(candidate_prot$UniprotID, run_mediation_P1,
                             .options = furrr_options(seed = 16)) %>%
  mutate(ACME_fdr = p.adjust(ACME_p, method = "BH")) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  arrange(ACME_fdr) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, TargetFullName, ACME_fdr, everything())

plan(sequential)

write.csv(res_P1_all, file.path(table_output_dir, "plasma_mediation_prot_AD.csv"), row.names = FALSE)

## Pathway 2 (downstream): APOE4 --> AD --> Protein ----
run_mediation_P2 <- function(protein) {
  df <- NoMCI_dat %>%
    dplyr::select(diagnosis_bin, apoe4, outcome_prot = all_of(protein), age_at_visit, msex)
  
  med.fit <- glm(diagnosis_bin ~ apoe4 + age_at_visit + msex,
                 data = df, family = binomial(link = "logit"))
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

res_P2_all <- future_map_dfr(candidate_prot$UniprotID, run_mediation_P2,
                             .options = furrr_options(seed = 16)) %>%
  mutate(ACME_fdr = p.adjust(ACME_p, method = "BH")) %>%
  left_join(protein_meta_clean, by = c("UniprotID" = "UniProt")) %>%
  arrange(ACME_fdr) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, TargetFullName, ACME_fdr, everything())

plan(sequential)

write.csv(res_P2_all, file.path(table_output_dir, "plasma_mediation_AD_prot.csv"), row.names = FALSE)

## Forest plot ----
mediation_plot_dat <- bind_rows(
  res_P1_all %>% mutate(pathway = "P1: APOE ε4 -> Protein -> AD"),
  res_P2_all %>% mutate(pathway = "P2: APOE ε4 -> AD -> Protein")
) %>%
  filter(ACME_fdr < 0.05) %>%
  dplyr::select(UniprotID, EntrezGeneSymbol, pathway, ACME_est, ACME_ci_low, ACME_ci_up, ACME_fdr) %>%
  mutate(direction = ifelse(ACME_est > 0, "Positive", "Negative"))

sig_proteins <- mediation_plot_dat %>% distinct(UniprotID, EntrezGeneSymbol)

plot_grid <- expand.grid(
  UniprotID = sig_proteins$UniprotID,
  pathway = c("P1: APOE ε4 -> Protein -> AD", "P2: APOE ε4 -> AD -> Protein"),
  stringsAsFactors = FALSE
) %>%
  left_join(sig_proteins, by = "UniprotID") %>%
  left_join(mediation_plot_dat, by = c("UniprotID", "EntrezGeneSymbol", "pathway"))

classification <- mediation_plot_dat %>%
  group_by(UniprotID, EntrezGeneSymbol) %>%
  summarise(n_pathways = n(),
            signs = paste(sort(unique(sign(ACME_est))), collapse = ","),
            .groups = "drop") %>%
  mutate(group = case_when(
    n_pathways == 2 & signs == "-1,1" ~ "Both (competitive)",
    n_pathways == 2 & signs %in% c("1", "-1") ~ "Both (complementary)",
    n_pathways == 1 & UniprotID %in%
      (res_P1_all %>% filter(ACME_fdr < 0.05) %>% pull(UniprotID)) ~ "Upstream only",
    n_pathways == 1 ~ "Downstream only"
  ))

plot_grid <- plot_grid %>%
  left_join(classification %>% dplyr::select(UniprotID, group), by = "UniprotID") %>%
  mutate(group = factor(group, levels = c("Both (complementary)", "Both (competitive)",
                                          "Upstream only", "Downstream only")),
         protein_label = forcats::fct_reorder2(EntrezGeneSymbol, group, ACME_est, .na_rm = FALSE))

p_mediation <- ggplot(plot_grid, aes(x = ACME_est, y = protein_label)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_errorbarh(aes(xmin = ACME_ci_low, xmax = ACME_ci_up, color = direction),
                 height = 0.25, na.rm = TRUE) +
  geom_point(aes(color = direction), size = 2.5, na.rm = TRUE) +
  scale_color_manual(values = c("Positive" = "indianred", "Negative" = "steelblue"),
                     na.translate = FALSE, name = "ACME direction") +
  facet_grid(group ~ pathway, scales = "free_y", space = "free_y", switch = "y") +
  labs(title = "Mediation analysis",
       x = "Average Causal Mediation Effect (ACME)",
       y = NULL) +
  theme_classic() +
  theme(strip.placement = "outside",
        strip.background = element_rect(fill = "grey90", color = NA),
        strip.text.y.left = element_text(angle = 0, hjust = 1),
        panel.spacing = unit(0.3, "lines"),
        legend.position = "bottom")

ggsave(file.path(pic_output_dir, "plasma_mediation_forestplot.tiff"), p_mediation, width = 7, height = 6, dpi = 300, compression = "lzw")


# Enrichment Analysis ----
## Export gene lists for NetworkAnalyst ----
background_universe <- all_prot_med %>% distinct(UniprotID) %>% pull(UniprotID)
dep_gene_list <- candidate_prot %>% distinct(UniprotID) %>% pull(UniprotID)

writeLines(background_universe, file.path(table_output_dir, "plasma_networkanalyst_background.txt"))
writeLines(dep_gene_list, file.path(table_output_dir, "plasma_networkanalyst_dep.txt"))

## GO:BP results ----
GO_res <- read.csv(file.path(table_output_dir, "plasma_GOBP.csv")) %>%
  filter(FDR < 0.05) %>%
  distinct(Pathway, .keep_all = TRUE) %>%
  mutate(negLog10FDR = -log10(FDR),
         Pathway = forcats::fct_reorder(Pathway, negLog10FDR))

p_GO <- ggplot(GO_res, aes(x = negLog10FDR, y = Pathway)) +
  geom_col(fill = "indianred") +
  labs(title = "GO: biological processes",
       x = "-log10 FDR",
       y = NULL) +
  theme_classic()

ggsave(file.path(pic_output_dir, "plasma_GO_barplot.tiff"), p_GO, width = 7, height = 6, dpi = 300, compression = "lzw")

## Reactome results ----
reactome_res <- read.csv(file.path(table_output_dir, "plasma_Reactome.csv")) %>%
  filter(FDR < 0.05) %>%
  mutate(negLog10FDR = -log10(FDR),
         Pathway = forcats::fct_reorder(Pathway, negLog10FDR))

p_reactome <- ggplot(reactome_res, aes(x = negLog10FDR, y = Pathway)) +
  geom_col(fill = "steelblue") +
  labs(title = "Reactome",
       x = "-log10 FDR",
       y = NULL) +
  theme_classic()

ggsave(file.path(pic_output_dir, "plasma_Reactome_barplot.tiff"), p_reactome, width = 7, height = 6, dpi = 300, compression = "lzw")
