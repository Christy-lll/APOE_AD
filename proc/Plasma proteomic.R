# Set up ----
suppressMessages(library(tidyverse))
library(here)

dat_output_dir <- here("cleaned_dat")
pic_output_dir <- here("output", "pics")
table_output_dir <- here("output", "tables")

proteomic_dat <- read.csv(here("raw", "OhNM2025_ROSMAP_plasma_Soma7k_protein_level_ANML_log10.csv"), 
                          header = T, strip.white = T, na.strings = "")
sample_meta <- read.csv(here("raw", "OhNM2025_ROSMAP_plasma_Soma7k_sample_metadata_1_.csv"), 
                        header = T, strip.white = T, na.strings = "") 
protein_meta <- read.csv(here("raw", "OhNM2025_ROSMAP_plasma_Soma7k_protein_metadata_1_.csv"), 
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
  select(UniProt, EntrezGeneSymbol, TargetFullName) %>%
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
  filter(Diagnosis %in% c("AD", "MCI", "NCI")) %>%# drop 30 samples 
  mutate(Diagnosis = factor(Diagnosis, 
                            levels = c("NCI", "MCI", "AD")), 
         apoe4 = factor(grepl("4", apoe_genotype),
                        levels = c(FALSE, TRUE),
                        labels = c("APOE4-", "APOE4+")),
         diag_apoe = factor(paste0(apoe4, " ", Diagnosis),
                            levels = c("APOE4- NCI", "APOE4+ NCI",
                                       "APOE4- MCI", "APOE4+ MCI",
                                       "APOE4- AD", "APOE4+ AD")),
         msex = factor(msex,
                       levels = c("0", "1"),
                       labels = c("Female", "Male")))

## Define variable groups ----
protein_vars <- colnames(proteomic_dat)[colnames(proteomic_dat) != "projid_visit"]
meta_vars <- setdiff(colnames(dat), protein_vars)

# final n = 793, proteins = 6402
# write.csv(dat, file.path(dat_output_dir, "plasma_dat_clean.csv"), row.names = FALSE)

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
  bind_cols(dat %>% select(projid_visit, diag_apoe)) 

p_pca_all <- ggplot(pca_all_df, aes(x = PC1, y = PC2, color = diag_apoe)) +
  geom_point() +
  labs(title = "PCA with 6402 plasma proteins",
       x = "PC1",
       y = "PC2",
       color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "plasma_allproteins_pca.png"), p_pca_all, width = 7, height = 6)

# Limma ----
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
  bind_cols(dat %>% select(projid_visit, diag_apoe)) 

p_pca_limma <- ggplot(pca_limma_df, aes(x = PC1, y = PC2, color = diag_apoe)) +
  geom_point() +
  labs(title = sprintf("PCA with %d differentially expressed plasma proteins (limma)", length(limma_sig_proteins)),
       x = "PC1",
       y = "PC2",
       color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "plasma_limma_pca.png"), p_pca_limma, width = 7, height = 6)

# Machine Learning ----
## Feature selection (on NCI) ----
mi_dat <- dat %>% 
  filter(Diagnosis == "NCI") %>%
  select(apoe4, all_of(protein_vars))

mi_scores_all <- FSelectorRcpp::information_gain(apoe4 ~ ., data = mi_dat, type = "infogain") %>%
  arrange(desc(importance))

mi_scores <- mi_scores_all %>%
  filter(importance > 0.08) %>% # threshold selected by sensitivity analysis --> 14 proteins
  left_join(protein_meta_clean, by = c("attributes" = "UniProt")) %>%
  select(attributes, EntrezGeneSymbol, TargetFullName, importance)

mi_sig_proteins <- mi_scores$attributes

# intersect(mi_sig_proteins, sig_res_limma$UniprotID)
# write.csv(mi_scores, file.path(table_output_dir, "plasma_MI.csv"), row.names = FALSE)

## PCA check ----
pca_MI <- prcomp(dat[, mi_sig_proteins], center = TRUE, scale. = TRUE)
pca_MI_df <- as.data.frame(pca_MI$x[, 1:2]) %>%
  bind_cols(dat %>% select(projid_visit, diag_apoe)) 

p_pca_MI <- ggplot(pca_MI_df, aes(x = PC1, y = PC2, color = diag_apoe)) +
  geom_point() +
  labs(title = sprintf("PCA with %d MI-selected plasma proteins", length(mi_sig_proteins)),
       x = "PC1",
       y = "PC2",
       color = "") +
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
  Sensitivity = mci_cm$byClass["Sensitivity"],
  Specificity = mci_cm$byClass["Specificity"],
  PPV = mci_cm$byClass["Pos Pred Value"],
  NPV = mci_cm$byClass["Neg Pred Value"],
  AUC = as.numeric(mci_roc$auc)
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
  Sensitivity = ad_cm$byClass["Sensitivity"],
  Specificity = ad_cm$byClass["Specificity"],
  PPV = ad_cm$byClass["Pos Pred Value"],
  NPV = ad_cm$byClass["Neg Pred Value"],
  AUC = as.numeric(ad_roc$auc)
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
  Sensitivity = male_cm$byClass["Sensitivity"],
  Specificity = male_cm$byClass["Specificity"],
  PPV = male_cm$byClass["Pos Pred Value"],
  NPV = male_cm$byClass["Neg Pred Value"],
  AUC = as.numeric(male_roc$auc)
)

## ML Metrics ----
all_metrics <- bind_rows(mci_metrics, ad_metrics, male_metrics)

# write.csv(all_metrics, file.path(table_output_dir, "plasma_ML_metrics.csv"), row.names = FALSE)

