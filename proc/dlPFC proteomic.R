# Set up ----
suppressMessages(library(tidyverse))
library(here)

dat_output_dir <- here("cleaned_dat")
pic_output_dir <- here("output", "pics")
table_output_dir <- here("output", "tables")

proteomic_dat <- read.csv(here("raw", "A1.PD-RAW_normalized_abundance-11672x500TMTchannels-no_batch_correction_1_.csv"), 
                          header = T, strip.white = T, na.strings = c("", "NA"), row.names = 1)
protein_meta <- read.csv(here("raw", "ROSMAP_assay_proteomics_TMTquantitation_metadata.csv"), 
                        header = T, strip.white = T, na.strings = c("", "NA"))
sample_meta <- read.csv(here("raw", "ROSMAP_clinical_1_.csv"), 
                        header = T, strip.white = T, na.strings = c("", "NA")) 

# Metadata ----
## Preprocessing  ----
meta_dat <- protein_meta %>%
  filter(!isAssayControl) %>%
  mutate(individualID = str_extract(specimenID, "[^.]+$")) %>%
  select(specimenID, individualID, batch, batchChannel) %>%
  left_join(sample_meta, by = "individualID") %>%
  mutate(
    apoe4 = factor(grepl("4", apoe_genotype), labels = c("APOE4-", "APOE4+")),
    cogdx = recode_values(cogdx, 1 ~ "NCI", 2 ~ "MCI", 4 ~ "AD"), # diagnosis at death
    cogdx = factor(cogdx, levels = c("NCI", "MCI", "AD")),
    msex = factor(msex, levels = c("0", "1"), labels = c("Female", "Male")),
    age_death_90plus = as.integer(age_death == "90+"), # age at death (numerical + indicator vars)
    age_death_num = as.numeric(if_else(age_death == "90+", "90", age_death))
    ) %>%
  filter(!is.na(cogdx)) %>%
  mutate(diag_apoe = factor(paste0(apoe4, " ", cogdx),
                            levels = c("APOE4- NCI", "APOE4+ NCI",
                                       "APOE4- MCI", "APOE4+ MCI",
                                       "APOE4- AD", "APOE4+ AD")))

## Demographic table ----
demographic <- meta_dat %>%
  select(diag_apoe, age_death_num, age_death_90plus, msex, educ, apoe_genotype) %>%
  gtsummary::tbl_summary(
    by = diag_apoe,
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age_death_num ~ "Age at Death (max 90)",
      age_death_90plus ~ "Age over 90",
      msex ~ "Sex",
      educ ~ "Education (years)",
      apoe_genotype ~ "APOE genotype"
    ),
    missing = "no"
  ) %>%
  add_p() %>%         
  add_overall(last = TRUE) %>%   
  bold_labels()

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "dlPFC demographic table.png"))

# Limma ----
## Preprocessing (log transform, impute) ----
expr_mat <- log2(as.matrix(proteomic_dat)[, meta_dat$batchChannel] + 1) # 11672 proteins x 374 individuals

median_impute <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

expr_mat <- expr_mat %>%
  .[rowMeans(is.na(.)) < 0.3, ] %>% # only keep proteins with <30% NA
  apply(1, median_impute) %>% # median imputation
  t() # final dim: 8007 proteins x 374 individuals

## Model fitting ----
library(limma)
# adjust for batch, sex, age at death and diagnosis
design_mat <- model.matrix(~ batch + msex + age_death_num + age_death_90plus + cogdx + apoe4, data = meta_dat) 
fit <- lmFit(expr_mat, design_mat)
fit <- eBayes(fit)

p_val_threshold <- 0.05
logfc_threshold <- 0.1

res_limma <- topTable(fit, coef = "apoe4APOE4+", number = Inf) %>%
  rownames_to_column(var = "Protein") %>%
  mutate(negLog10FDR = -log10(adj.P.Val),
         `Direction of change` = case_when(adj.P.Val < p_val_threshold & logFC > logfc_threshold ~ "Upregulated",
                                           adj.P.Val < p_val_threshold & logFC < -logfc_threshold ~ "Downregulated",
                                           TRUE ~ "Not significant"))

sig_res_limma <- res_limma %>%
  filter(adj.P.Val < p_val_threshold, abs(logFC) > logfc_threshold) 

# write.csv(sig_res_limma, file.path(table_output_dir, "dlPFC_limma.csv"), row.names = FALSE)

## Volcano plot ----
p_limma <- ggplot(res_limma, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = `Direction of change`), alpha = 0.6, size = 1.2) +
  ggrepel::geom_text_repel(data = sig_res_limma,
                  aes(label = Protein),
                  size = 3,
                  color = "black") +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
             linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(p_val_threshold), 
             linetype = "dashed", color = "grey") +
  labs(title = "Differentially Expressed dlPCF Proteins by APOE4 Status (limma)",
       subtitle = "Thresholds: adj_p < 0.05 & |log2 FC| > 0.1",
       x = "log2 fold change",
       y = "-log10 adjusted P",
       color = "Direction of change"
  ) +
  theme_classic() +
  theme(legend.position = "bottom")

# ggsave(file.path(pic_output_dir, "dlPFC_limma_volcano.png"), p_limma, width = 8, height = 6)

# Machine Learning ----
## Preprocessing (split, log transform, impute, batch correct) ----
# 70/30 stratified split
library(caret)
expr_mat <- log2(as.matrix(proteomic_dat)[, meta_dat$batchChannel] + 1) # 11672 proteins x 374 individuals

set.seed(16)
train_idx <- createDataPartition(meta_dat$apoe4, p = 0.7, list = FALSE)
meta_train <- meta_dat[train_idx, ]
meta_test <- meta_dat[-train_idx, ]
expr_train <- expr_mat[, train_idx]
expr_test <- expr_mat[, -train_idx]

# Median imputation: computed on train, applied to both
keep_proteins <- rownames(expr_train)[rowMeans(is.na(expr_train)) < 0.3]
expr_train <- expr_train[keep_proteins, ] # 7965 proteins x 263 individuals
expr_test <- expr_test[keep_proteins, ] # 7965 proteins x 111 individuals

protein_medians <- matrixStats::rowMedians(expr_train, na.rm = TRUE)

impute_with_reference <- function(mat, medians) {
  na_idx <- is.na(mat)
  mat[na_idx] <- medians[row(mat)[na_idx]]
  mat
}

expr_train <- impute_with_reference(expr_train, protein_medians)
expr_test <- impute_with_reference(expr_test,  protein_medians)

# Batch correction: fit on train, applied to both
design_train <- model.matrix(~ batch + apoe4 + cogdx + msex, data = meta_train) # double check the variables to be included
train_batch_fit <- lmFit(expr_train, design_train)

batch_levels <- levels(factor(meta_train$batch))
batch_coef_cols <- paste0("batch", batch_levels[-1])
batch_coefs <- coef(train_batch_fit)[, batch_coef_cols, drop = FALSE]

apply_batch_correction <- function(mat, meta, batch_coefs, batch_levels) {
  batch_factor <- factor(meta$batch, levels = batch_levels)
  batch_design <- model.matrix(~ batch_factor)[, -1, drop = FALSE]
  colnames(batch_design) <- colnames(batch_coefs)
  batch_effect <- batch_coefs %*% t(batch_design)
  mat - batch_effect
}

bc_train <- apply_batch_correction(expr_train, meta_train, batch_coefs, batch_levels) %>%
  t() %>% as.data.frame() %>% mutate(apoe4 = meta_train$apoe4)

bc_test <- apply_batch_correction(expr_test, meta_test, batch_coefs, batch_levels) %>%
  t() %>% as.data.frame() %>% mutate(apoe4 = meta_test$apoe4)

## Feature Selection (on train only) ----
mi_scores_all <- FSelectorRcpp::information_gain(apoe4 ~ ., data = bc_train, type = "infogain") %>%
  arrange(desc(importance))

mi_scores <- mi_scores_all %>% filter(importance > 0.035) # double check with the threshold
mi_sig_proteins <- mi_scores$attributes 
# write.csv(mi_scores, file.path(table_output_dir, "dlPFC_MI.csv"), row.names = FALSE)

## PCA check (train set) ----
pca_MI <- prcomp(bc_train[, mi_sig_proteins], center = TRUE, scale. = TRUE)
pca_MI_df <- as.data.frame(pca_MI$x[, 1:2]) %>%
  bind_cols(meta_train %>% select(apoe_genotype, apoe4, batch, cogdx))

p_pca_MI <- ggplot(pca_MI_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  # facet_wrap(~ cogdx, nrow = 1) +
  labs(title = sprintf("PCA with %d MI-selected proteins (Train, n = %d, MI > 0.035)", 
                       length(mi_sig_proteins), nrow(bc_train)),
       x = "PC1", y = "PC2", color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "dlPFC_MI_pca.png"), p_pca_MI, width = 7, height = 6)

## Model building (on 70% train) ----
train_dat <- bc_train %>%
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

## Test on 30% held-out set ----
test_dat <- bc_test %>%
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

test_class <- predict(model_fit, newdata = test_dat)
test_probs <- predict(model_fit, newdata = test_dat, type = "prob")[, "pos"]
test_cm <- confusionMatrix(test_class, test_dat$class, positive = "pos")
test_roc <- pROC::roc(response = test_dat$class, predictor = test_probs, levels = c("neg", "pos"))

test_metrics <- tibble(
  Model = "Random Forest (70/30 split, all diagnoses)",
  N_train = nrow(train_dat),
  N_test = nrow(test_dat),
  Sensitivity = test_cm$byClass["Sensitivity"],
  Specificity = test_cm$byClass["Specificity"],
  PPV = test_cm$byClass["Pos Pred Value"],
  NPV = test_cm$byClass["Neg Pred Value"],
  AUC = as.numeric(test_roc$auc)
)
