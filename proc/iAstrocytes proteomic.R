# Set up ----
suppressMessages(library(tidyverse))
library(here)

pic_output_dir <- here("output", "pics")
table_output_dir <- here("output", "tables")

proteomic_dat <- read.csv(here("raw", "ROSMAP-iA-TMT-MS.csv"), 
                          header = T, strip.white = T, na.strings = c("", "NA"), row.names = 1)
meta_dat <- read.csv(here("raw", "assay_meta_ROSMAP_iA_TMT-MS.csv"), 
                         header = T, strip.white = T, na.strings = c("", "NA"))

# Preprocessing  ----
## Metadata ----
meta_dat <- meta_dat %>%
  mutate(
    apoe4 = factor(grepl("4", apoe_genotype), labels = c("APOE4-", "APOE4+")),
    apoe_genotype = factor(apoe_genotype, levels = c("22", "23", "33", "24", "34", "44")),
    sex = factor(sex, levels = c("f", "m"), labels = c("Female", "Male")),
    dx = factor(dx, levels = c("LPNCI", "HPNCI", "AD")),
    pmAD = factor(pmAD, levels = c("0", "1"), labels = c("NCI", "AD")),
    AD_apoe = factor(paste0(apoe4, " ", pmAD), levels = c("APOE4- NCI", "APOE4+ NCI", "APOE4- AD", "APOE4+ AD"))
  ) 
    
## Demographic table ----
library(gtsummary)
demographic <- meta_dat %>%
  select(AD_apoe, sex, pmi, amyloid, tangles, dx, apoe_genotype) %>%
  tbl_summary(
    by = AD_apoe,
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      sex ~ "Sex",
      pmi ~ "PMI",
      amyloid ~ "Amyloid",
      tangles ~ "Tangles",
      dx ~ "Pathology",
      apoe_genotype ~ "APOE genotype"
    ),
    missing = "no"
  ) %>%
  add_p() %>%         
  add_overall(last = TRUE) %>%   
  bold_labels()

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "iAstrocytes demographic table.png"))

## Proteomic data ----
gene_symbols <- proteomic_dat[, 1, drop = FALSE] %>%
  rownames_to_column("UniprotID") %>%
  rename(GeneSymbol = symbol)

expr_mat <- as.matrix(proteomic_dat[, -1]) # 8346 proteins x 44 samples

# Log-transformation and normalisation check
anyNA(expr_mat) # no NA
summary(colMeans(expr_mat, na.rm = TRUE)) # comparable scale
boxplot(expr_mat, las = 2, cex.axis = 0.5, outline  = FALSE) # centered near 0

pca_all <- prcomp(t(expr_mat), center = TRUE, scale. = TRUE) # no significant batch effect
pca_all_df <- as.data.frame(pca_all$x[, 1:2]) %>%
  bind_cols(meta_dat %>% select(apoe_genotype, apoe4, batch))
p_pca_all <- ggplot(pca_all_df, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  labs(title = "PCA with 8346 plasma proteins",
       x = "PC1",
       y = "PC2",
       color = "") +
  theme_classic()

# Variance partition 
meta_vp <- meta_dat
colnames(expr_mat) <- meta_vp$SampleID
rownames(meta_vp) <- meta_vp$SampleID
# all(colnames(expr_mat) == rownames(meta_vp))
vp_form <- ~ (1 | batch) + (1 | apoe4) + (1 | pmAD) + (1 | sex) + pmi
vp_fit <- variancePartition::fitExtractVarPartModel(expr_mat, vp_form, meta_vp)

factor_labels <- c(apoe4 = "APOE4 Status",
                   batch = "Batch",
                   pmAD = "Diagnosis",
                   pmi = "PMI",
                   Residuals = "Residuals",
                   sex = "Sex")

vp_medians <- as.data.frame(vp_fit) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
  mutate(Factor = factor_labels[Factor]) %>%
  arrange(desc(Median_VarExplained))

p_varpart <- variancePartition::plotVarPart(vp_fit, main = "Variance Partition: dlPFC Proteomics")
p_varpart$data <- p_varpart$data %>%
  mutate(variable = factor_labels[as.character(variable)])

# write.csv(vp_medians, file.path(table_output_dir, "iAstrocytes_varpart_medians.csv"), row.names = FALSE)
# ggsave(file.path(pic_output_dir, "iAstrocytes_varpart_violin.png"), p_varpart, width = 8, height = 5)

# Limma ----
## Model fitting ---- 
library(limma)
design_mat <- model.matrix(~ apoe4 + sex + pmAD + batch, data = meta_dat)
fit <- lmFit(expr_mat, design_mat)
fit <- eBayes(fit)

p_val_threshold <- 0.05
logfc_threshold <- 0.5

res_limma <- topTable(fit, coef = "apoe4APOE4+", adjust = "fdr", number = Inf) %>%
  rownames_to_column(var = "UniprotID") %>%
  mutate(negLog10FDR = -log10(adj.P.Val),
         `Direction of change` = case_when(adj.P.Val < p_val_threshold & logFC > logfc_threshold ~ "Upregulated",
                                           adj.P.Val < p_val_threshold & logFC < -logfc_threshold ~ "Downregulated",
                                           TRUE ~ "Not significant")) %>%
  left_join(gene_symbols, by = "UniprotID") %>%
  select(UniprotID, GeneSymbol, everything())

sigfc_res_limma <- res_limma %>%
  filter(abs(logFC) > logfc_threshold) ## 5 proteins with |log FC| > 0.5

## Volcano plot ----
p_limma <- ggplot(res_limma, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = `Direction of change`), alpha = 0.6, size = 1) +
  ggrepel::geom_text_repel(data = sigfc_res_limma,
                           aes(label = GeneSymbol),
                           size = 3,
                           color = "black") +
  scale_color_manual(values = c("Upregulated" = "red",
                                "Downregulated" = "blue",
                                "Not significant" = "grey")) + 
  geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
             linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(p_val_threshold), 
             linetype = "dashed", color = "grey") +
  labs(title = "Differentially Expressed iAstrocytes Proteins by APOE4 Status (limma)",
       subtitle = "Thresholds: adj_p < 0.05 & |log2 FC| > 0.5",
       x = "log2 fold change",
       y = "-log10 adjusted P",
       color = "Direction of change") +
  theme_classic() +
  theme(legend.position = "bottom")

# ggsave(file.path(pic_output_dir, "iAstrocytes_limma_volcano.png"), p_limma, width = 10, height = 6)

# Machine learning ----
## Split train&test sets
library(caret)
set.seed(49)
train_idx <- createDataPartition(meta_dat$apoe4, p = 0.7, list = FALSE)
meta_train <- meta_dat[train_idx, ]
meta_test <- meta_dat[-train_idx, ]
dat_train <- expr_mat[, train_idx] %>% t() %>% as.data.frame() %>% mutate(apoe4 = meta_train$apoe4)
dat_test <- expr_mat[, -train_idx] %>% t() %>% as.data.frame() %>% mutate(apoe4 = meta_test$apoe4)

## Feature selection (train set) ----
mi_scores_all <- FSelectorRcpp::information_gain(apoe4 ~ ., data = dat_train, type = "infogain") %>%
  arrange(desc(importance))

mi_scores <- mi_scores_all %>% 
  filter(importance > 0.2) %>%
  left_join(gene_symbols, by = c("attributes" = "UniprotID")) %>%
  select(attributes, GeneSymbol, importance) 

mi_sig_proteins <- mi_scores$attributes 
# write.csv(mi_scores, file.path(table_output_dir, "iAstrocytes_MI.csv"), row.names = FALSE)

## PCA check (train set) ----
pca_MI <- prcomp(dat_train[, mi_sig_proteins], center = TRUE, scale. = TRUE)
pca_MI_df <- as.data.frame(pca_MI$x[, 1:2]) %>%
  bind_cols(meta_train %>% select(apoe_genotype, apoe4, batch, dx))

p_pca_MI <- ggplot(pca_MI_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA with %d MI-selected proteins (Train, n = %d, MI > 0.2)", 
                       length(mi_sig_proteins), nrow(dat_train)),
       x = "PC1", y = "PC2", color = "") +
  theme_classic()
# ggsave(file.path(pic_output_dir, "iAstrocyte_MI_pca.png"), p_pca_MI, width = 7, height = 6)

## Random Forest Model ----
# Model building (train set)
rf_train <- dat_train %>%
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

cv_ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(49)
model_fit <- train(
  class ~ .,
  data = rf_train,
  method = "rf",
  metric = "ROC",
  trControl = cv_ctrl,
  tuneGrid = expand.grid(.mtry = c(2, 3, 4, 5)),
  importance = TRUE,
  nodesize = 5, 
  classwt = c(pos = 21, neg = 11) # table(meta_train$apoe4)
)

# Test on 30% held-out set 
rf_test <- dat_test %>%
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

test_class <- predict(model_fit, newdata = rf_test)
test_probs <- predict(model_fit, newdata = rf_test, type = "prob")[, "pos"]
test_cm <- confusionMatrix(test_class, rf_test$class, positive = "pos")
test_roc <- pROC::roc(response = rf_test$class, predictor = test_probs, levels = c("neg", "pos"))

rf_metrics <- tibble(
  Model = "Random Forest",
  N_train = nrow(rf_train),
  N_test = nrow(rf_test),
  Sensitivity = test_cm$byClass["Sensitivity"],
  Specificity = test_cm$byClass["Specificity"],
  PPV = test_cm$byClass["Pos Pred Value"],
  NPV = test_cm$byClass["Neg Pred Value"],
  AUC = as.numeric(test_roc$auc)
)

# Check for overfitting
rf_cv_auc <- max(model_fit$results$ROC)
rf_test_auc <- as.numeric(test_roc$auc)
cat("Random forest - CV:", rf_cv_auc, "Test:", rf_test_auc, 
    "Gap:", rf_cv_auc - rf_test_auc, "\n")

## Lasso Model ----
# Model building (train set)
lasso_train <- dat_train %>%
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

cv_ctrl <- trainControl(
  method = "repeatedcv",
  number = 10,
  repeats = 10,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(49)
model_fit <- train(
  class ~ .,
  data = lasso_train,
  method = "glmnet",
  metric = "ROC",
  trControl = cv_ctrl,
  tuneGrid = expand.grid(alpha = 1, lambda = 10^seq(-4, 1, length = 50)),
  importance = TRUE,
  classwt = c(pos = 21, neg = 11) # table(meta_train$apoe4)
)

# Test on 30% held-out set
lasso_test <- dat_test %>%
  transmute(class = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("pos", "neg")),
            across(all_of(mi_sig_proteins)))

test_class <- predict(model_fit, newdata = lasso_test)
test_probs <- predict(model_fit, newdata = lasso_test, type = "prob")[, "pos"]
test_cm <- confusionMatrix(test_class, lasso_test$class, positive = "pos")
test_roc <- pROC::roc(response = lasso_test$class, predictor = test_probs, levels = c("neg", "pos"))

lasso_metrics <- tibble(
  Model = "LASSO",
  N_train = nrow(rf_train),
  N_test = nrow(rf_test),
  Sensitivity = test_cm$byClass["Sensitivity"],
  Specificity = test_cm$byClass["Specificity"],
  PPV = test_cm$byClass["Pos Pred Value"],
  NPV = test_cm$byClass["Neg Pred Value"],
  AUC = as.numeric(test_roc$auc)
)

# Check for overfitting
lasso_cv_auc <- max(model_fit$results$ROC)
lasso_test_auc <- as.numeric(test_roc$auc)
cat("LASSO - CV:", lasso_cv_auc, "Test:", lasso_test_auc, 
    "Gap:", lasso_cv_auc - lasso_test_auc, "\n")

## ML Metrics ----
all_metrics <- bind_rows(rf_metrics, lasso_metrics)

# write.csv(all_metrics, file.path(table_output_dir, "iAstrocytes_ML_metrics.csv"), row.names = FALSE)
