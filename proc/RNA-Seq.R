# Set up ----
suppressMessages(library(tidyverse))
library(here)

pic_output_dir <- here("output", "pics")
table_output_dir <- here("output", "tables")

expr_dat <- read_tsv(here("raw", "ROSMAP_Normalized_counts__CQN_.tsv")) %>%
  column_to_rownames(var = "feature")
count_dat <- read_tsv(here("raw", "ROSMAP_Filtered_counts__greater_than_1cpm_.tsv")) %>%
  column_to_rownames(var = "feature")

ID_meta <- read.csv(here("raw", "ROSMAP_biospecimen_metadata.csv"), 
                    header = T, strip.white = T, na.strings = "NA") # map specimen ID to individual ID
rna_meta <- read.csv(here("raw", "ROSMAP_assay_rnaSeq_metadata.csv"), 
                     header = T, strip.white = T, na.strings = "NA") # specimen ID, batch, RIN, etc
sample_meta <- read.csv(here("raw", "ROSMAP_clinical_1_.csv"), 
                        header = T, strip.white = T, na.strings = "NA") # individual ID, demographic 

# Preprocessing ----
## Metadata ----
meta_dat <- rna_meta %>% 
  left_join(ID_meta[, c("specimenID", "individualID", "tissue", "assay")], by = "specimenID") %>%
  filter(assay == "rnaSeq") %>%
  left_join(sample_meta, by = "individualID") %>%
  filter(specimenID %in% colnames(count_dat)) 

meta_dat <- meta_dat %>%
  mutate(apoe4 = factor(grepl("4", apoe_genotype), labels = c("APOE4-", "APOE4+")),
         apoe = factor(ifelse(apoe4 == "APOE4+", "pos", "neg"), levels = c("neg", "pos")),
         apoe_genotype = factor(apoe_genotype, levels = c("22", "23", "33", "24", "34", "44")),
         cogdx = recode_values(cogdx, 1 ~ "NCI", 2 ~ "MCI", 4 ~ "AD"), # diagnosis at death
         cogdx = factor(cogdx, levels = c("NCI", "MCI", "AD")),
         msex = factor(msex, levels = c("0", "1"), labels = c("Female", "Male")),
         age_death_90plus = factor(age_death == "90+", labels = c("No", "Yes")), # age at death (numerical + indicator vars)
         age_death_num = as.numeric(if_else(age_death == "90+", "90", age_death)),
         sequencingBatch = factor(sequencingBatch)) %>% 
  filter(!is.na(cogdx), # drop unknown diagnosis
         sequencingBatch != "0, 6, 7") %>% # drop 1 obs. with weird batch info
  mutate(sequencingBatch = droplevels(sequencingBatch),
         diag_apoe = factor(paste0(apoe4, " ", cogdx),
                            levels = c("APOE4- NCI", "APOE4- MCI", "APOE4- AD", 
                                       "APOE4+ NCI", "APOE4+ MCI", "APOE4+ AD")))

# length(unique(meta_dat$individualID)) --> total 2242 specimens from 1029 unique samples

tissue_labels <- c("dorsolateral prefrontal cortex" = "Dorsolateral Prefrontal Cortex",
                   "Head of caudate nucleus" = "Head of Caudate Nucleus",
                   "posterior cingulate cortex" = "Posterior Cingulate Cortex")

tissue_file_names <- c("dorsolateral prefrontal cortex" = "dlPFC",
                       "Head of caudate nucleus" = "hCN",
                       "posterior cingulate cortex" = "PCC")

## Split data by tissue ----
meta_by_tissue <- split(meta_dat, meta_dat$tissue)

count_by_tissue <- lapply(meta_by_tissue, function(meta) {
  count_dat[, meta$specimenID, drop = FALSE]
})

expr_by_tissue <- lapply(meta_by_tissue, function(meta) {
  expr_dat[, meta$specimenID, drop = FALSE]
})

## Demographic table ----
library(gtsummary)
demographic <- lapply(meta_by_tissue, function(meta) {
  meta %>%
  dplyr::select(apoe4, cogdx, msex, age_death_num, age_death_90plus, educ, pmi, RIN, apoe_genotype) %>%
  tbl_summary(
    by = apoe4,
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      age_death_num ~ "Age at Death (max 90)",
      age_death_90plus ~ "Age over 90",
      cogdx ~ "Diagnosis",
      educ ~ "Education (years)",
      msex ~ "Sex",
      pmi ~ "PMI",
      RIN ~ "RIN",
      apoe_genotype ~ "APOE genotype"
    ),
    missing = "no"
  ) %>%
  add_p() %>%         
  add_overall(last = TRUE) %>%   
  bold_labels()
})

# Variance Partition Analysis ---- 
## Analysis ----
library(BiocParallel)
param <- SnowParam(workers = parallel::detectCores() - 1, type = "SOCK")
register(param)

vp_form <- ~ (1|cogdx) + (1|apoe4) + (1|sequencingBatch) 

vp_by_tissue <- mapply(function(expr, meta) {
  rownames(meta) <- meta$specimenID
  expr <- expr[, meta$specimenID]
  variancePartition::fitExtractVarPartModel(expr, vp_form, meta, BPPARAM = param)
}, expr_by_tissue, meta_by_tissue, SIMPLIFY = FALSE)

## Result ----
factor_labels <- c(apoe4 = "APOE4 Status",
                   cogdx = "Diagnosis",
                   sequencingBatch = "Sequencing Batch",
                   Residuals= "Residuals")

vp_medians <- imap_dfr(vp_by_tissue, function(vp_fit, tissue_name) {
  as.data.frame(vp_fit) %>%
    summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
    pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
    mutate(Tissue = tissue_name,
           Factor = factor_labels[Factor],
           Median_VarExplained = round(Median_VarExplained, 4))
}) %>%
  arrange(Tissue, desc(Median_VarExplained))

# write.csv(vp_medians, file.path(table_output_dir, "RNASeq_varpart_medians.csv"), row.names = FALSE)

# DESeq2 ----
suppressMessages(library(DESeq2))
suppressMessages(library(biomaRt))

## Analysis ----
design_formula <- ~ cogdx + apoe + msex + age_death_90plus + sequencingBatch

min_group_size <- lapply(meta_by_tissue, function(meta) 
  min(table(meta$apoe)))

dds_by_tissue <- mapply(function(meta, count, min_n) {
  count <- count[, meta$specimenID]
  stopifnot(all(colnames(count) == meta$specimenID))
  dds <- DESeqDataSetFromMatrix(countData = count, 
                                colData = meta, 
                                design = design_formula)
  dds <- dds[rowSums(counts(dds) >= 10) >= min_n, ]  # filter: >=10 counts in at least min group size samples
  DESeq(dds)
}, meta_by_tissue, count_by_tissue, min_group_size, SIMPLIFY = FALSE)

## Results ----
p_val_threshold <- 0.05
logfc_threshold <- 0.26 

res_by_tissue <- lapply(dds_by_tissue, function(dds) {
  results(dds, contrast = c("apoe", "pos", "neg"), alpha = 0.05) %>% # use apoe4- as base level
    as.data.frame() %>%
    arrange(padj) %>%
    rownames_to_column("ensembl_gene_id") %>%
    mutate(negLog10FDR = -log10(padj),
           direction = case_when(padj < p_val_threshold & log2FoldChange >  logfc_threshold ~ "Upregulated",
                                 padj < p_val_threshold & log2FoldChange < -logfc_threshold ~ "Downregulated",
                                 TRUE ~ "Not significant"))
})

## Gene ID mapping ----
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

all_ensembl_ids <- lapply(res_by_tissue, `[[`, "ensembl_gene_id") %>%
  unlist() %>%
  unique()

id_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = "ensembl_gene_id",
                values = all_ensembl_ids,
                mart = mart) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  mutate(hgnc_symbol = if_else(hgnc_symbol == "", ensembl_gene_id, hgnc_symbol)) 

res_by_tissue <- lapply(res_by_tissue, function(res) 
  left_join(res, id_map, by = "ensembl_gene_id"))

sig_res_by_tissue <- lapply(res_by_tissue, function(res) 
  filter(res, direction != "Not significant")) # dlPFC has 0 sig gene; hCN has 7 and PCC has 4

## Volcano plots ----
make_volcano <- function(res, sig_res, tissue_name) {
  ggplot(res, aes(x = log2FoldChange, y = negLog10FDR)) +
    geom_point(aes(color = direction), alpha = 0.6, size = 1.2) +
    ggrepel::geom_text_repel(data = sig_res, aes(label = hgnc_symbol), 
                             size = 3, color = "black") +
    scale_color_manual(values = c("Upregulated" = "red", 
                                  "Downregulated" = "blue", 
                                  "Not significant" = "grey")) +
    geom_vline(xintercept = c(-logfc_threshold, logfc_threshold), 
               linetype = "dashed", color = "grey") +
    geom_hline(yintercept = -log10(p_val_threshold), 
               linetype="dashed", color="grey") +
    labs(title = paste("DE Genes by APOE4 Status in", tissue_labels[tissue_name]),
         subtitle = "Thresholds: adj_p < 0.05 & |log2 FC| > 0.26",
         x = "log2 fold change", y = "-log10 adjusted P", colour = "Direction of change") +
    theme_classic() +
    theme(legend.position = "bottom")
}

volcano_plots <- mapply(make_volcano,
                        res = res_by_tissue, 
                        sig_res = sig_res_by_tissue,
                        tissue_name = names(res_by_tissue), 
                        SIMPLIFY = FALSE)

## Save results ----
for (tissue in names(res_by_tissue)) {
  prefix <- tissue_file_names[tissue]
  write.csv(sig_res_by_tissue[[tissue]],
            file.path(table_output_dir, paste0(prefix, "_deseq2.csv")), row.names = FALSE)
  ggsave(file.path(pic_output_dir, paste0(prefix, "_deseq2_volcano.png")),
         plot = volcano_plots[[tissue]], width = 8, height = 6)
}

# Machine learning ----
## Batch correction ----
bc_expr_by_tissue <- mapply(function(expr, meta) {
  meta <- droplevels(meta)
  sva::ComBat(
    dat = as.matrix(expr[, meta$specimenID]),
    batch = meta$sequencingBatch,
    mod = model.matrix(~ apoe + cogdx + msex + age_death_90plus, data = meta)
  )
}, expr_by_tissue, meta_by_tissue, SIMPLIFY = FALSE)

# PCA before batch correction
pca_raw_by_tissue <- imap(expr_by_tissue, function(expr, tissue) {
  meta <- meta_by_tissue[[tissue]]

  pca <- prcomp(t(expr[, meta$specimenID]), center = TRUE, scale. = TRUE)
  as.data.frame(pca$x[, 1:2]) %>%
    bind_cols(meta %>% dplyr::select(sequencingBatch)) %>%
    ggplot(aes(PC1, PC2, color = sequencingBatch)) +
    geom_point() +
    labs(title = paste0(tissue_labels[tissue], ": PCA before batch correction"),
         color = "Sequencing batch") +
    theme_classic()
})

# PCA after batch correction
pca_bc_by_tissue <- imap(bc_expr_by_tissue, function(expr, tissue) {
  meta <- meta_by_tissue[[tissue]]
  
  pca <- prcomp(t(expr[, meta$specimenID]), center = TRUE, scale. = TRUE)
  as.data.frame(pca$x[, 1:2]) %>%
    bind_cols(meta %>% dplyr::select(sequencingBatch)) %>%
    ggplot(aes(PC1, PC2, color = sequencingBatch)) +
    geom_point() +
    labs(title = paste0(tissue_labels[tissue], ": PCA after batch correction"),
         color = "Sequencing batch") +
    theme_classic()
})

## Train/test split by tissue ----
library(caret)
set.seed(16)

split_by_tissue <- lapply(meta_by_tissue, function(meta) {
  train_idx <- createDataPartition(meta$apoe, p = 0.7, list = FALSE)
  list(train = meta[train_idx, ], test = meta[-train_idx, ])
})

train_meta_by_tissue <- lapply(split_by_tissue, `[[`, "train")
test_meta_by_tissue  <- lapply(split_by_tissue, `[[`, "test")

bc_train_by_tissue <- mapply(function(expr, meta) {
  expr[, meta$specimenID] %>% t() %>% as.data.frame() %>% mutate(apoe = meta$apoe)
}, bc_expr_by_tissue, train_meta_by_tissue, SIMPLIFY = FALSE)

bc_test_by_tissue <- mapply(function(expr, meta) {
  expr[, meta$specimenID] %>% t() %>% as.data.frame() %>% mutate(apoe = meta$apoe)
}, bc_expr_by_tissue, test_meta_by_tissue, SIMPLIFY = FALSE)


## Feature selection ----
# only calculate MI on 8000 proteins with highest variance
top_var_genes_by_tissue <- lapply(bc_train_by_tissue, function(df) {
  gene_cols  <- setdiff(colnames(df), "apoe")
  var_scores <- sapply(df[, gene_cols], var)
  names(sort(var_scores, decreasing = TRUE))[1:8000]
})

# MI calculation
mi_by_tissue <- imap(top_var_genes_by_tissue, function(genes, tissue) {
  df <- bc_train_by_tissue[[tissue]][, c(genes, "apoe")]
  result <- FSelectorRcpp::information_gain(apoe ~ ., data = df, type = "infogain") %>%
    arrange(desc(importance))
  colnames(result)[colnames(result) == "attributes"] <- "ensembl_gene_id"
  result
})

# Gene ID mapping
mart <- useEnsembl(biomart = "genes", dataset = "hsapiens_gene_ensembl")

all_ensembl_ids <- lapply(mi_by_tissue, `[[`, "ensembl_gene_id") %>%
  unlist() %>%
  unique()

id_map <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                filters = "ensembl_gene_id",
                values = all_ensembl_ids,
                mart = mart) %>%
  distinct(ensembl_gene_id, .keep_all = TRUE) %>%
  mutate(hgnc_symbol = if_else(hgnc_symbol == "", ensembl_gene_id, hgnc_symbol)) 

mi_by_tissue <- lapply(mi_by_tissue, function(res)
 left_join(res, id_map, by = "ensembl_gene_id"))

selected_genes_by_tissue <- lapply(mi_by_tissue, function(mi) 
  filter(mi, importance > 0)$hgnc_symbol)

# PCA check with MI > 0 proteins
pca_MI_by_tissue <- imap(selected_genes_by_tissue, function(genes, tissue) {
  meta <- train_meta_by_tissue[[tissue]]
  df <- bc_train_by_tissue[[tissue]]
  pca <- prcomp(df[, genes], center = TRUE, scale. = TRUE)
  
  as.data.frame(pca$x[, 1:2]) %>%
    bind_cols(meta %>% dplyr::select(apoe_genotype, apoe4, sequencingBatch, cogdx)) %>%
    ggplot(aes(x = PC1, y = PC2, color = apoe4)) +
    geom_point() +
    labs(title = sprintf("%s: PCA on %d MI-selected genes (Train, n = %d)",
                            tissue_labels[tissue], length(genes), nrow(df)),
         subtitle = "MI threshold > 0",
         color = "APOE status") +
    theme_classic()
})

## Model building ----
cv_ctrl <- trainControl(
  method = "repeatedcv",
  number = 5,
  repeats = 5,
  classProbs = TRUE,
  summaryFunction = twoClassSummary,
  savePredictions = "final"
)

set.seed(16)
rf_by_tissue <- imap(selected_genes_by_tissue, function(genes, tissue) {
  train_dat <- bc_train_by_tissue[[tissue]][, genes] %>%
    mutate(class = factor(bc_train_by_tissue[[tissue]]$apoe, levels = c("pos", "neg")))
  caret::train(class ~ .,
               data = train_dat,
               method = "rf",
               metric = "ROC",
               trControl = cv_ctrl,
               tuneLength = 10,
               importance = TRUE)
})

## Evaluation on held-out test set ----
test_results_by_tissue <- imap(rf_by_tissue, function(model, tissue) {
  genes <- selected_genes_by_tissue[[tissue]]
  test_dat <- bc_test_by_tissue[[tissue]][, genes] %>%
    mutate(class = factor(bc_test_by_tissue[[tissue]]$apoe, levels = c("pos", "neg")))
  
  pred_class <- predict(model, newdata = test_dat)
  pred_prob <- predict(model, newdata = test_dat, type = "prob")[, "pos"]
  cm <- confusionMatrix(pred_class, test_dat$class, positive = "pos")
  roc_obj <- pROC::roc(test_dat$class, pred_prob, levels = c("neg", "pos"))
  
  list(metrics = tibble(Tissue = tissue_labels[tissue],
                        N_train = nrow(model$trainingData),
                        N_test = nrow(test_dat),
                        N_features = length(genes),
                        Sensitivity = cm$byClass["Sensitivity"],
                        Specificity = cm$byClass["Specificity"],
                        PPV = cm$byClass["Pos Pred Value"],
                        NPV = cm$byClass["Neg Pred Value"],
                        AUC = as.numeric(roc_obj$auc),
                        CV_AUC = max(model$results$ROC),
                        AUC_gap = max(model$results$ROC) - as.numeric(roc_obj$auc)))
})

metrics_combined <- map_dfr(test_results_by_tissue, `[[`, "metrics")

## Save results ----
write.csv(metrics_combined, file.path(table_output_dir, "RNASeq_ML_metrics.csv"), row.names = FALSE)

for (tissue in names(rf_by_tissue)) {
  prefix <- tissue_file_names[tissue]
  write.csv(mi_by_tissue[[tissue]] %>% filter(importance > 0),
            file.path(table_output_dir, paste0(prefix, "_RNA_MI.csv")), row.names = FALSE)
  ggsave(file.path(pic_output_dir, paste0(prefix, "_MI_pca.png")),
         plot = pca_MI_by_tissue[[tissue]], width = 7, height = 6)
}


