# Set up ----
suppressMessages(library(tidyverse))
library(here)

pic_output_dir <- here("output", "pics")
table_output_dir <- here("output", "tables")

proteomic_dat <- read.csv(here("raw", "data_proteomics.csv"), 
                          header = T, strip.white = T, na.strings = c("", "NA"), row.names = 1) 
iA_meta <- read.csv(here("raw", "assay_meta_ROSMAP_iA_TMT-MS.csv"), 
                     header = T, strip.white = T, na.strings = c("", "NA"))

# Preprocessing ----
## Match proteomic data and metadata ----
colnames(proteomic_dat) <- str_extract(colnames(proteomic_dat), "BR\\d+")
shared_BRID <- intersect(colnames(proteomic_dat), iA_meta$BRID)
proteomic_dat <- proteomic_dat[, shared_BRID] %>%
  mutate(across(everything(), ~ na_if(., 0))) # 2728 proteins x 31 samples

meta_dat <- iA_meta %>% # borrow sex, diagnosis, apoe_genotype from iA meta
  mutate(apoe4 = factor(grepl("4", apoe_genotype), labels = c("APOE4-", "APOE4+")),
         apoe_genotype = factor(apoe_genotype, levels = c("22", "23", "33", "24", "34", "44")),
         sex = factor(sex, levels = c("f", "m"), labels = c("Female", "Male")),
         dx = factor(dx, levels = c("LPNCI", "HPNCI", "AD")),
         pmAD = factor(pmAD, levels = c("0", "1"), labels = c("NCI", "AD")),
         AD_apoe = factor(paste0(apoe4, " ", pmAD), levels = c("APOE4- NCI", "APOE4- AD", "APOE4+ NCI", "APOE4+ AD"))) %>%
  filter(BRID %in% shared_BRID) %>%
  arrange(match(BRID, shared_BRID)) # all(colnames(proteomic_dat) == meta_dat$BRID)
rownames(meta_dat) <- meta_dat$BRID 

## Demographic table ----
library(gtsummary)
demographic <- meta_dat %>%
  select(AD_apoe, sex, apoe_genotype) %>%
  tbl_summary(
    by = AD_apoe,
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      sex ~ "Sex",
      apoe_genotype ~ "APOE genotype"
    ),
    missing = "no"
  ) %>%
  add_p() %>%         
  add_overall(last = TRUE) %>%   
  bold_labels()

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "iNeurons demographic table.png"))

## Process proteomic data ----
# Imputation & Log transformation
median_impute <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

expr_mat <- as.matrix(proteomic_dat) 
expr_mat <- expr_mat[rowMeans(is.na(expr_mat)) < 0.3, ] %>% # only keep proteins with <30% NA
  apply(1, median_impute) %>% # median imputation
  t() %>%
  log2() # log2 transformation 
# final dim: 1312 proteins x 31 samples

# Median Normalization
expr_mat_norm <- sweep(expr_mat, 2, apply(expr_mat, 2, median), FUN = "-")

boxplot(expr_mat, las = 2, col = "lavender", outline = FALSE, cex.axis = 0.8,
        main = "Before Normalisation", ylab = "log2 LFQ") # before
boxplot(expr_mat_norm, las = 2, col = "mistyrose", outline = FALSE, cex.axis = 0.8,
        main = "After Normalisation", ylab = "log2 LFQ (normalised)") # after
summary(colMeans(expr_mat_norm, na.rm = TRUE))

## Variance partition ----
vp_form <- ~ (1 | apoe4) + (1 | pmAD) + (1 | sex) 
vp_fit <- variancePartition::fitExtractVarPartModel(expr_mat_norm, vp_form, meta_dat)

factor_labels <- c(apoe4 = "APOE4 Status",
                   pmAD = "Diagnosis",
                   Residuals = "Residuals",
                   sex = "Sex")

vp_medians <- as.data.frame(vp_fit) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
  mutate(Factor = factor_labels[Factor]) %>%
  arrange(desc(Median_VarExplained))

p_varpart <- variancePartition::plotVarPart(vp_fit, main = "Variance Partition: iNeurons Proteomics")
p_varpart$data <- p_varpart$data %>%
  mutate(variable = factor_labels[as.character(variable)])

# write.csv(vp_medians, file.path(table_output_dir, "iNeurons_varpart_medians.csv"), row.names = FALSE)
# ggsave(file.path(pic_output_dir, "iNeurons_varpart_violin.png"), p_varpart, width = 7, height = 5)

# Mutual Information Analysis ----
## Map Uniprot ID to Gene symbol ----
library(org.Hs.eg.db)
uniprot_ids <- rownames(expr_mat_norm)
gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = uniprot_ids,
                                  columns = c("SYMBOL"),
                                  keytype = "UNIPROT") %>%
  distinct(UNIPROT, .keep_all = TRUE)   

## Compute MI ----
mi_dat <- expr_mat_norm %>% t() %>% as.data.frame() %>% mutate(apoe4 = meta_dat$apoe4)

mi_scores_all <- FSelectorRcpp::information_gain(apoe4 ~ ., data = mi_dat, type = "infogain") %>%
  arrange(desc(importance))

mi_scores <- mi_scores_all %>% 
  filter(importance > 0) %>% # 43 proteins with MI > 0 
  left_join(gene_map, by = c("attributes" = "UNIPROT")) %>%
  select(attributes, SYMBOL, importance) %>%
  mutate(label = if_else(is.na(SYMBOL), attributes, paste0(attributes, " | ", SYMBOL)))

# write.csv(mi_scores, file.path(table_output_dir, "iNeurons_MI.csv"), row.names = FALSE)

## Barplot ----
p_mi_bar <- ggplot(mi_scores, aes(x = importance, 
                                  y = forcats::fct_reorder(label, importance))) +
  geom_col(width = 0.7, alpha = 0.9, fill = "steelblue") +
  labs(title= sprintf("%d Proteins with Mutual Information > 0", nrow(mi_scores)),
       x = "MI score",
       y = NULL) +
  theme_classic() +
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05))

# ggsave(file.path(pic_output_dir, "iNeurons_MI_barplot.png"), p_mi_bar, width = 10, height = 8)

## PCA check ----
# all proteins
pca_all <- prcomp(t(expr_mat_norm), center = TRUE, scale. = TRUE) 
pca_all_df <- as.data.frame(pca_all$x[, 1:2]) %>%
  bind_cols(meta_dat %>% select(apoe_genotype, apoe4))
p_pca_all <- ggplot(pca_all_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = "PCA with all 1312 proteins (n = 31)", x = "PC1", y = "PC2", color = "") +
  theme_classic()

# MI selected proteins
top_proteins <- mi_scores_all %>%
  slice_head(n = 13) %>% # mi > 0.24
  pull(attributes)

pca_MI <- prcomp(mi_dat[, top_proteins], center = TRUE, scale. = TRUE)
pca_MI_df <- as.data.frame(pca_MI$x[, 1:2]) %>%
  bind_cols(meta_dat %>% select(apoe_genotype, apoe4, pmAD))

p_pca_MI <- ggplot(pca_MI_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA with %d MI > 0.24 proteins (n = %d)", 
                       length(top_proteins), nrow(mi_dat)),
       x = "PC1", y = "PC2", color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "iNeurons_all_pca.png"), p_pca_all, width = 7, height = 6)
# ggsave(file.path(pic_output_dir, "iNeurons_MI_pca.png"), p_pca_MI, width = 7, height = 6)
