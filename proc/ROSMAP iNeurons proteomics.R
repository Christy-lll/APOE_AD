# Set up ----
suppressMessages(library(gtsummary))
suppressMessages(library(here))
suppressMessages(library(limma))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(tidyverse))

pic_output_dir <- here("output", "ROSMAP", "pics")
table_output_dir <- here("output", "ROSMAP", "tables")

proteomic_dat <- read.csv(here("raw", "ROSMAP", "data_proteomics.csv"),
                          header = T, strip.white = T, na.strings = c("", "NA"), row.names = 1)
meta_dat <- read.csv(here("raw", "ROSMAP", "assay_meta_ROSMAP_iA_TMT-MS.csv"),
                     header = T, strip.white = T, na.strings = c("", "NA"))

# Preprocessing ----
## Match metadata to proteomics data ----
colnames(proteomic_dat) <- str_extract(colnames(proteomic_dat), "BR\\d+")
shared_BRID <- intersect(colnames(proteomic_dat), meta_dat$BRID)
proteomic_dat <- proteomic_dat[, shared_BRID] %>%
  mutate(across(everything(), ~ na_if(., 0))) 

meta_dat <- meta_dat %>%
  dplyr::select(apoe_genotype, sex, pmAD, BRID, SubjectID) %>% # borrow shared individual info from iAstrocytes meta
  mutate(
    apoe4 = factor(grepl("4", apoe_genotype), labels = c("APOE4-", "APOE4+")),
    apoe_genotype = factor(apoe_genotype, levels = c("22", "23", "33", "24", "34", "44")),
    sex = factor(sex, levels = c("f", "m"), labels = c("Female", "Male")),
    pmAD = factor(pmAD, levels = c("0", "1"), labels = c("NCI", "AD"))
  ) %>%
  filter(BRID %in% shared_BRID) %>%
  arrange(match(BRID, shared_BRID))

rownames(meta_dat) <- meta_dat$BRID
# all(colnames(proteomic_dat) == meta_dat$BRID)

## Impute & log transform & normalize proteomics data ----
expr_mat <- as.matrix(proteomic_dat) 

median_impute <- function(x) {
  x[is.na(x)] <- median(x, na.rm = TRUE)
  x
}

expr_mat <- expr_mat[rowMeans(is.na(expr_mat)) < 0.3, ] %>% # only keep proteins with <30% NA
  apply(1, median_impute) %>% # median imputation
  t() %>%
  log2() # log2 transformation 

expr_mat_norm <- sweep(expr_mat, 2, apply(expr_mat, 2, median), FUN = "-") # median normalisation


# Demographic ----
demographic <- meta_dat %>%
  dplyr::select(pmAD, sex, apoe_genotype) %>%
  tbl_summary(
    by = pmAD,
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


# Variance Partition Analysis ----
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

# write.csv(vp_medians, file.path(table_output_dir, "iNeurons_varpart.csv"), row.names = FALSE)


# Mutual Information Analysis ----
## Map Uniprot ID to Gene symbol ----
uniprot_ids <- rownames(expr_mat_norm)
gene_map <- AnnotationDbi::select(org.Hs.eg.db,
                                  keys = uniprot_ids,
                                  columns = c("SYMBOL"),
                                  keytype = "UNIPROT") %>%
  distinct(UNIPROT, .keep_all = TRUE)   

## Compute MI ----
mi_dat <- expr_mat_norm %>% t() %>% as.data.frame() %>% mutate(apoe4 = meta_dat$apoe4)

mi_scores_all <- FSelectorRcpp::information_gain(apoe4 ~ ., data = mi_dat, type = "infogain") %>%
  arrange(desc(importance)) %>% 
  left_join(gene_map, by = c("attributes" = "UNIPROT")) %>%
  dplyr::select(attributes, SYMBOL, importance) %>%
  mutate(label = if_else(is.na(SYMBOL), attributes, paste0(attributes, " | ", SYMBOL)))

mi_scores <- mi_scores_all %>% 
  filter(importance > 0)
  
# write.csv(mi_scores, file.path(table_output_dir, "iNeurons_MI.csv"), row.names = FALSE)


## Barplot ----
p_mi_bar <- ggplot(mi_scores %>% mutate(selected = importance > 0.24), 
                   aes(x = importance, y = forcats::fct_reorder(label, importance), fill = selected)) +
  geom_col(alpha = 0.9) +
  geom_vline(xintercept = 0.24, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey"),
                    name = "Selected for PCA") +
  labs(title = sprintf("%d Proteins with Mutual Information > 0", nrow(mi_scores)),
       subtitle = "13 proteins with MI > 0.24",
       x = "MI score",
       y = NULL) +
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05)) +
  theme_classic()

# ggsave(file.path(pic_output_dir, "iNeurons_MI_barplot.png"), p_mi_bar, width = 10, height = 8)

## PCA ----
top_proteins <- mi_scores_all %>%
  filter(importance > 0.24) %>%
  pull(attributes)

pca_MI <- prcomp(mi_dat[, top_proteins], center = TRUE, scale. = TRUE)
pca_MI_df <- as.data.frame(pca_MI$x[, 1:2]) %>%
  bind_cols(meta_dat %>% dplyr::select(apoe_genotype, apoe4, pmAD))

p_pca_MI <- ggplot(pca_MI_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA with %d MI > 0.24 proteins (n = %d)", 
                       length(top_proteins), nrow(mi_dat)),
       color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "iNeurons_MI_pca.png"), p_pca_MI, width = 7, height = 6)
