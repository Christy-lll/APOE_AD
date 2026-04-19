# Set up ----
suppressMessages(library(tidyverse))
library(here)

pic_output_dir <- here("output", "pics")
table_output_dir <- here("output", "tables")

proteomic_dat <- read.csv(here("raw", "ROSMAP-iA-TMT-MS.csv"), 
                          header = T, strip.white = T, na.strings = c("", "NA"), row.names = 1)
meta_dat <- read.csv(here("raw", "assay_meta_ROSMAP_iA_TMT-MS.csv"), 
                         header = T, strip.white = T, na.strings = c("", "NA"))

# Preprocessing ----
## Metadata ----
meta_dat <- meta_dat %>%
  mutate(
    apoe4 = factor(grepl("4", apoe_genotype), labels = c("APOE4-", "APOE4+")),
    apoe_genotype = factor(apoe_genotype, levels = c("22", "23", "33", "24", "34", "44")),
    sex = factor(sex, levels = c("f", "m"), labels = c("Female", "Male")),
    dx = factor(dx, levels = c("LPNCI", "HPNCI", "AD")),
    pmAD = factor(pmAD, levels = c("0", "1"), labels = c("NCI", "AD")),
    AD_apoe = factor(paste0(apoe4, " ", pmAD), levels = c("APOE4- NCI", "APOE4- AD", "APOE4+ NCI", "APOE4+ AD"))
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
colnames(expr_mat) <- meta_dat$SampleID
rownames(meta_dat) <- meta_dat$SampleID # all(colnames(expr_mat) == rownames(meta_dat))

# Log-transformation and normalisation check
anyNA(expr_mat) # no NA
summary(colMeans(expr_mat, na.rm = TRUE)) # comparable scale
boxplot(expr_mat, las = 2, cex.axis = 0.5, outline  = FALSE) # centered near 0

# Variance partition 
vp_form <- ~ (1 | batch) + (1 | apoe4) + (1 | pmAD) + (1 | sex) 
vp_fit <- variancePartition::fitExtractVarPartModel(expr_mat, vp_form, meta_dat)

factor_labels <- c(apoe4 = "APOE4 Status",
                   batch = "Batch",
                   pmAD = "Diagnosis",
                   Residuals = "Residuals",
                   sex = "Sex")

vp_medians <- as.data.frame(vp_fit) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
  mutate(Factor = factor_labels[Factor]) %>%
  arrange(desc(Median_VarExplained))

p_varpart <- variancePartition::plotVarPart(vp_fit, main = "Variance Partition: iAstrocytes Proteomics")
p_varpart$data <- p_varpart$data %>%
  mutate(variable = factor_labels[as.character(variable)])

# write.csv(vp_medians, file.path(table_output_dir, "iAstrocytes_varpart_medians.csv"), row.names = FALSE)
# ggsave(file.path(pic_output_dir, "iAstrocytes_varpart_violin.png"), p_varpart, width = 7, height = 5)

# Batch correction & PCA
design_bc <- model.matrix(~ apoe4 + sex + pmAD, data = meta_dat)
expr_mat_bc <- limma::removeBatchEffect(x = expr_mat, batch = meta_dat$batch, design = design_bc)

pca_before <- prcomp(t(expr_mat), center = TRUE, scale. = TRUE) # no significant batch effect
pca_before_df <- as.data.frame(pca_before$x[, 1:2]) %>%
  bind_cols(meta_dat %>% select(apoe_genotype, apoe4, batch))
ggplot(pca_before_df, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  labs(title = "Before batch correction", x = "PC1", y = "PC2", color = "") +
  theme_classic()

pca_after <- prcomp(t(expr_mat_bc), center = TRUE, scale. = TRUE) # no significant batch effect
pca_after_df <- as.data.frame(pca_after$x[, 1:2]) %>%
  bind_cols(meta_dat %>% select(apoe_genotype, apoe4, batch))
ggplot(pca_after_df, aes(x = PC1, y = PC2, color = batch)) +
  geom_point() +
  labs(title = "After batch correction", x = "PC1", y = "PC2", color = "") +
  theme_classic()

# Mutual Information Analysis ----
## Compute MI ----
mi_dat <- expr_mat_bc %>% t() %>% as.data.frame() %>% mutate(apoe4 = meta_dat$apoe4)

mi_scores_all <- FSelectorRcpp::information_gain(apoe4 ~ ., data = mi_dat, type = "infogain") %>%
  arrange(desc(importance)) %>%
  left_join(gene_symbols, by = c("attributes" = "UniprotID")) %>%
  select(attributes, GeneSymbol, importance) 

mi_scores <- mi_scores_all %>% 
  filter(importance > 0)

# write.csv(mi_scores, file.path(table_output_dir, "iAstrocytes_MI.csv"), row.names = FALSE)

## Barplot ----
p_mi_bar <- ggplot(mi_scores, aes(x = importance, 
                                  y = forcats::fct_reorder(GeneSymbol, importance))) +
  geom_col(width = 0.7, alpha = 0.9, fill = "steelblue") +
  labs(title= sprintf("%d Proteins with Mutual Information > 0", nrow(mi_scores)),
       x = "MI score",
       y = NULL) +
  theme_classic() +
  scale_x_continuous(limits = c(0, 0.25))

# ggsave(file.path(pic_output_dir, "iAstrocytes_MI_barplot.png"), p_mi_bar, width = 8, height = 6)

## PCA check ----
top_proteins <- mi_scores_all %>%
  slice_head(n = 13) %>% # mi > 0.15
  pull(attributes)

pca_MI <- prcomp(mi_dat[, top_proteins], center = TRUE, scale. = TRUE)
pca_MI_df <- as.data.frame(pca_MI$x[, 1:2]) %>%
  bind_cols(meta_dat %>% select(apoe_genotype, apoe4, dx))

p_pca_MI <- ggplot(pca_MI_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA with %d MI > 0.15 proteins (n = %d)", 
                       length(top_proteins), nrow(mi_dat)),
       x = "PC1", y = "PC2", color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "iAstrocytes_MI_pca.png"), p_pca_MI, width = 7, height = 6)