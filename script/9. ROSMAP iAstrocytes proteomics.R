# Set up ----
suppressMessages(library(gtsummary))
suppressMessages(library(here))
suppressMessages(library(limma))
suppressMessages(library(org.Hs.eg.db))
suppressMessages(library(tidyverse))

pic_output_dir <- here("output", "ROSMAP", "pics")
table_output_dir <- here("output", "ROSMAP", "tables")

proteomic_dat <- read.csv(here("raw", "ROSMAP", "ROSMAP-iA-TMT-MS.csv"),
                          header = T, strip.white = T, na.strings = c("", "NA"), row.names = 1)
meta_dat <- read.csv(here("raw", "ROSMAP", "assay_meta_ROSMAP_iA_TMT-MS.csv"),
                     header = T, strip.white = T, na.strings = c("", "NA"))

# Preprocessing ----
## Metadata ----
meta_dat <- meta_dat %>%
  mutate(
    apoe4 = factor(grepl("4", apoe_genotype), labels = c("APOE ε4-", "APOE ε4+")),
    apoe_genotype = factor(apoe_genotype, levels = c("22", "23", "33", "24", "34", "44")),
    sex = factor(sex, levels = c("f", "m"), labels = c("Female", "Male")),
    pmAD = factor(pmAD, levels = c("0", "1"), labels = c("NCI", "AD"))
  )

## Proteomics data ----
gene_symbols <- proteomic_dat[, 1, drop = FALSE] %>%
  rownames_to_column("UniprotID") %>%
  dplyr::rename(GeneSymbol = symbol)

expr_mat <- as.matrix(proteomic_dat[, -1])
colnames(expr_mat) <- meta_dat$SampleID
rownames(meta_dat) <- meta_dat$SampleID

# no missing values; data already log-transformed and centred


# Demographic ----
demographic <- meta_dat %>%
  dplyr::select(sex, pmi, pmAD, apoe_genotype) %>%
  tbl_summary(
    by = pmAD,
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      sex ~ "Sex",
      pmi ~ "PMI",
      apoe_genotype ~ "APOE genotype"
    ),
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels()

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "iAstrocytes demographic table.png"))


# Variance Partition Analysis ----
vp_form <- ~ (1 | batch) + (1 | apoe4) + (1 | pmAD) + (1 | sex)
vp_fit <- variancePartition::fitExtractVarPartModel(expr_mat, vp_form, meta_dat)

factor_labels <- c(apoe4 = "APOE ε4 Status",
                   batch = "Batch",
                   pmAD = "Diagnosis",
                   Residuals = "Residuals",
                   sex = "Sex")

vp_medians <- as.data.frame(vp_fit) %>%
  summarise(across(everything(), ~ median(.x, na.rm = TRUE))) %>%
  pivot_longer(everything(), names_to = "Factor", values_to = "Median_VarExplained") %>%
  mutate(Factor = factor_labels[Factor]) %>%
  arrange(desc(Median_VarExplained))


# Limma ----
## Model fitting & results ----
p_val_threshold <- 0.05

design_mat <- model.matrix(~ apoe4 + pmAD + sex, data = meta_dat)
fit <- lmFit(expr_mat, design_mat) %>% eBayes()

res_limma <- topTable(fit, coef = "apoe4APOE ε4+", number = Inf) %>%
  rownames_to_column(var = "UniprotID") %>%
  left_join(gene_symbols, by = "UniprotID") %>%
  mutate(negLog10FDR = -log10(adj.P.Val),
         `Direction of change` = case_when(
           adj.P.Val < p_val_threshold & logFC > 0 ~ "Upregulated",
           adj.P.Val < p_val_threshold & logFC < 0 ~ "Downregulated",
           TRUE ~ "Not significant"))

sig_res_limma <- res_limma %>% filter(adj.P.Val < p_val_threshold)

## Volcano plot ----
n_up <- sum(sig_res_limma$`Direction of change` == "Upregulated")
n_down <- sum(sig_res_limma$`Direction of change` == "Downregulated")

p_limma_volcano <- ggplot(res_limma, aes(x = logFC, y = negLog10FDR)) +
  geom_point(aes(color = `Direction of change`), size = 1.2) +
  ggrepel::geom_text_repel(data = sig_res_limma,
                           aes(label = GeneSymbol),
                           size = 3, color = "black") +
  scale_color_manual(values = c("Upregulated" = "indianred",
                                "Downregulated" = "steelblue",
                                "Not significant" = "grey")) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey") +
  geom_hline(yintercept = -log10(p_val_threshold), linetype = "dashed", color = "grey") +
  labs(title = sprintf("%d differentially abundant iAstrocytes proteins by APOE ε4 carriage",
                       nrow(sig_res_limma)),
       subtitle = sprintf("Up: %d  Down: %d", n_up, n_down),
       x = "log2 fold change",
       y = "-log10 FDR",
       color = NULL) +
  theme_classic() +
  theme(legend.position = "bottom")

# ggsave(file.path(pic_output_dir, "iAstrocytes_limma_volcano.tiff"), p_limma_volcano, width = 7, height = 6, dpi = 300, compression = "lzw")


# Mutual Information Analysis ----
## Batch correction ----
design_bc <- model.matrix(~ apoe4 + sex + pmAD, data = meta_dat)
expr_mat_bc <- removeBatchEffect(x = expr_mat, batch = meta_dat$batch, design = design_bc)

## Compute MI ----
mi_dat <- expr_mat_bc %>% t() %>% as.data.frame() %>% mutate(apoe4 = meta_dat$apoe4)

mi_scores_all <- FSelectorRcpp::information_gain(apoe4 ~ ., data = mi_dat, type = "infogain") %>%
  arrange(desc(importance)) %>%
  left_join(gene_symbols, by = c("attributes" = "UniprotID")) %>%
  dplyr::rename(UniprotID = attributes) %>%
  dplyr::select(UniprotID, GeneSymbol, importance)

mi_scores <- mi_scores_all %>%
  filter(importance > 0)

# write.csv(mi_scores, file.path(table_output_dir, "iAstrocytes_MI.csv"), row.names = FALSE)

## Barplot ----
p_mi_bar <- ggplot(mi_scores %>% mutate(selected = importance > 0.15),
                   aes(x = importance, y = forcats::fct_reorder(GeneSymbol, importance), fill = selected)) +
  geom_col(alpha = 0.9) +
  geom_vline(xintercept = 0.15, linetype = "dashed", color = "black", linewidth = 0.6) +
  scale_fill_manual(values = c("TRUE" = "steelblue", "FALSE" = "grey"),
                    name = "Selected for PCA") +
  labs(title = sprintf("%d iAstrocytes proteins with non-zero mutual information with APOE ε4", nrow(mi_scores)),
       x = "MI score",
       y = NULL) +
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05)) +
  theme_classic() +
  theme(legend.position = "bottom")

# ggsave(file.path(pic_output_dir, "iAstrocytes_MI_barplot.tiff"), p_mi_bar, width = 7, height = 6, dpi = 300, compression = "lzw")

## PCA ----
top_proteins <- mi_scores_all %>%
  filter(importance > 0.15) %>%
  pull(UniprotID)

pca_MI <- prcomp(mi_dat[, top_proteins], center = TRUE, scale. = TRUE)
pca_MI_df <- as.data.frame(pca_MI$x[, 1:2]) %>%
  bind_cols(meta_dat %>% dplyr::select(apoe_genotype, apoe4))

p_pca_MI <- ggplot(pca_MI_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA of %d MI-selected iAstrocytes proteins by APOE ε4 carriage",
                       length(top_proteins)),
       color = "") +
  theme_classic() +
  theme(legend.position = "bottom")
  
# ggsave(file.path(pic_output_dir, "iAstrocytes_MI_pca.tiff"), p_pca_MI, width = 7, height = 6, dpi = 300, compression = "lzw")


# Enrichment Analysis ----
## Export gene lists for NetworkAnalyst ----
background_universe <- mi_scores_all %>% distinct(UniprotID) %>% pull(UniprotID)
dep_gene_list <- mi_scores %>% distinct(UniprotID) %>% pull(UniprotID)

# writeLines(background_universe, file.path(table_output_dir, "iAstrocytes_networkanalyst_background.txt"))
# writeLines(dep_gene_list, file.path(table_output_dir, "iAstrocytes_networkanalyst_dep.txt"))
