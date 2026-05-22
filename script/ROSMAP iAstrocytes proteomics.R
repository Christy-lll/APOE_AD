# Set up ----
suppressMessages(library(tidyverse))
suppressMessages(library(gtsummary))
suppressMessages(library(limma))
suppressMessages(library(here))

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
    apoe4 = factor(grepl("4", apoe_genotype), labels = c("APOE4-", "APOE4+")),
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

## QC checks ----
# anyNA(expr_mat) # no missing values
# boxplot(expr_mat, las = 2, cex.axis = 0.5, outline = FALSE) # data already log-transformed and centred 


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

# write.csv(vp_medians, file.path(table_output_dir, "iAstrocytes_varpart.csv"), row.names = FALSE)


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
  labs(title = sprintf("%d Proteins with Mutual Information > 0", nrow(mi_scores)),
       subtitle = "13 proteins with MI > 0.15",
       x = "MI score",
       y = NULL) +
  scale_x_continuous(breaks = seq(0, 0.3, by = 0.05)) +
  theme_classic()

# ggsave(file.path(pic_output_dir, "iAstrocytes_MI_barplot.png"), p_mi_bar, width = 8, height = 6)

## PCA ----
top_proteins <- mi_scores_all %>%
  filter(importance > 0.15) %>%
  pull(UniprotID)

pca_MI <- prcomp(mi_dat[, top_proteins], center = TRUE, scale. = TRUE)
pca_MI_df <- as.data.frame(pca_MI$x[, 1:2]) %>%
  bind_cols(meta_dat %>% dplyr::select(apoe_genotype, apoe4))

p_pca_MI <- ggplot(pca_MI_df, aes(x = PC1, y = PC2, color = apoe4)) +
  geom_point() +
  labs(title = sprintf("PCA with %d MI > 0.15 proteins (n = %d)",
                       length(top_proteins), nrow(mi_dat)),
       color = "") +
  theme_classic()

# ggsave(file.path(pic_output_dir, "iAstrocytes_MI_pca.png"), p_pca_MI, width = 7, height = 6)


# Enrichment Analysis ----
## Pull UniprotID for NetworkAnalyst analysis ----
background_universe <- mi_scores_all %>%
  distinct(UniprotID) %>%
  pull(UniprotID)

dep_gene_list <- mi_scores %>%
  distinct(UniprotID) %>%
  pull(UniprotID)

# writeLines(background_universe, file.path(table_output_dir, "iAstrocytes_networkanalyst_background.txt"))
# writeLines(dep_gene_list, file.path(table_output_dir, "iAstrocytes_networkanalyst_dep.txt"))

## KEGG results ----
kegg_res <- read.csv(file.path(table_output_dir, "iAstrocytes_KEGG.csv")) %>%
  filter(FDR < 0.05) %>%
  mutate(negLog10FDR = -log10(FDR),
         Pathway = forcats::fct_reorder(Pathway, negLog10FDR))

kegg_plot <- ggplot(kegg_res, aes(x = negLog10FDR, y = Pathway)) +
  geom_col(fill = "steelblue") +
  labs(title = "KEGG Pathway Enrichment",
       x = "-log10 FDR", y = NULL) +
  theme_classic()

# ggsave(file.path(pic_output_dir, "iAstrocytes_KEGG_barplot.png"), kegg_plot, width = 13, height = 13)

## PANTHER results ----
panther_res <- read.csv(file.path(table_output_dir, "iAstrocytes_PANTHER.csv")) %>%
  filter(FDR < 0.05) %>%
  mutate(negLog10FDR = -log10(FDR),
         Pathway = forcats::fct_reorder(Pathway, negLog10FDR))

panther_plot <- ggplot(panther_res, aes(x = negLog10FDR, y = Pathway)) +
  geom_col(fill = "indianred") +
  labs(title = "PANTHER Biological Processess Enrichment",
       x = "-log10 FDR", y = NULL) +
  theme_classic()

# ggsave(file.path(pic_output_dir, "iAstrocytes_PANTHER_barplot.png"), panther_plot, width = 8, height = 6)

## Reactome results ----
reactome_res <- read.csv(file.path(table_output_dir, "iAstrocytes_Reactome.csv")) %>%
  filter(FDR < 0.05) %>%
  mutate(negLog10FDR = -log10(FDR),
         Pathway = forcats::fct_reorder(Pathway, negLog10FDR))

reactome_plot <- ggplot(reactome_res, aes(x = negLog10FDR, y = Pathway)) +
  geom_col(fill = "#ABDDA4") +
  labs(title = "Reactome Enrichment",
       x = "-log10 FDR", y = NULL) +
  theme_classic()

# ggsave(file.path(pic_output_dir, "iAstrocytes_Reactome_barplot.png"), reactome_plot, width = 15, height = 13)
