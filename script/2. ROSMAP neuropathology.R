# Set up ----
suppressMessages(library(brant))
suppressMessages(library(gtsummary))
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(MASS))

pic_output_dir <- here("output", "ROSMAP", "pics")
table_output_dir <- here("output", "ROSMAP", "tables")

clinical_dat <- read.csv(here("raw", "ROSMAP", "ROSMAP_clinical_1_.csv"),
                         header = T, strip.white = T, na.strings = c("", "NA"))

# Preprocessing ----
clinical_dat <- clinical_dat %>%
  filter(!is.na(apoe_genotype)) %>%
  mutate(apoe4 = factor(apoe_genotype %in% c("24", "34", "44"), labels = c("APOE ε4-", "APOE ε4+")),
         apoe_genotype = factor(apoe_genotype, levels = c("22", "23", "33", "24", "34", "44")),
         braaksc = factor(braaksc,
                          levels = c("0", "1", "2", "3", "4", "5", "6"),
                          labels = c("None", "Stage I", "Stage II", "Stage III", "Stage IV", "Stage V", "Stage VI")),
         braaksc_grouped = factor(
           case_when(braaksc %in% c("None", "Stage I", "Stage II") ~ "None-Stage II",
                     braaksc %in% c("Stage III", "Stage IV") ~ "Stages III-IV",
                     braaksc %in% c("Stage V", "Stage VI") ~ "Stages V-VI",
                     TRUE ~ NA_character_),
           levels = c("None-Stage II", "Stages III-IV", "Stages V-VI")),
         ceradsc = factor(ceradsc,
                          levels = c("4", "3", "2", "1"),
                          labels = c("None/No AD", "Sparse/Possible", "Moderate/Probable", "Frequent/Definite")),
         cogdx = factor(
           case_when(cogdx == "1" ~ "NCI",
                     cogdx == "2" ~ "MCI",
                     cogdx == "4" ~ "AD",
                     TRUE ~ NA_character_),
           levels = c("NCI", "MCI", "AD")),
         msex = factor(msex,
                       levels = c("0", "1"),
                       labels = c("Female", "Male")),
         race = factor(race,
                       levels = 1:7,
                       labels = c("White", "Black or African American", "American Indian or Alaska Native",
                                  "Native Hawaiian or Other Pacific Islander", "Asian", "Other", "Unknown")),
         age_death_num = as.numeric(if_else(age_death == "90+", "90", age_death))) %>%
  filter(!is.na(braaksc_grouped), !is.na(ceradsc), !is.na(cogdx), !is.na(msex), !is.na(age_death_num)) %>%
  dplyr::select(projid, apoe4, apoe_genotype, braaksc, braaksc_grouped, ceradsc, cogdx, msex, age_death_num, educ, pmi, race)

any(duplicated(clinical_dat$projid))


# Demographic ----
demographic <- clinical_dat %>%
  dplyr::select(cogdx, age_death_num, msex, apoe_genotype, race,pmi) %>%
  tbl_summary(
    by = cogdx,
    statistic = list(
      all_continuous() ~ "{mean} \u00b1 {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      msex ~ "Sex",
      age_death_num ~ "Age at Death (top-coded 90)",
      apoe_genotype ~ "APOE genotype",
      pmi ~ "Post Mortem Interval (hours)",
      race ~ "Race"
    ),
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels()


# Visualisations ----
## Braak stage ----
braak_plot_dat <- clinical_dat %>%
  count(cogdx, apoe4, braaksc_grouped) %>%
  group_by(cogdx, apoe4) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

braak <- ggplot(braak_plot_dat, aes(x = apoe4, y = prop, fill = braaksc_grouped)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ cogdx) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = 1) +
  labs(title = "Braak stage distribution in ROSMAP cohort",
       x = NULL, y = "Proportion", fill = "Braak Stage") +
  theme_classic()

ggsave(file.path(pic_output_dir, "clinical_braakstage_barplot.tiff"), braak, width = 7, height = 6, dpi = 300, compression = "lzw")

## CERAD score ----
cerad_plot_dat <- clinical_dat %>%
  count(cogdx, apoe4, ceradsc) %>%
  group_by(cogdx, apoe4) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

cerad <- ggplot(cerad_plot_dat, aes(x = apoe4, y = prop, fill = ceradsc)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ cogdx) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = 5) +
  labs(title = "CERAD score distribution in ROSMAP cohort",
       x = NULL, y = "Proportion", fill = "CERAD Score") +
  theme_classic()

ggsave(file.path(pic_output_dir, "clinical_ceradscore_barplot.tiff"), cerad, width = 7, height = 6, dpi = 300, compression = "lzw")


# Regression Analysis ----
polr_or_table <- function(fit, conf_level = 0.95) {
  n_coef <- length(coef(fit))
  ci <- confint(fit, level = conf_level)[1:n_coef, , drop = FALSE]
  exp(cbind(OR = coef(fit), lower = ci[, 1], upper = ci[, 2])) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term")
}

## Braak stage ----
for (dx in c("NCI", "MCI", "AD")) {
  sub <- filter(clinical_dat, cogdx == dx)
  if (nrow(sub) < 30) next
  
  fit <- polr(braaksc_grouped ~ apoe4 + age_death_num + msex,
              data = sub, Hess = TRUE, method = "logistic")
  
  message("--- Braak | group: ", dx, " (n=", nrow(fit$model), ") ---")
  print(polr_or_table(fit))
  
  tryCatch(
    print(brant(fit)),
    warning = function(w) message(" Brant test warning: ", conditionMessage(w))
  )
}

## CERAD score ----
for (dx in c("NCI", "MCI", "AD")) {
  sub <- filter(clinical_dat, cogdx == dx)
  if (nrow(sub) < 30) next
  
  fit <- polr(ceradsc ~ apoe4 + age_death_num + msex,
              data = sub, Hess = TRUE, method = "logistic")
  
  message("--- CERAD | group: ", dx, " (n=", nrow(fit$model), ") ---")
  print(polr_or_table(fit))
  
  tryCatch(
    print(brant(fit)),
    warning = function(w) message(" Brant test warning: ", conditionMessage(w))
  )
}