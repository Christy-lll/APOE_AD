# Set up ----
suppressMessages(library(brant))
suppressMessages(library(gtsummary))
suppressMessages(library(here))
suppressMessages(library(tidyverse))
suppressMessages(library(MASS))

pic_output_dir <- here("output", "DiverseCohorts", "pics")
table_output_dir <- here("output", "DiverseCohorts", "tables")

clinical_dat <- read.csv(here("raw", "DiverseCohorts", "AMP-AD_DiverseCohorts_individual_metadata.csv"), 
                         header = T, strip.white = T, na.strings = c("", "NA", "missing or unknown"))

# Preprocessing ----
clinical_dat <- clinical_dat %>%
  filter(ADoutcome %in% c("Control", "AD")) %>%
  mutate(apoe4 = factor(apoeGenotype %in% c("24", "34", "44"), labels = c("APOE4-", "APOE4+")),
         apoeGenotype = factor(apoeGenotype, levels = c("22", "23", "33", "24", "34", "44")),
         ADoutcome = factor(ADoutcome, levels = c("Control", "AD")),
         ageDeath_num = as.numeric(if_else(ageDeath == "90+", "90", ageDeath)),
         sex = factor(sex),
         Braak = factor(Braak, 
                        levels = c("None", "Stage I", "Stage II", "Stage III", "Stage IV", "Stage V", "Stage VI")),
         amyCerad = factor(amyCerad, 
                           levels = c("None/No AD/C0", "Sparse/Possible/C1", "Moderate/Probable/C2", "Frequent/Definite/C3"))) %>%
  filter(!is.na(Braak), !is.na(amyCerad), !is.na(apoe4), !is.na(ageDeath_num), !is.na(sex)) %>%
  dplyr::select(individualID, apoe4, apoeGenotype, ADoutcome, sex, ageDeath_num, Braak, amyCerad)

# any(duplicated(clinical_dat$individualID))


# Demographic ----
demographic <- clinical_dat %>%
  dplyr::select(sex, ADoutcome, ageDeath_num, apoeGenotype) %>%
  tbl_summary(
    by = ADoutcome,
    statistic = list(
      all_continuous() ~ "{mean} ± {sd}",
      all_categorical() ~ "{n} ({p}%)"
    ),
    digits = all_continuous() ~ 1,
    label = list(
      sex ~ "Sex",
      ageDeath_num ~ "Age at Death (top-coded 90)",
      apoeGenotype ~ "APOE genotype"
    ),
    missing = "no"
  ) %>%
  add_overall(last = TRUE) %>%
  bold_labels()

# demographic %>% as_gt() %>% gt::gtsave(filename = file.path(table_output_dir, "clinical demographic table.png"))


# Visualisations ----
## Braak stage ----
braak_plot_dat <- clinical_dat %>%
  count(ADoutcome, apoe4, Braak) %>%
  group_by(ADoutcome, apoe4) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

braak <- ggplot(braak_plot_dat, aes(x = apoe4, y = prop, fill = Braak)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ ADoutcome) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = 1) +
  labs(title = "Braak Stage Distribution by APOE4 Status",
       x = "APOE4 Status", y = "Proportion", fill = "Braak Stage") +
  theme_classic() 

# ggsave(file.path(pic_output_dir, "clinical_braakstage_barplot.png"), braak, width = 7, height = 6)

## CERAD score ----
cerad_plot_dat <- clinical_dat %>%
  count(ADoutcome, apoe4, amyCerad) %>%
  group_by(ADoutcome, apoe4) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

cerad <- ggplot(cerad_plot_dat, aes(x = apoe4, y = prop, fill = amyCerad)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ ADoutcome) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = 5) +
  labs(title = "CERAD Score Distribution by APOE4 Status",
       x = "APOE4 Status", y = "Proportion", fill = "CERAD Score") +
  theme_classic() 

# ggsave(file.path(pic_output_dir, "clinical_ceradscore_barplot.png"), cerad, width = 7, height = 6)


# Regression Analysis ---- 
## Braak stage ----
polr_or_table <- function(fit, conf_level = 0.95) {
  n_coef <- length(coef(fit))
  ci <- confint(fit, level = conf_level)[1:n_coef, , drop = FALSE]
  exp(cbind(OR = coef(fit), lower = ci[, 1], upper = ci[, 2])) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term")
}

for (dx in c("Control", "AD")) {
  sub <- clinical_dat %>%
    filter(ADoutcome == dx) %>%
    droplevels() # only keep categories with observations
  
  if (nrow(sub) < 30) next
  
  fit <- polr(Braak ~ apoe4 + ageDeath_num + sex, # adjust for age and sex
              data = sub, Hess = TRUE, method = "logistic")
  
  message("--- Braak | group: ", dx, " (n=", nrow(fit$model), ", levels: ",
          paste(levels(sub$Braak), collapse = " < "), ") ---")
  print(polr_or_table(fit))
  
  tryCatch(
    print(brant(fit)),
    warning = function(w) message(" Brant test warning: ", conditionMessage(w))
  )
}

## CERAD score ----
glm_or_table <- function(fit, conf_level = 0.95) {
  ci <- confint(fit, level = conf_level)
  exp(cbind(OR = coef(fit), lower = ci[, 1], upper = ci[, 2])) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term") %>%
    filter(term != "(Intercept)")
}

for (dx in c("Control", "AD")) {
  sub <- clinical_dat %>%
    filter(ADoutcome == dx) %>%
    droplevels() %>% # only keep categories with observations
    mutate(amyCerad_bin = as.integer(amyCerad == levels(amyCerad)[2]))
  
  if (nrow(sub) < 30) next
  
  fit <- glm(amyCerad_bin ~ apoe4 + ageDeath_num + sex,
             data = sub, family = binomial(link = "logit"))
  
  message("--- CERAD | group: ", dx, " | binary logistic (n=", nrow(sub), ", 0 = ",
          levels(sub$amyCerad)[1], ", 1 = ", levels(sub$amyCerad)[2], ") ---")
  print(glm_or_table(fit))
}
