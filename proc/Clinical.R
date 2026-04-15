# Set up ----
suppressMessages(library(tidyverse))
library(here)

dat_output_dir <- here("cleaned_dat")
pic_output_dir <- here("output", "pics")
table_output_dir <- here("output", "tables")

clinical_dat <- read.csv(here("raw", "ROSMAP_clinical_1_.csv"), 
                         header = T, strip.white = T, na.strings = c("", "NA"))

# Preprocessing ----
clinical_dat <- clinical_dat %>%
  mutate(apoe4 = factor(grepl("4", apoe_genotype),
                        labels = c("APOE4-", "APOE4+")),
         braaksc = factor(braaksc,
                          levels = c("0", "1", "2", "3", "4", "5", "6"),
                          labels = c("0", "I", "II", "III", "IV", "V", "VI")),
         braaksc_grouped = case_when(braaksc %in% c("0", "I", "II") ~ "none/early",
                                     braaksc %in% c("III", "IV") ~ "moderate",
                                     braaksc %in% c("V", "VI") ~ "late",
                                     TRUE ~ NA_character_),
         braaksc_grouped = factor(braaksc_grouped, levels = c("none/early", "moderate", "late")),
         ceradsc = factor(ceradsc, 
                          levels = c("4", "3", "2", "1"),
                          labels = c("Not AD", "Possible", "Probable", "Definite")),
         cogdx = case_when(cogdx == "1" ~ "NCI", # diagnosis at death
                           cogdx == "2" ~ "MCI",
                           cogdx == "4" ~ "AD",
                           TRUE ~ NA_character_),
         dcfdx_lv = case_when(dcfdx_lv == "1" ~ "NCI", # diagnosis at last visit
                              dcfdx_lv == "2" ~ "MCI",
                              dcfdx_lv == "4" ~ "AD",
                              TRUE ~ NA_character_),
         msex = factor(msex,
                       levels = c("0", "1"),
                       labels = c("Female", "Male")),
         age_death_90plus = as.integer(age_death == "90+"), # age at death (numerical + indicator vars)
         age_death_num = as.numeric(if_else(age_death == "90+", "90", age_death)),
         age_lv_90plus = as.integer(age_at_visit_max == "90+"), # age at last visit (numerical + indicator vars)
         age_lv_num = as.numeric(if_else(age_at_visit_max == "90+", "90", age_at_visit_max)))
       

# Regression Analysis ---- 
library(MASS)
library(broom)

polr_or_table <- function(fit, conf_level = 0.95) {
  n_coef <- length(coef(fit))
  ci <- confint(fit, level = conf_level)[1:n_coef, , drop = FALSE]
  exp(cbind(OR = coef(fit),
            lower = ci[, 1],
            upper = ci[, 2])) %>%
    as.data.frame() %>%
    tibble::rownames_to_column("term")
}

## Braak stage ----
for (dx in c("NCI", "MCI", "AD")) {
  sub <- filter(clinical_dat, cogdx == dx)
  if (nrow(sub) < 30) next
  
  fit <- polr(braaksc_grouped ~ apoe4 + age_death_num + age_death_90plus + msex,
              data = sub, Hess = TRUE, method = "logistic")
  
  n_model <- nrow(fit$model)
  message("--- group: ", dx, " (n=", n_model, ") ---")
  print(polr_or_table(fit))
  tryCatch(
    print(brant::brant(fit)),
    warning = function(w) message("Brant test warning: ", conditionMessage(w))
  )
}

## CERAD score ----
for (dx in c("NCI", "MCI", "AD")) {
  sub <- filter(clinical_dat, cogdx == dx)
  if (nrow(sub) < 30) next
  
  fit <- polr(ceradsc ~ apoe4 + age_death_num + age_death_90plus + msex,
              data = sub, Hess = TRUE, method = "logistic")
  
  n_model <- nrow(fit$model)
  message("--- group: ", dx, " (n=", n_model, ") ---")
  print(polr_or_table(fit))
  print(brant::brant(fit))
}

## MMSE score ----
for (dx in c("NCI", "MCI", "AD")) {
  sub <- filter(clinical_dat, cogdx == dx)
  if (nrow(sub) < 30) next
  
  fit <- lm(cts_mmse30_lv ~ apoe4 + age_lv_num + age_lv_90plus + msex,
            data = sub)
  
  n_model <- nrow(fit$model)
  message("--- group: ", dx, " (n=", n_model, ") ---")
  
  print(summary(fit)$coefficients)
  print(confint(fit))
}

# Traditional Statistical Tests ----
## Braak stage - Chi-square test ----
for (dx in c("NCI", "MCI", "AD")) {
  sub <- filter(clinical_dat, cogdx == dx) %>%
    filter(!is.na(braaksc_grouped), !is.na(apoe4))
  if (nrow(sub) < 30) next
  
  message("--- group: ", dx, " (n=", nrow(sub), ") ---")
  tbl <- table(sub$apoe4, sub$braaksc_grouped)
  print(tbl)
  print(chisq.test(tbl))
}

## CERAD score - Chi-square test ----
for (dx in c("NCI", "MCI", "AD")) {
  sub <- filter(clinical_dat, cogdx == dx) %>%
    filter(!is.na(ceradsc), !is.na(apoe4))
  if (nrow(sub) < 30) next
  
  message("--- group: ", dx, " (n=", nrow(sub), ") ---")
  tbl <- table(sub$apoe4, sub$ceradsc)
  print(tbl)
  print(chisq.test(tbl))
}

## MMSE score - Mann-Whitney test ----
hist(clinical_dat$cts_mmse30_lv, breaks = 30)

for (dx in c("NCI", "MCI", "AD")) {
  sub <- filter(clinical_dat, cogdx == dx) %>%
    filter(!is.na(cts_mmse30_lv), !is.na(apoe4))
  if (nrow(sub) < 30) next
  
  message("--- group: ", dx, " (n=", nrow(sub), ") ---")
  print(wilcox.test(cts_mmse30_lv ~ apoe4, data = sub))
}

# Visual analysis ----
## Braak stage ----
braak_plot_dat <- clinical_dat %>%
  filter(!is.na(braaksc_grouped), !is.na(apoe4), !is.na(cogdx)) %>%
  count(cogdx, apoe4, braaksc_grouped) %>%
  group_by(cogdx, apoe4) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(braak_plot_dat, aes(x = apoe4, y = prop, fill = braaksc_grouped)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ cogdx) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = 1) +
  labs(title = "Braak Stage Distribution by APOE4 Status",
       x = "APOE4 Status", y = "Proportion", fill = "Braak Stage") +
  theme_classic() 

## CERAD score ----
cerad_plot_dat <- clinical_dat %>%
  filter(!is.na(ceradsc), !is.na(apoe4), !is.na(cogdx)) %>%
  count(cogdx, apoe4, ceradsc) %>%
  group_by(cogdx, apoe4) %>%
  mutate(prop = n / sum(n)) %>%
  ungroup()

ggplot(cerad_plot_dat, aes(x = apoe4, y = prop, fill = ceradsc)) +
  geom_bar(stat = "identity", position = "fill") +
  facet_wrap(~ cogdx) +
  scale_y_continuous(labels = scales::percent) +
  scale_fill_brewer(palette = 5) +
  labs(title = "CERAD Score Distribution by APOE4 Status",
       x = "APOE4 Status", y = "Proportion", fill = "CERAD Score") +
  theme_classic() 

## MMSE scores ----
mmse_plot_dat <- clinical_dat %>%
  filter(!is.na(cts_mmse30_lv), !is.na(apoe4), !is.na(cogdx))

ggplot(mmse_plot_dat, aes(x = apoe4, y = cts_mmse30_lv, fill = apoe4)) +
  geom_boxplot(width = 0.15, alpha = 0.8) +
  facet_wrap(~ cogdx) +
  scale_fill_manual(values = c("APOE4-" = "#4393c3",
                               "APOE4+" = "#d6604d"),
                    guide = "none") +
  labs(title = "MMSE Score by APOE4 Status",
       x = "APOE4 Status", y = "MMSE Score (last visit)") +
  theme_classic() 