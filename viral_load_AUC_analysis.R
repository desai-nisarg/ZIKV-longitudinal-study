###############################################################################
## Viral load (AUC) analysis
##
## This script includes examination of (a) SEX DIFFERENCES in
## viral replication and (b) ASSOCIATION OF VIREMIA with the behavioral /
## neuroimaging outcomes.
##
## We use the area under the viremia curve (AUC) as a single
## per-animal summary of viral-load BURDEN (magnitude x duration). AUC exists
## only for the 12 ZIKV-infected infants (6 F / 6 M), so every analysis here is
## WITHIN the ZIKV group. All modelling follows the manuscript's Bayesian
## conventions: brms, weakly-informative priors, and posterior summaries
## reported as the median with the 89% highest-density interval (HDI) and the
## probability of direction (pd). 
##
## Two analyses:
##   A. Sex difference in viral burden:      log(AUC) ~ Sex
##   B. Viremia-outcome dose-response (within ZIKV): outcome ~ AUC_z, in the
##      SAME model family each outcome uses in the manuscript.
##         - Neuroimaging  : Gaussian, ICV-adjusted, outcome standardized
##         - Acute stress  : negative binomial (log-rate)
##         - Attachment    : negative binomial (log-rate)
##         - INAS domains  : ordinal cumulative-logit over age (weeks 2-9)
##
## Targeted to the principal findings only, to limit multiplicity at n = 12. 
## SEX HANDLING: within-sex for the sex-specific outcomes, pooled otherwise.
##
## AUC is highly right-skewed and one male is an extreme
## low outlier (~250k, >2 log10 units below the next animal), so everything is
## on the log scale and every model is refit EXCLUDING that animal as a
## sensitivity check.
###############################################################################

## ------------------------------------------------------------------ set-up ---
if (!requireNamespace("pacman", quietly = TRUE)) install.packages("pacman")
pacman::p_load(brms, dplyr, tidyr, readr, readxl, janitor, stringr, purrr, tibble)

set.seed(111)

## Sampler settings (small n -> higher adapt_delta for stability)
CHAINS      <- 4
ITER        <- 4000
WARMUP      <- 2000
ADAPT_DELTA <- 0.95
CORES       <- 4

## File locations (run with this repo folder as the working directory) --------
AUC_XLSX   <- "ZIKV-R01_LongitudinalStudy_animals&sex_8-14-26.xlsx"  # place the xlsx here
MRI_CSV    <- "ZIKV_3moMRI.csv"
HI_CSV     <- "human_intruder.csv"
ATTACH_CSV <- "attachment.csv"
INAS_CSV   <- "SNIP_cleaned.csv"

## ------------------------------------------------------- helper functions ---

## 89% highest-density interval from a posterior draw vector
hdi89 <- function(x, prob = 0.89) {
  x <- sort(x[is.finite(x)])
  n <- length(x)
  if (n < 2) return(c(lower = NA_real_, upper = NA_real_))
  m <- max(1, floor(prob * n))
  widths <- x[(m):n] - x[1:(n - m + 1)]
  i <- which.min(widths)
  c(lower = x[i], upper = x[i + m - 1])
}

## probability of direction: posterior mass on the dominant side of 0
pdir <- function(x) {
  x <- x[is.finite(x)]
  max(mean(x > 0), mean(x < 0))
}

## Summarize a posterior draw vector into a one-row tibble
summarise_effect <- function(draws, analysis, outcome, subset_label, term,
                             scale_note) {
  hd <- hdi89(draws)
  tibble(
    Analysis   = analysis,
    Outcome    = outcome,
    Subset     = subset_label,
    Term       = term,
    Scale      = scale_note,
    Median     = round(median(draws), 3),
    `HDI 2.5`  = round(unname(hd["lower"]), 3),  # 89% HDI lower
    `HDI 97.5` = round(unname(hd["upper"]), 3),  # 89% HDI upper
    pd         = round(pdir(draws), 3),
    Credible   = ifelse(hd["lower"] > 0 | hd["upper"] < 0, "Yes", "No")
  )
}

## Pull a single population-level coefficient's draws by name
coef_draws <- function(fit, varname) {
  d <- as.data.frame(as_draws_df(fit))
  if (!varname %in% names(d)) {
    stop(sprintf("Coefficient '%s' not found. Available: %s",
                 varname, paste(names(d), collapse = ", ")))
  }
  d[[varname]]
}

## Weakly-informative prior sets by family ------------------------------------
priors_gaussian_std <- c(
  prior(normal(0, 1), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "sigma")
)
priors_negbin <- c(
  prior(normal(0, 5), class = "Intercept"),
  prior(normal(0, 1), class = "b"),
  prior(exponential(1), class = "shape")
)
priors_ordinal <- c(
  prior(normal(0, 1), class = "b")
)

## ------------------------------------------------------------- load AUC -----
auc_raw <- readxl::read_excel(AUC_XLSX, sheet = 1) %>% clean_names()
nm <- names(auc_raw)
pub_col    <- nm[str_detect(nm, "publication|pub")][1]
auc_col    <- nm[str_detect(nm, "area|auc|curve")][1]
sex_col    <- nm[str_detect(nm, "^sex")][1]
trt_col    <- nm[str_detect(nm, "treatment|group")][1]
cohort_col <- nm[str_detect(nm, "cohort")][1]

auc <- auc_raw %>%
  transmute(
    pub    = as.integer(.data[[pub_col]]),
    auc    = as.numeric(.data[[auc_col]]),
    sex    = factor(ifelse(as.numeric(.data[[sex_col]]) == 1, "Male", "Female"),
                    levels = c("Female", "Male")),
    trt    = as.numeric(.data[[trt_col]]),
    cohort = if (!is.na(cohort_col)) factor(.data[[cohort_col]]) else factor(NA)
  ) %>%
  filter(trt == 2, !is.na(auc)) %>%          # ZIKV only
  mutate(
    log_auc     = log(auc),
    outlier_low = pub == 2                    # Publication order 2 = extreme low
  )

## Global z-score of log(AUC): one unit = 1 SD of log burden among infected
auc$auc_z <- as.numeric(scale(auc$log_auc))
## Sensitivity: z recomputed among the 11 animals retained after dropping outlier
auc$auc_z_no_out <- NA_real_
keep <- !auc$outlier_low
auc$auc_z_no_out[keep] <- as.numeric(scale(auc$log_auc[keep]))

cat("\n=== ZIKV AUC data (n =", nrow(auc), ") ===\n")
print(auc %>% select(pub, sex, auc, log_auc, auc_z, outlier_low))

## Merge helper: attach AUC columns to an outcome data frame keyed on `pub`
attach_auc <- function(dat, pub_field) {
  dat$pub <- as.integer(dat[[pub_field]])
  dat %>% left_join(auc %>% select(pub, auc_z, auc_z_no_out, auc_sex = sex),
                    by = "pub")
}

## Generic runner: fit main + outlier-excluded, summarize one coefficient ------
run_assoc <- function(build_data, formula, family, priors, coef_name,
                      analysis, outcome, subset_label, scale_note) {
  out <- list()
  for (mode in c("main", "no_outlier")) {
    d <- build_data(mode)
    if (is.null(d) || nrow(d) < 4) next
    fit <- tryCatch(
      brm(formula = formula, data = d, family = family, prior = priors,
          chains = CHAINS, iter = ITER, warmup = WARMUP, cores = CORES,
          seed = 123, control = list(adapt_delta = ADAPT_DELTA),
          refresh = 0, silent = 2),
      error = function(e) { message("  fit failed (", mode, "): ", e$message); NULL })
    if (is.null(fit)) next
    draws <- coef_draws(fit, coef_name)
    lab <- if (mode == "main") subset_label else paste0(subset_label, " [excl. outlier]")
    out[[mode]] <- summarise_effect(draws, analysis, outcome, lab,
                                    coef_name, scale_note)
  }
  bind_rows(out)
}

results <- list()

###############################################################################
## ANALYSIS A -- Sex difference in viral burden
##   log(AUC) ~ Sex          (Gaussian on log scale)
##   b_sexMale = male - female difference in log(AUC); exp() = fold-change.
###############################################################################
cat("\n\n################ ANALYSIS A: sex difference in viral burden ########\n")

priors_A <- c(
  prior(normal(0, 5),  class = "Intercept"),
  prior(normal(0, 2),  class = "b"),          # weakly-informative on log scale
  prior(exponential(1), class = "sigma")
)

runA <- function(mode) {
  d <- if (mode == "main") auc else auc %>% filter(!outlier_low)
  d
}
for (mode in c("main", "no_outlier")) {
  d <- runA(mode)
  fitA <- brm(log_auc ~ sex, data = d, family = gaussian(), prior = priors_A,
              chains = CHAINS, iter = ITER, warmup = WARMUP, cores = CORES,
              seed = 123, control = list(adapt_delta = ADAPT_DELTA),
              refresh = 0, silent = 2)
  dr <- coef_draws(fitA, "b_sexMale")
  lab <- if (mode == "main") "All ZIKV (6F/6M)" else "ZIKV [excl. outlier] (6F/5M)"
  results[[paste0("A_", mode)]] <-
    summarise_effect(dr, "A. Sex diff in burden", "log(AUC)", lab,
                     "b_sexMale (M - F)", "log scale; exp() = fold-change")
  ## Also report the fold-change directly
  cat(sprintf("[%s] Male vs Female fold-change in AUC = %.2f  (89%% HDI %.2f, %.2f), pd = %.3f\n",
              lab, exp(median(dr)), exp(hdi89(dr)["lower"]),
              exp(hdi89(dr)["upper"]), pdir(dr)))
}

## Optional cohort-adjusted sensitivity (only if cohort present & varies)
if (!all(is.na(auc$cohort)) && dplyr::n_distinct(auc$cohort) > 1) {
  fitA_c <- brm(log_auc ~ sex + cohort, data = auc, family = gaussian(),
                prior = priors_A, chains = CHAINS, iter = ITER, warmup = WARMUP,
                cores = CORES, seed = 123,
                control = list(adapt_delta = ADAPT_DELTA), refresh = 0, silent = 2)
  dr <- coef_draws(fitA_c, "b_sexMale")
  results[["A_cohort"]] <-
    summarise_effect(dr, "A. Sex diff in burden", "log(AUC)",
                     "All ZIKV, cohort-adjusted", "b_sexMale (M - F)",
                     "log scale; exp() = fold-change")
}

###############################################################################
## ANALYSIS B1 -- Neuroimaging (Gaussian, ICV-adjusted, standardized outcome)
##   Within the relevant sex for each sex-specific region:
##     scale(volume) ~ auc_z + scale(ICV)
##   b_auc_z = standardized change in regional volume per 1 SD of log(AUC).
###############################################################################
cat("\n\n################ ANALYSIS B1: neuroimaging associations ############\n")

mri <- read_csv(MRI_CSV, show_col_types = FALSE) %>% clean_names()
mri_sexcol <- names(mri)[str_detect(names(mri), "^sex")][1]
mri_pubcol <- names(mri)[str_detect(names(mri), "pub")][1]
mri_icvcol <- names(mri)[str_detect(names(mri), "intracranial")][1]
mri <- mri %>%
  mutate(sex = factor(ifelse(.data[[mri_sexcol]] == 1, "Male", "Female"),
                      levels = c("Female", "Male"))) %>%
  attach_auc(mri_pubcol) %>%
  filter(!is.na(auc_z))                       # ZIKV animals only (have AUC)

## principal neuroimaging findings: region clean-name, label, relevant sex
neuro_targets <- tribble(
  ~col,                        ~label,                       ~which_sex,
  "total_amygdala_mm3",        "Amygdala",                   "Male",
  "total_csf_volume_mm3",      "Total CSF",                  "Female",
  "lateral_ventricle_mm3",     "Lateral ventricle",          "Female",
  "rt_temporal_limbic_gm_wm",  "Right temporal-limbic",      "Female",
  "lt_temporal_limbic_gm_wm",  "Left temporal-limbic",       "Female",
  "lt_temporal_visual_gm_wm",  "Left temporal-visual",       "Male"
)

for (i in seq_len(nrow(neuro_targets))) {
  tg <- neuro_targets[i, ]
  if (!tg$col %in% names(mri)) { message("  missing region col: ", tg$col); next }
  cat("\n--- Neuro:", tg$label, "(", tg$which_sex, ") ---\n")

  build <- function(mode) {
    z <- if (mode == "main") "auc_z" else "auc_z_no_out"
    d <- mri %>%
      filter(sex == tg$which_sex, !is.na(.data[[tg$col]]), !is.na(.data[[z]])) %>%
      transmute(
        y     = as.numeric(scale(.data[[tg$col]])),
        auc_z = .data[[z]],
        icv_z = as.numeric(scale(.data[[mri_icvcol]]))
      )
    if (nrow(d) < 4) return(NULL)
    d
  }
  res <- run_assoc(build, bf(y ~ auc_z + icv_z), gaussian(), priors_gaussian_std,
                   "b_auc_z", "B. Viremia-outcome", tg$label,
                   paste0(tg$which_sex, " ZIKV, ICV-adj"),
                   "std. slope (per SD log-AUC)")
  results[[paste0("B_neuro_", tg$col)]] <- res
}

###############################################################################
## ANALYSIS B2 -- Acute stress (negative binomial, log-rate)
##   Hostile behavior (Alone, Profile conditions), pooled across sex
##   Female threat barks (within females, condition-adjusted)
###############################################################################
cat("\n\n################ ANALYSIS B2: acute-stress associations ############\n")

hi <- read_csv(HI_CSV, show_col_types = FALSE) %>% clean_names()
hi_sexcol  <- names(hi)[str_detect(names(hi), "^sex")][1]
hi_pubcol  <- names(hi)[str_detect(names(hi), "pub")][1]
hi_condcol <- names(hi)[str_detect(names(hi), "condition")][1]
hi_hostcol <- names(hi)[str_detect(names(hi), "hostile")][1]
hi_barkcol <- names(hi)[str_detect(names(hi), "threat_barks")][1]

hi <- hi %>%
  mutate(
    sex        = factor(ifelse(.data[[hi_sexcol]] == 1, "Male", "Female"),
                        levels = c("Female", "Male")),
    condition  = factor(.data[[hi_condcol]], levels = c(0, 1, 2),
                        labels = c("Alone", "Profile", "Stare")),
    animal     = factor(.data[[hi_pubcol]]),
    hostile    = as.integer(round(.data[[hi_hostcol]])),
    threatbark = as.integer(round(.data[[hi_barkcol]]))
  ) %>%
  attach_auc(hi_pubcol) %>%
  filter(!is.na(auc_z))

## Hostile behaviour, per credible condition (Alone, Profile), pooled sexes
for (cond in c("Alone", "Profile")) {
  cat("\n--- Acute stress: Hostile behaviour (", cond, ") ---\n")
  build <- function(mode) {
    z <- if (mode == "main") "auc_z" else "auc_z_no_out"
    d <- hi %>% filter(condition == cond, !is.na(.data[[z]])) %>%
      transmute(hostile, auc_z = .data[[z]])
    if (nrow(d) < 4) return(NULL)
    d
  }
  results[[paste0("B_host_", cond)]] <-
    run_assoc(build, bf(hostile ~ auc_z), negbinomial(), priors_negbin,
              "b_auc_z", "B. Viremia-outcome",
              paste0("Hostile (", cond, ")"), "ZIKV, both sexes",
              "log-rate slope (per SD log-AUC)")
}

## Threat barks, within females (condition-adjusted, animal random intercept)
cat("\n--- Acute stress: Threat barks (females) ---\n")
build_bark <- function(mode) {
  z <- if (mode == "main") "auc_z" else "auc_z_no_out"
  d <- hi %>% filter(sex == "Female", !is.na(.data[[z]])) %>%
    transmute(threatbark, auc_z = .data[[z]], condition, animal)
  if (nrow(d) < 4) return(NULL)
  d
}
results[["B_bark_F"]] <-
  run_assoc(build_bark, bf(threatbark ~ auc_z + condition + (1 | animal)),
            negbinomial(), priors_negbin, "b_auc_z", "B. Viremia-outcome",
            "Threat barks", "Female ZIKV, condition-adj",
            "log-rate slope (per SD log-AUC)")

###############################################################################
## ANALYSIS B3 -- Attachment (negative binomial, log-rate)
##   Male scream vocalizations
###############################################################################
cat("\n\n################ ANALYSIS B3: attachment association ###############\n")

att <- read_csv(ATTACH_CSV, show_col_types = FALSE) %>% clean_names()
att_sexcol <- names(att)[str_detect(names(att), "^sex")][1]
att_pubcol <- names(att)[str_detect(names(att), "pub")][1]
att_scrcol <- names(att)[str_detect(names(att), "scream")][1]

att <- att %>%
  mutate(sex    = factor(ifelse(.data[[att_sexcol]] == 1, "Male", "Female"),
                         levels = c("Female", "Male")),
         scream = as.integer(round(.data[[att_scrcol]]))) %>%
  attach_auc(att_pubcol) %>%
  filter(!is.na(auc_z))

cat("\n--- Attachment: Scream vocalisations (males) ---\n")
build_scream <- function(mode) {
  z <- if (mode == "main") "auc_z" else "auc_z_no_out"
  d <- att %>% filter(sex == "Male", !is.na(scream), !is.na(.data[[z]])) %>%
    transmute(scream, auc_z = .data[[z]])
  if (nrow(d) < 4) return(NULL)
  d
}
results[["B_scream_M"]] <-
  run_assoc(build_scream, bf(scream ~ auc_z), negbinomial(), priors_negbin,
            "b_auc_z", "B. Viremia-outcome", "Scream (attachment)",
            "Male ZIKV", "log-rate slope (per SD log-AUC)")

###############################################################################
## ANALYSIS B4 -- INAS domains (ordinal cumulative-logit over age)
##   State control (= Temperament), Neuromotor, Orientation
##   DV_ord ~ auc_z * Age_c + (1 + Age_c | animal)   within ZIKV
##   (random age SLOPE per animal so the auc_z:Age_c interaction is not
##    overconfident from within-animal pseudo-replication across weeks/items)
##   b_auc_z            = latent-scale shift per SD log-AUC (avg age)
##   b_auc_z:Age_c      = does burden modulate the developmental trajectory
###############################################################################
cat("\n\n################ ANALYSIS B4: INAS associations ####################\n")

inas <- read_csv(INAS_CSV, show_col_types = FALSE)
inas <- inas %>%
  mutate(
    pub     = as.integer(.data[["Pub_#"]]),
    animal  = factor(Animal_ID),
    Age_c   = Age_Weeks - mean(unique(Age_Weeks))
  ) %>%
  left_join(auc %>% select(pub, auc_z, auc_z_no_out), by = "pub") %>%
  filter(!is.na(auc_z))                       # ZIKV animals only

inas_targets <- tribble(
  ~col,                    ~label,
  "Temperament_Scores_02", "State control (INAS)",
  "Neuromotor_Scores_02",  "Neuromotor (INAS)",
  "Orientation_Scores_02", "Orientation (INAS)"
)

for (i in seq_len(nrow(inas_targets))) {
  tg <- inas_targets[i, ]
  cat("\n--- INAS:", tg$label, "---\n")
  build <- function(mode) {
    z <- if (mode == "main") "auc_z" else "auc_z_no_out"
    d <- inas %>%
      filter(!is.na(.data[[tg$col]]), !is.na(.data[[z]])) %>%
      transmute(DV_ord = ordered(.data[[tg$col]]),
                auc_z = .data[[z]], Age_c, animal)
    if (nrow(d) < 20 || nlevels(d$DV_ord) < 2) return(NULL)
    d
  }
  ## main AUC effect
  results[[paste0("B_inas_", tg$col)]] <-
    run_assoc(build, bf(DV_ord ~ auc_z * Age_c + (1 + Age_c | animal)),
              cumulative("logit"), priors_ordinal, "b_auc_z",
              "B. Viremia-outcome", tg$label, "ZIKV, age-adj",
              "latent-logit slope (per SD log-AUC)")
  ## interaction with age (burden x trajectory)
  results[[paste0("B_inasIx_", tg$col)]] <-
    run_assoc(build, bf(DV_ord ~ auc_z * Age_c + (1 + Age_c | animal)),
              cumulative("logit"), priors_ordinal, "b_auc_z:Age_c",
              "B. Viremia-outcome", paste0(tg$label, " x age"), "ZIKV",
              "latent-logit slope change per week")
}

###############################################################################
## COLLATE + WRITE  ->  Supplementary Table S13
###############################################################################
S13 <- bind_rows(results) %>%
  rename(`89% HDI lower` = `HDI 2.5`, `89% HDI upper` = `HDI 97.5`)

cat("\n\n===================== RESULTS (Supplementary Table S13) ===========\n")
print(as.data.frame(S13), row.names = FALSE)

write_csv(S13 %>% filter(str_detect(Analysis, "^A")),
          "viral_load_AUC_sex_burden.csv")
write_csv(S13 %>% filter(str_detect(Analysis, "^B")),
          "viral_load_AUC_outcome_associations.csv")
write_csv(S13, "supplementary_table_S13_AUC_analysis.csv")

cat("\nWrote:\n  viral_load_AUC_sex_burden.csv\n",
    " viral_load_AUC_outcome_associations.csv\n",
    " supplementary_table_S13_AUC_analysis.csv\n")
###############################################################################
## END
###############################################################################
