################################################################################
# Cytokine analysis: innate immune activation in Poly-IC vs uninfected controls
#
# Consistent with the paper's Bayesian approach (brms): each analyte is modeled
# with a Bayesian regression using weakly-informative priors, and effects are
# summarized by the posterior median with an 89% highest-density interval (HDI)
# and the probability of direction (pd).
#
# The ZIKV group has no cytokine measurements (limited blood volume; priority
# given to ZIKV viral quantification), so the informative comparison is
# Poly-IC (PIC) versus uninfected control (UIC). Values are log10 fold change,
# so a Gaussian likelihood on that scale is appropriate; a difference of 0
# corresponds to no change relative to baseline.
#
# Input : ZIKV-R01_LongitudinalStudy_animals&sex_8-14-26.xlsx (sheet All_subjects)
# Output: cytokine_PIC_vs_UIC_contrasts.csv
################################################################################

library(brms)
library(readxl)
library(dplyr)
library(tidyr)
library(tidybayes)

# --- Data ---------------------------------------------------------------------
dat <- read_excel("ZIKV-R01_LongitudinalStudy_animals&sex_8-14-26.xlsx",
                  sheet = "All_subjects")
names(dat) <- c("pub", "cohort", "sex", "group", "vrna_auc",
                "IFNa2a", "IFNg", "IL1Ra", "IP10")

dat <- dat %>%
  filter(group %in% c(0, 1)) %>%                         # 0 = UIC, 1 = PIC (2 = ZIKV has no cytokines)
  mutate(Group = factor(group, levels = c(0, 1), labels = c("UIC", "PIC")))

# analyte column -> display label
cytokines <- c(IFNa2a = "IFN-alpha2a", IFNg = "IFN-gamma",
               IL1Ra = "IL-1Ra",       IP10 = "IP-10")

# --- Weakly-informative priors on the log10 fold-change scale ------------------
# GroupPIC is the PIC - UIC difference (UIC is the reference level).
priors <- c(
  set_prior("normal(0, 2)",       class = "Intercept"),
  set_prior("normal(0, 2)",       class = "b"),
  set_prior("student_t(3, 0, 2.5)", class = "sigma")
)

# --- Fit one model per analyte and extract the PIC - UIC contrast --------------
results <- lapply(names(cytokines), function(col) {
  f   <- as.formula(paste0("`", col, "` ~ Group"))
  fit <- brm(f, data = dat, family = gaussian(), prior = priors,
             iter = 4000, warmup = 1000, chains = 4, cores = 4, seed = 123,
             control = list(adapt_delta = 0.95),
             file = paste0("brms_cytokine_", col))

  draws <- as_draws_df(fit)[["b_GroupPIC"]]          # PIC - UIC
  hdi   <- median_hdci(draws, .width = 0.89)

  tibble(
    cytokine                 = cytokines[[col]],
    estimate_PIC_minus_UIC   = median(draws),
    lower_89_HDI             = hdi$ymin,
    upper_89_HDI             = hdi$ymax,
    prob_PIC_greater_than_UIC = mean(draws > 0),
    pd                       = max(mean(draws > 0), mean(draws < 0)),
    credible_89              = (hdi$ymin > 0) | (hdi$ymax < 0)
  )
}) %>% bind_rows()

print(as.data.frame(results), digits = 3)
write.csv(results, "cytokine_PIC_vs_UIC_contrasts.csv", row.names = FALSE)

cat("\nValues are Poly-IC minus uninfected-control differences in log10 fold change,\n",
    "with the 89% HDI and the posterior probability that Poly-IC exceeds control.\n", sep = "")
cat("--- CYTOKINE ANALYSIS FINISHED ---\n")
