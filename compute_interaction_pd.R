# =====================================================================
# Recompute the neuroimaging group x sex INTERACTION with holistic
# Bayesian summaries (no refitting -- loads the saved per-region models).
#
# For each region it reproduces the interaction contrast exactly as in
# neuroimaging_assessment.R:
#     interaction = (Zika_male - Control_male) - (Zika_female - Control_female)
# standardized by the same pooled SD, and adds:
#     - posterior mean
#     - 89% equal-tailed interval (ETI) and 89% highest-density interval (HDI)
#     - prob_gt0 and pd (probability of direction)
#     - prob_meaningful = P(|d| > 0.2)   (matches MEANINGFUL_THRESHOLD_STD)
#
# Run this from inside the ZIKV-longitudinal-study folder.
# Output: interaction_results_full.csv
# =====================================================================

suppressPackageStartupMessages({
  library(brms); library(dplyr); library(tidyr); library(tidybayes)
  library(posterior); library(janitor); library(stringr)
})

MEANINGFUL_THRESHOLD_STD <- 0.2   # same threshold used in the main script

# --- pooled SDs from the raw data (same standardizer as the main script) ---
dat <- janitor::clean_names(read.csv("ZIKV_3moMRI.csv", check.names = TRUE))

# region keys come from the saved neuroimaging model filenames
rds_files <- list.files(".", pattern = "^brms_3group_.*\\.rds$")
rds_files <- rds_files[!grepl("attach|_va_", rds_files)]
regions   <- sub("\\.rds$", "", sub("^brms_3group_", "", rds_files))

get_sd <- function(r) {
  if (!r %in% names(dat)) stop(paste0("Region column not found in data: '", r,
    "'. Available: ", paste(names(dat), collapse = ", ")))
  sd(dat[[r]], na.rm = TRUE)
}

# 89% HDI helper 
hdi89 <- function(x, ci = 0.89) {
  x <- sort(x); n <- length(x); k <- floor(ci * n)
  w <- x[(k + 1):n] - x[1:(n - k)]; i <- which.min(w)
  c(lower = x[i], upper = x[i + k])
}

results <- lapply(seq_along(regions), function(j) {
  reg <- regions[j]
  cat("Processing:", reg, "\n")
  fit <- readRDS(rds_files[j])
  region_sd <- get_sd(reg)

  d <- fit %>%
    spread_draws(`b_groupZika:sexMale`,   `b_groupControl:sexMale`,
                 `b_groupZika:sexFemale`, `b_groupControl:sexFemale`) %>%
    mutate(
      zika_male   = `b_groupZika:sexMale`   - `b_groupControl:sexMale`,
      zika_female = `b_groupZika:sexFemale` - `b_groupControl:sexFemale`,
      inter       = zika_male - zika_female,
      inter_std   = inter / region_sd
    )

  h <- hdi89(d$inter_std)
  tibble(
    region          = reg,
    estimate_std    = mean(d$inter_std),
    eti_lower_std   = as.numeric(quantile(d$inter_std, 0.055)),
    eti_upper_std   = as.numeric(quantile(d$inter_std, 0.945)),
    hdi_lower_std   = as.numeric(h["lower"]),
    hdi_upper_std   = as.numeric(h["upper"]),
    prob_gt0        = mean(d$inter_std > 0),
    pd              = max(mean(d$inter_std > 0), mean(d$inter_std < 0)),
    prob_meaningful = mean(abs(d$inter_std) > MEANINGFUL_THRESHOLD_STD)
  )
}) %>% bind_rows()

print(as.data.frame(results), digits = 3)
write.csv(results, "interaction_results_full.csv", row.names = FALSE)
cat("\nWrote interaction_results_full.csv\n")
