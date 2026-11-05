###############################################################################
#  Bayesian analysis script for ZIKV MRI data
#
# Requirements: brms, tidyverse, janitor, tidybayes, posterior, bayesplot
###############################################################################

## --- 0. PACKAGES & SETUP ---
if (!require(pacman)) install.packages("pacman")
pacman::p_load(brms, dplyr, janitor, ggplot2, tidyr, tidybayes, posterior, bayesplot, readr, stringr, emmeans)

# ==============================================
# CONFIGURATION SWITCHES FOR VOLUME ADJUSTMENT
# ==============================================
# Options: "none", "ICV", "TBV"
VOLUME_ADJUSTMENT <- "ICV"  # Change to "none", "ICV" or "TBV" as needed

# Validate the choice
if(!VOLUME_ADJUSTMENT %in% c("none", "ICV", "TBV")) {
  stop("VOLUME_ADJUSTMENT must be one of: 'none', 'ICV', 'TBV'")
}

set.seed(123)
options(mc.cores = parallel::detectCores())

cat("--- Packages loaded ---\n")
cat("--- Volume Adjustment:", VOLUME_ADJUSTMENT, "---\n")

## --- 1. LOAD DATA ---
zikv_path <- "ZIKV_3moMRI.csv"
prior_path <- "prior.csv"

if(!file.exists(zikv_path)) stop("ZIKV_3moMRI.csv not found in working directory.")
if(!file.exists(prior_path)) stop("prior.csv not found in working directory.")

zikv_data_raw <- read_csv(zikv_path, show_col_types = FALSE)
prior_data_raw <- read_csv(prior_path, skip = 1, show_col_types = FALSE) %>% clean_names()

zikv_data <- zikv_data_raw %>% clean_names()
cat("Loaded data: ", nrow(zikv_data), " rows experimental,", nrow(prior_data_raw), " rows prior\n")

## --- 2. STANDARDIZE CATEGORICAL NAMES & CREATE GROUP VARIABLES (ROBUST) ---
fcol <- function(df, patterns) {
  cols <- names(df)
  for(p in patterns) {
    m <- cols[str_detect(cols, regex(p, ignore_case = TRUE))]
    if(length(m)>0) return(m[1])
  }
  return(NA_character_)
}

# Use the UN-CLEANED names to find the columns first
sex_col_raw <- fcol(zikv_data_raw, c("^sex$", "sex \\(0=f, 1=m\\)"))
treat_col_raw <- fcol(zikv_data, c("treatment_group", "treatment_group_0", "tr", "trtmnt", "tr.*group", "tr...al"))

# Convert raw column names to clean names for use with the cleaned dataframe
sex_col <- make_clean_names(sex_col_raw)
treat_col <- make_clean_names(treat_col_raw)

cat("Detected columns: sex ->", sex_col, "treatment ->", treat_col, "\n")

if(is.na(sex_col) || is.na(treat_col)) {
  stop("Could not automatically detect sex or treatment_group columns.")
}

# Rename using the dynamically found clean names
zikv_data <- zikv_data %>%
  rename(
    sex_raw = all_of(sex_col),
    treatment_group_raw = all_of(treat_col)
  )

# Standardize sex column robustly
zikv_data <- zikv_data %>%
  mutate(
    sex_standardized = case_when(
      tolower(as.character(sex_raw)) %in% c("f", "female", "0") ~ "Female",
      tolower(as.character(sex_raw)) %in% c("m", "male", "1") ~ "Male",
      TRUE ~ NA_character_
    ),
    sex = factor(sex_standardized, levels = c("Female", "Male"))
  )

if(any(is.na(zikv_data$sex) & !is.na(zikv_data$sex_raw))) {
  warning("Some 'sex' values could not be mapped and were set to NA. Check 'sex_raw' column.")
}

zikv_data <- zikv_data %>%
  mutate(
    treatment_group = factor(case_when(
      treatment_group_raw == 0 ~ "UIC",
      treatment_group_raw == 1 ~ "PIC",
      treatment_group_raw == 2 ~ "Zika"
    ), levels = c("UIC", "PIC", "Zika")),
    group = factor(ifelse(treatment_group %in% c("UIC","PIC"), "Control", "Zika"), levels = c("Control","Zika"))
  )

# T-tests to compare ICV and TBV 
t.test(total_brain_vol_mm3 ~ treatment_group, 
       data = droplevels(subset(zikv_data, treatment_group %in% c("UIC", "PIC"))))

t.test(total_intracranial_vol_mm3 ~ treatment_group, 
       data = droplevels(subset(zikv_data, treatment_group %in% c("UIC", "PIC"))))

## NOTE: The t-tests above check UIC vs PIC only on the global measures (TBV, ICV). 
## To confirm if the ZIKV findings hold when the two control groups are analysed SEPARATELY and
## report UIC-vs-PIC *regional* volume differences, we perform the SENSITIVITY ANALYSIS in Section 10 at the end of this file.

## --- 3. DEFINE REGIONS LIST ---
regions_to_analyze <- c("total_hippocampus_mm3", "total_amygdala_mm3", "total_csf_volume_mm3", "lateral_ventricle_mm3",
                        "rt_temporal_auditory_gm_wm", "rt_temporal_visual_gm_wm", "rt_temporal_limbic_gm_wm", "rt_occipital_gm_wm",
                        "lt_temporal_auditory_gm_wm", "lt_temporal_visual_gm_wm", "lt_temporal_limbic_gm_wm", "lt_occipital_gm_wm")     

# Check which regions are available in the cleaned data
available_regions <- regions_to_analyze[regions_to_analyze %in% names(zikv_data)]
missing_regs <- regions_to_analyze[!regions_to_analyze %in% names(zikv_data)]

if(length(missing_regs) > 0) {
  warning("The following regions were not found and will be skipped: ", paste(missing_regs, collapse = ", "))
}


## --- 4. PROCESS PRIOR CSV ---
prior_cols <- names(prior_data_raw)

# Logic to find the Sex column in the Prior CSV
prior_sex_col <- fcol(prior_data_raw, c("sex", "sex.*female.*male", "sex_female_0_male_1"))
if(is.na(prior_sex_col)) {
  # Fallback: try to find a column containing "sex"
  prior_sex_col <- prior_cols[str_detect(prior_cols, "sex")][1]
}
cat("Prior Sex Column detected:", prior_sex_col, "\n")

find_prior_col <- function(prior_df_cols, target_reg) {
  # 1. Normalize the target region name for matching
  clean_key <- target_reg %>%
    str_remove("^total_") %>%
    str_remove("_mm3$") %>%
    str_remove("_gm_wm$") %>%
    str_remove("_gm$") %>%
    str_remove("_wm$")
  
  # 2. Handle specific mappings (ZIKV name part -> Prior name part)
  if(str_detect(target_reg, "intracranial")) clean_key <- "icv"
  if(str_detect(target_reg, "brain_vol")) clean_key <- "tbv"
  
  # 3. First pass: look for exact containment of the clean key
  candidate_cols <- prior_df_cols[str_detect(prior_df_cols, regex(clean_key, ignore_case = TRUE))]
  
  # Filter out 'sem' or secondary columns if possible (though exact match usually prefers the base name)
  # In the new file, duplicates have suffixes like _2, _3 or .1. We usually want the first one.
  if (length(candidate_cols) > 0) return(candidate_cols[1])
  
  # 4. Second pass: Try known keywords
  table_try <- c("hippocampus", "amygdala", "lateral_ventricle", "csf",
                 "temporal_auditory", "temporal_visual", "occipital_gm", "occipital_wm", 
                 "icv", "tbv")
  
  for (p in table_try) {
    if (str_detect(clean_key, p)) {
      cand2 <- prior_df_cols[str_detect(prior_df_cols, regex(p, ignore_case = TRUE))]
      if (length(cand2) > 0) return(cand2[1])
    }
  }
  return(NA_character_)
}


priors_list <- list()
cat("\n--- Mapping and Calculating SEX-SPECIFIC Priors ---\n")
for(reg in available_regions) {
  prior_col <- find_prior_col(prior_cols, reg)
  priors_list[[reg]] <- list(prior_col = NA, mu_f = NA, sd_f = NA, n_f = 0, mu_m = NA, sd_m = NA, n_m = 0)
  
  if(is.na(prior_col)) {
    cat("Prior not found for:", reg, "\n")
    next
  }
  
  priors_list[[reg]]$prior_col <- prior_col
  cat("Processing Region:", reg, "<-", prior_col, "\n")
  
  
  # Extract Females (0)
  vec_f <- prior_data_raw %>% 
    filter(!!sym(prior_sex_col) == 0) %>% 
    pull(!!prior_col) %>% 
    as.numeric()
  
  vec_f <- vec_f[!is.na(vec_f)]
  
  if(length(vec_f) >= 3) {
    priors_list[[reg]]$mu_f <- mean(vec_f)
    priors_list[[reg]]$sd_f <- sd(vec_f)
    priors_list[[reg]]$n_f <- length(vec_f)
    cat("  Female prior: n=", length(vec_f), ", mu=", round(mean(vec_f),2), ", sd=", round(sd(vec_f),2), "\n")
  } else {
    cat("  Female prior: Not enough data points (n=", length(vec_f), ")\n")
  }
  
  # Extract Males (1)
  vec_m <- prior_data_raw %>% 
    filter(!!sym(prior_sex_col) == 1) %>% 
    pull(!!prior_col) %>% 
    as.numeric()
  
  vec_m <- vec_m[!is.na(vec_m)]
  
  if(length(vec_m) >= 3) {
    priors_list[[reg]]$mu_m <- mean(vec_m)
    priors_list[[reg]]$sd_m <- sd(vec_m)
    priors_list[[reg]]$n_m <- length(vec_m)
    cat("  Male prior:   n=", length(vec_m), ", mu=", round(mean(vec_m),2), ", sd=", round(sd(vec_m),2), "\n")
  } else {
    cat("  Male prior:   Not enough data points (n=", length(vec_m), ")\n")
  }
}


## --- 5. CREATE MODELING DATASET ---
base_vars <- c("group", "sex", available_regions)

if(VOLUME_ADJUSTMENT == "ICV") {
  modeling_data <- zikv_data %>%
    select(all_of(base_vars), total_intracranial_vol_mm3) %>%
    filter(!is.na(sex), !is.na(group), !is.na(total_intracranial_vol_mm3)) %>%
    rename(covariate = total_intracranial_vol_mm3)
  cat("Created modeling dataset WITH ICV adjustment (n =", nrow(modeling_data), ")\n")
  
} else if(VOLUME_ADJUSTMENT == "TBV") {
  modeling_data <- zikv_data %>%
    select(all_of(base_vars), total_brain_vol_mm3) %>%
    filter(!is.na(sex), !is.na(group), !is.na(total_brain_vol_mm3)) %>%
    rename(covariate = total_brain_vol_mm3)
  cat("Created modeling dataset WITH TBV adjustment (n =", nrow(modeling_data), ")\n")
  
} else {  # "none"
  modeling_data <- zikv_data %>%
    select(all_of(base_vars)) %>%
    filter(!is.na(sex), !is.na(group))
  cat("Created modeling dataset WITHOUT volume adjustment (n =", nrow(modeling_data), ")\n")
}





## --- 6. FIT PER-REGION MODELS (RAW SCALE PRIORS) ---
per_region_models <- list()
all_posthoc_draws <- list()

for (reg in available_regions) {
  cat("\n--- Processing region:", reg, " ---\n")
  
  # Create region-specific dataset
  if(VOLUME_ADJUSTMENT != "none") {
    df <- modeling_data %>%
      select(group, sex, covariate, region_value = all_of(reg)) %>%
      filter(!is.na(region_value)) %>%
      droplevels()
  } else {
    df <- modeling_data %>%
      select(group, sex, region_value = all_of(reg)) %>%
      filter(!is.na(region_value)) %>%
      droplevels()
  }
  
  if (nrow(df) < 6) {
    warning("Skipping ", reg, ": fewer than 6 non-NA rows.")
    next
  }
  
  # Make sure factor levels are clean
  df$group <- factor(df$group, levels = c("Control", "Zika"))
  df$sex <- factor(df$sex, levels = c("Female", "Male"))
  
  # NO STANDARDIZATION - use raw values
  df <- df %>% mutate(outcome_raw = region_value)
  
  # Standardize covariate if needed (but not outcome)
  if(VOLUME_ADJUSTMENT != "none") {
    cov_mean <- mean(df$covariate, na.rm = TRUE)
    cov_sd <- sd(df$covariate, na.rm = TRUE)
    df <- df %>% mutate(covariate_z = (covariate - cov_mean) / cov_sd)
  }
  
  # Model formula - using RAW outcome
  if(VOLUME_ADJUSTMENT != "none") {
    form <- bf(outcome_raw ~ 0 + group:sex + covariate_z)
    cat("Using formula with", VOLUME_ADJUSTMENT, "adjustment:", deparse(form$formula), "\n")
  } else {
    form <- bf(outcome_raw ~ 0 + group:sex)
    cat("Using formula without volume adjustment:", deparse(form$formula), "\n")
  }
  
  # Check available parameter names
  cat("Checking parameter names for", reg, "...\n")
  prior_table <- get_prior(formula = form, data = df, family = gaussian())
  print(prior_table)
  
  # Get coefficient parameters and clean them
  coef_params <- prior_table$coef[prior_table$class == "b" & prior_table$coef != ""]
  coef_params <- trimws(coef_params)
  cat("Available coefficient parameters:", paste(coef_params, collapse = ", "), "\n")
  
  # Build priors list with RAW SCALE values
  priors_as_list <- list()
  
  # Process each coefficient parameter
  for(param in coef_params) {
    cat("Processing parameter:", param, "\n")
    
    # Skip covariate_z in this loop - we'll add it separately
    if(param == "covariate_z") {
      cat("  -> Skipping covariate_z (will add separately)\n")
      next
    }
    
    if(param == "groupControl:sexFemale") {
      if(!is.na(priors_list[[reg]]$mu_f)) {
        # Use RAW SCALE priors directly (no standardization)
        prior_spec <- sprintf("normal(%f, %f)", 
                              priors_list[[reg]]$mu_f, 
                              max(priors_list[[reg]]$sd_f, 10))  # minimum SD of 10 mm³
        cat("  -> Setting Control:Female prior (raw scale):", prior_spec, "\n")
        priors_as_list[[length(priors_as_list) + 1]] <- set_prior(prior_spec, class = "b", coef = param)
      } else {
        # Generic brain region prior in mm³
        priors_as_list[[length(priors_as_list) + 1]] <- set_prior("normal(2000, 500)", class = "b", coef = param)
      }
    } else if(param == "groupControl:sexMale") {
      if(!is.na(priors_list[[reg]]$mu_m)) {
        # Use RAW SCALE priors directly (no standardization)
        prior_spec <- sprintf("normal(%f, %f)", 
                              priors_list[[reg]]$mu_m, 
                              max(priors_list[[reg]]$sd_m, 10))  # minimum SD of 10 mm³
        cat("  -> Setting Control:Male prior (raw scale):", prior_spec, "\n")
        priors_as_list[[length(priors_as_list) + 1]] <- set_prior(prior_spec, class = "b", coef = param)
      } else {
        # Generic brain region prior in mm³
        priors_as_list[[length(priors_as_list) + 1]] <- set_prior("normal(2000, 500)", class = "b", coef = param)
      }
    } else {
      # Default prior for other parameters (Zika groups) - weakly informative
      cat("  -> Setting default prior for:", param, "\n")
      priors_as_list[[length(priors_as_list) + 1]] <- set_prior("normal(4000, 20000)", class = "b", coef = param)
    }
  }
  
  # Volume covariate prior (if needed)
  if(VOLUME_ADJUSTMENT != "none") {
    cat("Adding covariate_z prior for", VOLUME_ADJUSTMENT, "adjustment\n")
    priors_as_list[[length(priors_as_list) + 1]] <- set_prior("normal(0, 3000)", class = "b", coef = "covariate_z")
  }
  
  # Residual prior - needs to be much larger for raw scale data
  priors_as_list[[length(priors_as_list) + 1]] <- set_prior("student_t(3, 0, 500)", class = "sigma")
  
  # Combine priors
  priors_model <- do.call(c, priors_as_list)
  
  # Print final prior specification for debugging
  cat("Final prior specification:\n")
  print(priors_model)
  
  cat("Fitting brms model for", reg, "... (this may take a while)\n")
  fit <- brm(
    formula = form,
    data = df,
    prior = priors_model,
    family = gaussian(),
    iter = 4000, warmup = 1000, chains = 4,
    control = list(adapt_delta = 0.95, max_treedepth = 15),
    seed = 123,
    file = paste0("brms_per_region_", reg)
  )
  
  per_region_models[[reg]] <- fit
  
  group_eff_draws <- emmeans(fit, specs = pairwise ~ group | sex) %>% gather_emmeans_draws()
  sex_eff_draws   <- emmeans(fit, specs = pairwise ~ sex | group) %>% gather_emmeans_draws()
  
  all_posthoc_draws[[reg]] <- bind_rows(group_eff_draws, sex_eff_draws) %>% mutate(region = reg)
}

cat("\n=== ANALYSIS COMPLETE ===\n")
cat("Models fitted using RAW SCALE priors (no standardization)\n")
cat("Results are in mm³ units, not standardized units\n")




## --- 7. GENERATE CUSTOM PLOT (WITH STANDARDIZATION SWITCH) ---
cat("\n--- Preparing draws for plotting and summary ---\n")

# CONFIGURATION SWITCH FOR PLOTTING
# Set to TRUE for a standardized plot (Cohen's d) or FALSE for a raw plot (mm³)
PLOT_STANDARDIZED <- TRUE

# Combine all posthoc draws from the models
final_posthoc_draws <- bind_rows(all_posthoc_draws)

# Calculate the pooled standard deviation for each region to use as the standardizer
cat("Calculating pooled standard deviations for standardization...\n")
pooled_sds <- zikv_data %>%
  pivot_longer(
    cols = all_of(available_regions),
    names_to = "region",
    values_to = "volume"
  ) %>%
  filter(!is.na(volume)) %>%
  group_by(region) %>%
  summarise(pooled_sd = sd(volume, na.rm = TRUE))

# Add standardized values to the draws dataframe
processed_draws <- final_posthoc_draws %>%
  left_join(pooled_sds, by = "region") %>%
  mutate(
    standardized_value = -.value / pooled_sd # -ve for Zika- Control
  )

# --- Plotting Logic ---
cat("\n--- Generating custom plot for Zika vs. Control contrasts ---\n")

new_region_labels <- c(
  "total_csf_volume_mm3" = "Total CSF", 
  "lateral_ventricle_mm3" = "Lateral Ventricle",
  "total_amygdala_mm3" = "Amygdala", 
  "total_hippocampus_mm3" = "Hippocampus", 
  "lt_temporal_limbic_gm_wm" = "L. Temporal-Limbic", 
  "rt_temporal_limbic_gm_wm" = "R. Temporal-Limbic",
  "lt_temporal_auditory_gm_wm" = "L. Temporal-Auditory",
  "rt_temporal_auditory_gm_wm" = "R. Temporal-Auditory", 
  "lt_temporal_visual_gm_wm" = "L. Temporal-Visual", 
  "rt_temporal_visual_gm_wm" = "R. Temporal-Visual", 
  "lt_occipital_gm_wm" = "L. Occipital",
  "rt_occipital_gm_wm" = "R. Occipital")

# Define which labels should be bold based on their display name
bold_labels <- ifelse(new_region_labels %in% c("Amygdala", "Lat. Ventricle", "Total CSF", 
                                               "LT Auditory GM WM", "RT Visual GM WM", 
                                               "LT Visual GM WM", "RT Limbic GM WM", 
                                               "LT Limbic GM WM", "RT Occipital GM WM", 
                                               "LT Occipital GM WM"), 
                      "bold", "plain")

# Prepare data for plotting based on the switch
plot_data <- processed_draws %>%
  filter(contrast == "Control - Zika") %>%
  mutate(
    region = dplyr::recode(region, !!!new_region_labels),
    region = factor(region, levels = unique(as.character(new_region_labels))),
    sex = factor(sex, levels = c("Female", "Male"))
  )

# Dynamically set plot aesthetics and labels based on the switch
if (PLOT_STANDARDIZED) {
  y_aes <- aes(y = standardized_value, x = region, fill = sex)
  y_lab <- "Standardized Difference (Cohen's d)"
  plot_title <- "Neuroimaging Assessment: ZIKV - Control Contrasts (ICV-adjusted)"
  file_name <- "zika_effect_by_sex_plot_standardized.png"
  cat("Plotting on STANDARDIZED scale.\n")
} else {
  y_aes <- aes(y = -.value, x = region, fill = sex)
  y_lab <- "Estimated Difference (mm³)"
  plot_title <- "Zika - Control Contrast Across Brain Regions (Raw Scale)"
  file_name <- "zika_effect_by_sex_plot_raw.png"
  cat("Plotting on RAW (natural) scale.\n")
}

zika_effect_plot <- ggplot(plot_data, y_aes) +
  stat_halfeye(
    position = position_dodge(width = 0.75),
    alpha = 0.7,
    .width = c(0.66, 0.89),
    interval_color = "black",
    slab_color = NA
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  scale_fill_manual(values = c("Female" = "#fb8072", "Male" = "#80b1d3"), name = "Sex") +
  labs(
    title = plot_title,
    subtitle = "Positive values indicate Zika > Control",
    y = y_lab,
    x = "Brain Region"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1),# face = bold_labels),
    legend.title = element_text(size = 16),
    legend.text = element_text(size = 14),
    axis.title = element_text(size = 18),
    axis.text = element_text(size = 16),
    plot.title = element_text(size = 20, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40")
  ) +
  coord_cartesian(ylim = c(-3, 5))


print(zika_effect_plot)
ggsave(file_name, zika_effect_plot, width = 18, height = 10, dpi = 300)
cat("Custom Zika effect plot saved to:", file_name, "\n")


## --- 8. INTERACTION EFFECT ANALYSIS (RAW AND STANDARDIZED) ---
interaction_tests <- list()

for (reg in names(per_region_models)) {
  cat("\n--- Testing interaction for region:", reg, " ---\n")
  
  fit <- per_region_models[[reg]]
  
  # Get the standard deviation for the current region to standardize the interaction
  region_sd <- pooled_sds %>% filter(region == reg) %>% pull(pooled_sd)
  
  # Calculate interaction contrast on the RAW scale
  interaction_draws <- fit %>%
    spread_draws(`b_groupZika:sexMale`, `b_groupControl:sexMale`, 
                 `b_groupZika:sexFemale`, `b_groupControl:sexFemale`) %>%
    mutate(
      zika_effect_male = `b_groupZika:sexMale` - `b_groupControl:sexMale`,
      zika_effect_female = `b_groupZika:sexFemale` - `b_groupControl:sexFemale`,
      interaction_effect = zika_effect_male - zika_effect_female,
      # Also create the STANDARDIZED interaction effect
      interaction_effect_std = interaction_effect / region_sd
    )
  
  # Summarize interaction effect on BOTH scales
  interaction_summary <- interaction_draws %>%
    summarise(
      # Raw scale summary
      estimate_raw = mean(interaction_effect),
      lower_89_raw = quantile(interaction_effect, 0.055),
      upper_89_raw = quantile(interaction_effect, 0.945),
      
      # Standardized scale summary
      estimate_std = mean(interaction_effect_std),
      lower_89_std = quantile(interaction_effect_std, 0.055),
      upper_89_std = quantile(interaction_effect_std, 0.945),
      
      .groups = "drop"
    ) %>%
    mutate(region = reg)
  
  interaction_tests[[reg]] <- interaction_summary
  
  cat("Interaction Effect (Raw Scale): ", round(interaction_summary$estimate_raw, 2), 
      "mm³, 89% CI: [", round(interaction_summary$lower_89_raw, 2), ", ", 
      round(interaction_summary$upper_89_raw, 2), "]\n")
  
  cat("Interaction Effect (Std. Scale):", round(interaction_summary$estimate_std, 3), 
      ", 89% CI: [", round(interaction_summary$lower_89_std, 3), ", ", 
      round(interaction_summary$upper_89_std, 3), "]\n")
}

# Combine all interaction tests into a single data frame
interaction_results <- bind_rows(interaction_tests)
cat("\n\n--- Combined Interaction Analysis Results ---\n")
print(interaction_results)
write.csv(interaction_results, "interaction_results.csv")

## --- 9. CREATE CONTROL - ZIKA CONTRAST SUMMARY TABLES (RAW AND STANDARDIZED) ---

# --- Configuration for Meaningful Effects ---
# Define what constitutes a "meaningful" effect for summary tables.
# This makes it easy to change the thresholds later.
MEANINGFUL_THRESHOLD_RAW <- 50    # In absolute units (e.g., 50 mm³)
MEANINGFUL_THRESHOLD_STD <- 0.2   # In standard deviation units (for Cohen's d)


# Start with the processed draws containing both raw and standardized values
contrast_summary <- processed_draws %>%
  filter(contrast == "Control - Zika") %>%
  group_by(region, sex) %>%
  summarise(
    # Raw scale summaries
    estimate_raw = mean(-.value), # -ve for Zika - Control
    lower_89_raw = quantile(-.value, 0.055),
    upper_89_raw = quantile(-.value, 0.945),
    prob_negative_raw = mean(-.value < 0),
    prob_meaningful_raw = mean(abs(-.value) > MEANINGFUL_THRESHOLD_RAW),
    
    # Standardized scale summaries
    estimate_std = mean(standardized_value), # standardized_value includes Zika - Control
    lower_89_std = quantile(standardized_value, 0.055),
    upper_89_std = quantile(standardized_value, 0.945),
    prob_negative_std = mean(standardized_value < 0),
    prob_meaningful_std = mean(abs(standardized_value) > MEANINGFUL_THRESHOLD_STD),
    
    .groups = "drop"
  ) %>%
  mutate(
    # Formatted strings for easier reading
    ci_89_raw_formatted = sprintf("%.2f [%.2f, %.2f]", estimate_raw, lower_89_raw, upper_89_raw),
    ci_89_std_formatted = sprintf("%.3f [%.3f, %.3f]", estimate_std, lower_89_std, upper_89_std),
    
    # Clean region labels for presentation
    region_clean = dplyr::recode(region, !!!new_region_labels)
  ) %>%
  arrange(region_clean, sex)


# --- RAW SCALE SUMMARY TABLE ---
contrast_table_89_raw <- contrast_summary %>%
  select(region_clean, sex, estimate_raw, lower_89_raw, upper_89_raw, prob_negative_raw, prob_meaningful_raw) %>%
  rename(
    Region = region_clean,
    Sex = sex,
    `Estimate (mm³)` = estimate_raw,
    `Lower 89%` = lower_89_raw,
    `Upper 89%` = upper_89_raw,
    `P(Negative)` = prob_negative_raw,
    `P(|Effect| > THRESHOLD)` = prob_meaningful_raw
  ) %>%
  mutate(
    across(c(`Estimate (mm³)`, `Lower 89%`, `Upper 89%`), ~ round(.x, 2)),
    across(c(`P(Negative)`, `P(|Effect| > THRESHOLD)`), ~ round(.x, 3))
  )

cat("\n\n=== CONTROL - ZIKA CONTRAST SUMMARY (RAW SCALE, 89% CI) ===\n")
cat("Meaningful effect threshold: >", MEANINGFUL_THRESHOLD_RAW, "mm³\n")
print(contrast_table_89_raw, n = Inf)
write_csv(contrast_table_89_raw, "control_zika_contrasts_raw_89ci.csv")
cat("Raw scale table saved to: control_zika_contrasts_raw_89ci.csv\n")


# --- STANDARDIZED SCALE SUMMARY TABLE ---
contrast_table_89_std <- contrast_summary %>%
  select(region_clean, sex, estimate_std, lower_89_std, upper_89_std, prob_negative_std, prob_meaningful_std) %>%
  rename(
    Region = region_clean,
    Sex = sex,
    `Estimate (SDs)` = estimate_std,
    `Lower 89%` = lower_89_std,
    `Upper 89%` = upper_89_std,
    `P(Negative)` = prob_negative_std,
    `P(|Effect| > THRESHOLD)` = prob_meaningful_std
  ) %>%
  mutate(
    across(c(`Estimate (SDs)`, `Lower 89%`, `Upper 89%`), ~ round(.x, 3)),
    across(c(`P(Negative)`, `P(|Effect| > THRESHOLD)`), ~ round(.x, 3))
  )

cat("\n\n=== ZIKA - CONTROL CONTRAST SUMMARY (STANDARDIZED SCALE, 89% CI) ===\n")
cat("Meaningful effect threshold: >", MEANINGFUL_THRESHOLD_STD, "SDs\n")
print(contrast_table_89_std, n = Inf)
write_csv(contrast_table_89_std, "control_zika_contrasts_std_89ci.csv")
cat("Standardized scale table saved to: control_zika_contrasts_std_89ci.csv\n")


# --- WIDE FORMAT TABLE (STANDARDIZED) ---
# Wide format is often most useful with standardized effects for easy comparison
contrast_table_wide_std <- contrast_summary %>%
  select(region_clean, sex, ci_89_std_formatted, prob_meaningful_std) %>%
  pivot_wider(
    names_from = sex,
    values_from = c(ci_89_std_formatted, prob_meaningful_std),
    names_sep = "_"
  ) %>%
  rename(
    Region = region_clean,
    `Female: Estimate (SDs) [89% CI]` = ci_89_std_formatted_Female,
    `Male: Estimate (SDs) [89% CI]` = ci_89_std_formatted_Male,
    `Female: P(|Effect| > THRESHOLD)` = prob_meaningful_std_Female,
    `Male: P(|Effect| > THRESHOLD)` = prob_meaningful_std_Male
  )

cat("\n\n=== ZIKA - CONTROL CONTRAST SUMMARY (STANDARDIZED, WIDE FORMAT) ===\n")
cat("Meaningful effect threshold: >", MEANINGFUL_THRESHOLD_STD, "SDs\n")
print(contrast_table_wide_std, n = Inf)
write_csv(contrast_table_wide_std, "control_zika_contrasts_wide_std_89ci.csv") 
cat("Wide format standardized table saved to: control_zika_contrasts_wide_std_89ci.csv\n")


# --- Quick visual check using standardized effects ---
significant_effects <- contrast_summary %>%
  filter(prob_meaningful_std > 0.8) %>%  # High probability of meaningful standardized effect
  select(region_clean, sex, estimate_std, lower_89_std, upper_89_std, prob_meaningful_std) %>%
  arrange(desc(prob_meaningful_std))

if(nrow(significant_effects) > 0) {
  cat("\n\n=== REGIONS WITH LIKELY MEANINGFUL STANDARDIZED EFFECTS (P > 0.8) ===\n")
  print(significant_effects, n = Inf)
} else {
  cat("\n\nNo regions showed high probability (>0.8) of meaningful standardized effects.\n")
}

write.csv(significant_effects, "significant_regions.csv")



###############################################################################
## --- 10. SENSITIVITY ANALYSIS: CONTROL GROUPS ANALYSED SEPARATELY ----------
##
##  The primary analysis (Sections 5-9) pools UIC + PIC into a single "Control"
##  group. Here we re-fit every regional model keeping all THREE treatment
##  groups separate (UIC, PIC, Zika) so we can:
##    (a) confirm if the Zika effects are consistent whether compared against UIC
##        or against PIC (i.e. the pooled contrast is not an artefact of
##        pooling), and
##    (b) directly report whether the two control groups differ in regional
##        brain volume (PIC vs UIC).
##
##  Every contrast table now also carries the posterior MEDIAN, the 89% HDI, 
##  the probability of direction (pd) and P(|d| > 0.2), 
##  in addition to the original mean / 89% ETI / prob_gt0. 
##  A new subsection 10f reuses the 2-group models (per_region_models) to report
##  the group x sex INTERACTION the same holistic way.
##
##  MODEL RE-USE: the 3-group models are cached via brms `file=`, so re-running
##  loads them instead of refitting (numbers stay identical to the saved run).
##  If you load your saved manuscript workspace first, 10b will reuse any
##  per_region_models_3 already in memory, and 10f will reuse per_region_models.
##
##  Priors: the historical CONTROL priors (priors_list, computed in Section 4)
##  are applied to BOTH the UIC and PIC subgroups, because the historical data
##  are uninfected controls and do not distinguish a PIC-type control. The Zika
##  subgroups keep the same weakly-informative prior as the primary analysis.
##
##  CAUTION: split by sex, each control subgroup is very small (n ~ 3). These
##  contrasts are a consistency / sensitivity check; their intervals are
##  expected to be wide and should be interpreted as such.
##
##  Reuses objects created earlier: available_regions, VOLUME_ADJUSTMENT,
##  zikv_data, priors_list, pooled_sds, new_region_labels, per_region_models.
###############################################################################
cat("\n\n=== SENSITIVITY ANALYSIS: UIC, PIC, Zika kept separate ===\n")

## --- 10.0 Holistic-summary helpers (added for revision: pd, HDI, P(meaningful)) ---
# 89% highest-density interval (to match the median/HDI convention in the text)
hdi_bounds <- function(x, ci = 0.89) {
  x <- sort(x); n <- length(x); k <- floor(ci * n)
  if (k < 1L || k >= n) return(c(lower = NA_real_, upper = NA_real_))
  w <- x[(k + 1):n] - x[1:(n - k)]; i <- which.min(w)
  c(lower = x[i], upper = x[i + k])
}
pd_value <- function(x) max(mean(x > 0), mean(x < 0))   # probability of direction
# practical-significance threshold (matches Section 9); define if not already set
if (!exists("MEANINGFUL_THRESHOLD_STD")) MEANINGFUL_THRESHOLD_STD <- 0.2

## --- 10a. Build 3-group modeling dataset (mirrors Section 5) ---
base_vars_3 <- c("treatment_group", "sex", available_regions)
if (VOLUME_ADJUSTMENT == "ICV") {
  modeling_data_3 <- zikv_data %>%
    select(all_of(base_vars_3), total_intracranial_vol_mm3) %>%
    filter(!is.na(sex), !is.na(treatment_group), !is.na(total_intracranial_vol_mm3)) %>%
    rename(covariate = total_intracranial_vol_mm3)
} else if (VOLUME_ADJUSTMENT == "TBV") {
  modeling_data_3 <- zikv_data %>%
    select(all_of(base_vars_3), total_brain_vol_mm3) %>%
    filter(!is.na(sex), !is.na(treatment_group), !is.na(total_brain_vol_mm3)) %>%
    rename(covariate = total_brain_vol_mm3)
} else {  # "none"
  modeling_data_3 <- zikv_data %>%
    select(all_of(base_vars_3)) %>%
    filter(!is.na(sex), !is.na(treatment_group))
}
cat("3-group modeling dataset created (n =", nrow(modeling_data_3), ")\n")
cat("Per-subgroup counts (rows before per-region NA filtering):\n")
print(table(modeling_data_3$treatment_group, modeling_data_3$sex))

## --- 10b. Fit per-region models with treatment_group (3 levels) ---
# Guard the init so a per_region_models_3 loaded from your saved workspace is
# NOT wiped; the loop below reuses any model already present (or cached on disk).
if (!exists("per_region_models_3")) per_region_models_3 <- list()
all_posthoc_draws_3  <- list()

for (reg in available_regions) {
  cat("\n--- [3-group] Processing region:", reg, " ---\n")
  if (VOLUME_ADJUSTMENT != "none") {
    df <- modeling_data_3 %>%
      select(treatment_group, sex, covariate, region_value = all_of(reg)) %>%
      filter(!is.na(region_value)) %>%
      droplevels()
  } else {
    df <- modeling_data_3 %>%
      select(treatment_group, sex, region_value = all_of(reg)) %>%
      filter(!is.na(region_value)) %>%
      droplevels()
  }
  if (nrow(df) < 6) {
    warning("Skipping ", reg, " (3-group): fewer than 6 non-NA rows.")
    next
  }
  df$treatment_group <- factor(df$treatment_group, levels = c("UIC", "PIC", "Zika"))
  df$sex             <- factor(df$sex, levels = c("Female", "Male"))
  df <- df %>% mutate(outcome_raw = region_value)
  if (VOLUME_ADJUSTMENT != "none") {
    cov_mean <- mean(df$covariate, na.rm = TRUE)
    cov_sd   <- sd(df$covariate, na.rm = TRUE)
    df <- df %>% mutate(covariate_z = (covariate - cov_mean) / cov_sd)
    form <- bf(outcome_raw ~ 0 + treatment_group:sex + covariate_z)
  } else {
    form <- bf(outcome_raw ~ 0 + treatment_group:sex)
  }
  
  # Reuse a preloaded model if present; otherwise build priors and fit (brms
  # also caches to disk via file=, so a second run loads rather than refits).
  if (!is.null(per_region_models_3[[reg]])) {
    cat("  -> Reusing preloaded/cached 3-group model for", reg, "(no refit)\n")
    fit <- per_region_models_3[[reg]]
  } else {
    # Build priors for the six group:sex cell means (+ covariate + sigma),
    # matching the primary analysis conventions in Section 6.
    prior_table <- get_prior(formula = form, data = df, family = gaussian())
    coef_params <- trimws(prior_table$coef[prior_table$class == "b" & prior_table$coef != ""])
    priors_as_list <- list()
    for (param in coef_params) {
      if (param == "covariate_z") next  # added separately below
      is_female  <- str_detect(param, "sexFemale")
      is_male    <- str_detect(param, "sexMale")
      is_control <- str_detect(param, "treatment_groupUIC") ||
        str_detect(param, "treatment_groupPIC")
      if (is_control && is_female && !is.na(priors_list[[reg]]$mu_f)) {
        prior_spec <- sprintf("normal(%f, %f)", priors_list[[reg]]$mu_f,
                              max(priors_list[[reg]]$sd_f, 10))
        priors_as_list[[length(priors_as_list) + 1]] <-
          set_prior(prior_spec, class = "b", coef = param)
      } else if (is_control && is_male && !is.na(priors_list[[reg]]$mu_m)) {
        prior_spec <- sprintf("normal(%f, %f)", priors_list[[reg]]$mu_m,
                              max(priors_list[[reg]]$sd_m, 10))
        priors_as_list[[length(priors_as_list) + 1]] <-
          set_prior(prior_spec, class = "b", coef = param)
      } else if (is_control) {
        # control subgroup but no historical prior available -> generic (as Section 6)
        priors_as_list[[length(priors_as_list) + 1]] <-
          set_prior("normal(2000, 500)", class = "b", coef = param)
      } else {
        # Zika subgroups -> weakly informative (identical to primary analysis)
        priors_as_list[[length(priors_as_list) + 1]] <-
          set_prior("normal(4000, 20000)", class = "b", coef = param)
      }
    }
    if (VOLUME_ADJUSTMENT != "none") {
      priors_as_list[[length(priors_as_list) + 1]] <-
        set_prior("normal(0, 3000)", class = "b", coef = "covariate_z")
    }
    priors_as_list[[length(priors_as_list) + 1]] <-
      set_prior("student_t(3, 0, 500)", class = "sigma")
    priors_model <- do.call(c, priors_as_list)
    cat("Fitting 3-group brms model for", reg, "... (this may take a while)\n")
    fit <- brm(
      formula = form, data = df, prior = priors_model, family = gaussian(),
      iter = 4000, warmup = 1000, chains = 4,
      control = list(adapt_delta = 0.95, max_treedepth = 15),
      seed = 123,
      file = paste0("brms_3group_", reg)
    )
    per_region_models_3[[reg]] <- fit
  }
  
  # revpairwise gives "later - earlier" within each sex, i.e.:
  #   "PIC - UIC", "Zika - UIC", "Zika - PIC"
  # so positive = second-named (Zika, or PIC) larger -- same direction as the
  # manuscript's "Zika - Control" convention.
  emm <- emmeans(fit, specs = ~ treatment_group | sex)
  ctr_draws <- contrast(emm, method = "revpairwise") %>% gather_emmeans_draws()
  all_posthoc_draws_3[[reg]] <- ctr_draws %>% mutate(region = reg)
}

## --- 10c. Summarize all pairwise contrasts (raw + standardized) ---
# Original convention (mean + equal-tailed 89% interval) is kept; median, 89%
# HDI, pd and P(|d| > threshold) are ADDED for holistic Bayesian reporting.
posthoc_draws_3 <- bind_rows(all_posthoc_draws_3) %>%
  left_join(pooled_sds, by = "region") %>%
  mutate(standardized_value = .value / pooled_sd)  # already Zika/PIC - control direction
summary_3 <- posthoc_draws_3 %>%
  group_by(region, sex, contrast) %>%
  summarise(
    # posterior mean + equal-tailed 89% interval (unchanged)
    estimate_raw = mean(.value),
    lower_89_raw = quantile(.value, 0.055),
    upper_89_raw = quantile(.value, 0.945),
    estimate_std = mean(standardized_value),
    lower_89_std = quantile(standardized_value, 0.055),
    upper_89_std = quantile(standardized_value, 0.945),
    # ADDED: median + 89% HDI (matches the manuscript's median/HDI convention)
    median_std    = median(standardized_value),
    hdi_lower_std = hdi_bounds(standardized_value)[["lower"]],
    hdi_upper_std = hdi_bounds(standardized_value)[["upper"]],
    # ADDED: holistic direction + practical magnitude
    prob_gt0            = mean(.value > 0),
    pd                  = pd_value(standardized_value),
    prob_meaningful_std = mean(abs(standardized_value) > MEANINGFUL_THRESHOLD_STD),
    .groups = "drop"
  ) %>%
  mutate(
    region_clean = dplyr::recode(region, !!!new_region_labels),
    credible_89  = (lower_89_std > 0) | (upper_89_std < 0)  # 89% CI excludes 0
  ) %>%
  arrange(region_clean, sex, contrast)
cat("\n=== 3-GROUP: ALL PAIRWISE CONTRASTS (within sex, standardized) ===\n")
print(summary_3 %>%
        select(region_clean, sex, contrast, estimate_std, median_std,
               hdi_lower_std, hdi_upper_std, pd, prob_meaningful_std, credible_89),
      n = Inf)
write_csv(summary_3, "sensitivity_3group_all_contrasts.csv")

## --- 10d. Zika-vs-each-control consistency table ---
# If the pooled result is robust, "Zika - UIC" and "Zika - PIC" should agree in
# direction/magnitude with each other and with the pooled "Zika - Control".
zika_vs_controls <- summary_3 %>%
  filter(contrast %in% c("Zika - UIC", "Zika - PIC")) %>%
  select(region_clean, sex, contrast,
         estimate_std, median_std, lower_89_std, upper_89_std,
         hdi_lower_std, hdi_upper_std, prob_gt0, pd, prob_meaningful_std, credible_89) %>%
  arrange(region_clean, sex, contrast)
cat("\n=== ZIKA vs UIC and ZIKA vs PIC (compare with pooled Zika vs Control) ===\n")
print(zika_vs_controls, n = Inf)
write_csv(zika_vs_controls, "sensitivity_zika_vs_each_control.csv")

## --- 10e. UIC vs PIC regional comparison ---
# Pooling is justified for regions where this contrast's 89% CI includes 0.
pic_vs_uic <- summary_3 %>%
  filter(contrast == "PIC - UIC") %>%
  select(region_clean, sex,
         estimate_raw, lower_89_raw, upper_89_raw,
         estimate_std, median_std, lower_89_std, upper_89_std,
         hdi_lower_std, hdi_upper_std, prob_gt0, pd, prob_meaningful_std, credible_89) %>%
  arrange(region_clean, sex)
cat("\n=== PIC vs UIC by region (credible_89 == TRUE flags a control-group difference) ===\n")
print(pic_vs_uic, n = Inf)
write_csv(pic_vs_uic, "control_PIC_vs_UIC_by_region.csv")

## --- 10f. GROUP x SEX INTERACTION with holistic summaries --
# Reuses the 2-group models from Sections 6/8 (per_region_models) -- NO refit.
# Reproduces the Section 8 interaction contrast, i.e. the difference in the
# Zika-vs-Control effect between the sexes:
#     (Zika - Control | male) - (Zika - Control | female)
# standardized by the same pooled SD, and reports it holistically (median, 89%
# HDI, pd, P(|d| > threshold)) so the sex-DIFFERENCE is not collapsed to a
# binary credible/not-credible call.
if (exists("per_region_models") && length(per_region_models) > 0) {
  interaction_full <- lapply(names(per_region_models), function(reg) {
    fit       <- per_region_models[[reg]]
    region_sd <- pooled_sds %>% filter(region == reg) %>% pull(pooled_sd)
    d <- fit %>%
      spread_draws(`b_groupZika:sexMale`,   `b_groupControl:sexMale`,
                   `b_groupZika:sexFemale`, `b_groupControl:sexFemale`) %>%
      mutate(
        zika_male   = `b_groupZika:sexMale`   - `b_groupControl:sexMale`,
        zika_female = `b_groupZika:sexFemale` - `b_groupControl:sexFemale`,
        inter       = zika_male - zika_female,
        inter_std   = inter / region_sd
      )
    h <- hdi_bounds(d$inter_std)
    tibble(
      region              = reg,
      estimate_std        = mean(d$inter_std),
      median_std          = median(d$inter_std),
      lower_89_std        = as.numeric(quantile(d$inter_std, 0.055)),
      upper_89_std        = as.numeric(quantile(d$inter_std, 0.945)),
      hdi_lower_std       = as.numeric(h[["lower"]]),
      hdi_upper_std       = as.numeric(h[["upper"]]),
      prob_gt0            = mean(d$inter_std > 0),
      pd                  = pd_value(d$inter_std),
      prob_meaningful_std = mean(abs(d$inter_std) > MEANINGFUL_THRESHOLD_STD)
    )
  }) %>% bind_rows() %>%
    mutate(region_clean = dplyr::recode(region, !!!new_region_labels)) %>%
    arrange(region_clean)
  
  cat("\n=== GROUP x SEX INTERACTION (Zika effect: male - female), holistic ===\n")
  print(interaction_full %>%
          select(region_clean, estimate_std, median_std,
                 hdi_lower_std, hdi_upper_std, pd, prob_meaningful_std),
        n = Inf)
  write_csv(interaction_full, "interaction_results_full.csv")
  cat("Interaction summary written to: interaction_results_full.csv\n")
} else {
  cat("\nNOTE: per_region_models not in memory; skipping the interaction summary (10f).\n",
      "Run Sections 6-8 first, or load your saved workspace, then re-run Section 10.\n", sep = "")
}

cat("\n=== SENSITIVITY ANALYSIS COMPLETE ===\n")
cat("Outputs written: sensitivity_3group_all_contrasts.csv,",
    "sensitivity_zika_vs_each_control.csv, control_PIC_vs_UIC_by_region.csv,",
    "interaction_results_full.csv\n")