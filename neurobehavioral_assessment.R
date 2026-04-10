################################################################################
# SNIP Bayesian Analysis Script
#
# Description:
# This script is divided into two main parts.
# PART 1 fits the Bayesian ordinal regression models, a process that can be
# time-consuming. It creates model objects in your R environment.
# PART 2 can be run in the same session, or a new session after loading the
# .RData workspace, to perform post-hoc analyses on the fitted models.
#
################################################################################


# --- 0. Setup: Load Libraries and Set Options ---
# Make sure these are installed: install.packages(c("brms", "tidyverse", "readr", "emmeans", "tidybayes"))
library(brms)
library(tidyverse)
library(readr)
library(emmeans)
library(tidybayes)

# Set options for Stan compilation to use multiple cores if available
options(mc.cores = parallel::detectCores())
set.seed(42)

################################################################################
# PART 1: MODEL FITTING
# Run this section to fit the models. This may take time.
################################################################################

# --- 1.1. Data Loading and Preparation ---
# Load the cleaned experimental data
main_data_file <- "SNIP_cleaned.csv"
if (!file.exists(main_data_file)) {
  stop("Error: SNIP_cleaned.csv not found. Please ensure the file is in the working directory.")
}
data_full <- read_csv(main_data_file)

# Prepare factor variables
data_full <- data_full %>%
  mutate(
    Sex_factor = factor(Sex_0f_1m, labels = c("Female", "Male")),
    Treatment_factor = factor(Treatment_Group_0UIC_1PIC_2Zika, labels = c("UIC", "PIC", "Zika")),
    Animal_ID_factor = factor(Animal_ID)
  )

# Load and synthesize prior data from control-like studies
prior_data_file <- "SNIP_prior_data.csv"
if (!file.exists(prior_data_file)) {
  warning("Warning: SNIP_prior_data.csv not found. Proceeding without synthesized prior info.")
  prior_info <- NULL
} else {
  prior_data_raw <- read_csv(prior_data_file, skip = 1)
  colnames(prior_data_raw) <- c(
    "Publication", "sample_size", "females", "males", "Rearing", "age_weeks",
    "Orientation_Mean", "Orientation_SEM", "Motor_Mean", "Motor_SEM",
    "Neuromotor_Mean", "Neuromotor_SEM", "Temperament_Mean", "Temperament_SEM", "X15", "X16"
  )
  prior_data <- prior_data_raw %>%
    select(sample_size, age_weeks, ends_with("_Mean"), ends_with("_SEM")) %>%
    mutate(across(everything(), as.numeric)) %>%
    filter(!is.na(sample_size) & sample_size > 0)
  
  synthesize_prior_info <- function(dv_name) {
    mean_col <- paste0(dv_name, "_Mean"); sem_col <- paste0(dv_name, "_SEM")
    df_prior <- prior_data %>%
      select(n = sample_size, mean = all_of(mean_col), sem = all_of(sem_col)) %>%
      filter(!is.na(mean) & !is.na(sem) & !is.na(n)) %>%
      mutate(sd = sem * sqrt(n))
    if(nrow(df_prior) == 0) return(list(mean = NA, sd = NA))
    weighted_mean <- sum(df_prior$mean * df_prior$n, na.rm = TRUE) / sum(df_prior$n, na.rm = TRUE)
    total_n <- sum(df_prior$n); k_studies <- nrow(df_prior)
    if (total_n <= k_studies) { pooled_sd <- mean(df_prior$sd, na.rm = TRUE) }
    else { pooled_var <- sum(df_prior$sd^2 * (df_prior$n - 1), na.rm = TRUE) / (total_n - k_studies); pooled_sd <- sqrt(pooled_var) }
    return(list(mean = weighted_mean, sd = pooled_sd))
  }
  prior_info <- list(
    Orientation = synthesize_prior_info("Orientation"), Neuromotor = synthesize_prior_info("Neuromotor"),
    Motor = synthesize_prior_info("Motor"), Temperament = synthesize_prior_info("Temperament")
  )
  print("Synthesized Prior Knowledge for Control Group Outcomes:")
  print(prior_info)
}


# --- 1.2. Define Priors for Models ---
prior_orientation <- c(
  prior("normal(-1.5, 1.0)", class = "Intercept", coef = "1"), 
  prior("normal(-0.5, 1.0)", class = "Intercept", coef = "2"),
  prior("normal( 0.5, 1.0)", class = "Intercept", coef = "3"), 
  prior("normal( 1.5, 1.0)", class = "Intercept", coef = "4"),
  prior("normal(0, 0.8)", class = "b", coef = "Age_Weeks"),
  prior("normal(0, 0.4)", class = "b", coef = "Treatment_factorPIC"), 
  prior("normal(0, 0.8)", class = "b", coef = "Treatment_factorPIC:Age_Weeks"), 
  prior("normal(0, 1.5)", class = "b", coef = "Treatment_factorZika"), 
  prior("normal(0, 1.5)", class = "b", coef = "Treatment_factorZika:Age_Weeks"), 
  prior("student_t(3, 0, 2.5)", class = "sd")
)
prior_neuromotor <- prior_orientation
prior_motor <- prior_orientation
beta_sd_priors <- prior_orientation[!grepl("Intercept", prior_orientation$prior), ]
prior_temperament <- c(
  # Different Intercept priors based on synthesized mean for Temperament
  prior("normal(-1.0, 1.0)", class = "Intercept", coef = "1"), 
  prior("normal( 0.0, 1.0)", class = "Intercept", coef = "2"),
  prior("normal( 1.0, 1.0)", class = "Intercept", coef = "3"), 
  prior("normal( 2.0, 1.0)", class = "Intercept", coef = "4"),
  prior("normal( 3.0, 1.0)", class = "Intercept", coef = "5"),
  
  # The same beta and sd priors as the other models
  prior("normal(0, 0.8)", class = "b", coef = "Age_Weeks"),
  prior("normal(0, 0.4)", class = "b", coef = "Treatment_factorPIC"), 
  prior("normal(0, 0.8)", class = "b", coef = "Treatment_factorPIC:Age_Weeks"), 
  prior("normal(0, 1.5)", class = "b", coef = "Treatment_factorZika"), 
  prior("normal(0, 1.5)", class = "b", coef = "Treatment_factorZika:Age_Weeks"), 
  prior("student_t(3, 0, 2.5)", class = "sd")
)

# --- 1.3. Model Fitting Function ---
fit_brms_model <- function(data, dv_col, dv_levels, model_name, priors) {
  data_filtered <- data %>%
    filter(!is.na(.data[[dv_col]])) %>%
    mutate(DV_ord = factor(.data[[dv_col]], levels = dv_levels, ordered = TRUE))
  if (nrow(data_filtered) < 20 || nlevels(data_filtered$DV_ord) < 2) {
    print(paste("Skipping", model_name, "Model: Insufficient data or outcome levels."))
    return(NULL)
  }
  print(paste("--- Fitting", model_name, "Model with", nrow(data_filtered), "observations. ---"))
  formula_brms <- bf(DV_ord ~ Treatment_factor * Age_Weeks + (1 | Animal_ID_factor))
  model_fit <- brm(
    formula = formula_brms, data = data_filtered, family = cumulative("logit"),
    prior = priors, chains = 4, iter = 3000, warmup = 1000, seed = 12345,
    control = list(adapt_delta = 0.95), save_pars = save_pars(all = TRUE), silent = 1
  )
  return(model_fit)
}

# --- 1.4. Execute Model Fitting for All Variables ---
# This is the main time-consuming step.
model_orientation <- fit_brms_model(data_full, "Orientation_Scores_02", c(0,0.5,1,1.5,2), "Orientation Scores", prior_orientation)
summary(model_orientation, prob = 0.89)
#round(posterior_summary(model_orientation, probs = c(0.055, 0.945)), 3)

model_neuromotor  <- fit_brms_model(data_full, "Neuromotor_Scores_02", c(0,0.5,1,1.5,2), "Neuromotor Scores", prior_neuromotor)
summary(model_neuromotor, prob = 0.89)
#round(posterior_summary(model_orientation, probs = c(0.055, 0.945)), 3)

model_motor <- fit_brms_model(data_full, "Motor_Scores_02", c(0,0.5,1,1.5,2), "Motor Scores", prior_motor)
summary(model_motor, prob = 0.89)
#round(posterior_summary(model_orientation, probs = c(0.055, 0.945)), 3)

model_temperament <- fit_brms_model(data_full, "Temperament_Scores_02", c(0,0.5,1,1.5,2,2.5), "Temperament Scores", prior_temperament)
summary(model_temperament, prob = 0.89)
#round(posterior_summary(model_orientation, probs = c(0.055, 0.945)), 3)

# --- 1.5. Model Fitting Complete ---
# The model objects (e.g., model_orientation) now exist in your R environment.
# You can now save the entire workspace if you wish to close R and return later.
# For example:
# save.image("SNIP_analysis_workspace.RData")
# To load it back in a future session:
# load("SNIP_analysis_workspace.RData")
print("--- MODEL FITTING FINISHED ---")


################################################################################
# PART 2: POST-HOC CONTRAST ANALYSIS
# Run this section after Part 1 is complete. It assumes the model objects
# are loaded in the current R environment.
################################################################################

# --- 2.1. Contrast Analysis Function ---
perform_treatment_contrasts <- function(model, model_name, data_for_model) {
  if (is.null(model) || !inherits(model, "brmsfit")) {
    print(paste("Model object for", model_name, "is invalid or was not run. Skipping contrasts."))
    return(NULL)
  }
  print(paste("--- Performing Pairwise Contrasts for", model_name, "Model ---"))

  age_points <- unique(data_for_model$Age_Weeks)
  
  grid_df <- expand.grid(
    Treatment_factor = levels(data_for_model$Treatment_factor),
    Age_Weeks = age_points
  )
  
  linpred_draws_matrix <- posterior_linpred(model, newdata = grid_df, re_formula = NA)
  
  contrasts_df <- grid_df %>%
    mutate(draws = as.list(as.data.frame(linpred_draws_matrix))) %>%
    pivot_wider(names_from = Treatment_factor, values_from = draws) %>%
    ungroup() %>%
    mutate(
      `PIC - UIC` = map2(PIC, UIC, ~.x - .y),
      `Zika - UIC` = map2(Zika, UIC, ~.x - .y),
      `Zika - PIC` = map2(Zika, PIC, ~.x - .y)
    ) %>%
    select(Age_Weeks, `PIC - UIC`, `Zika - UIC`, `Zika - PIC`) %>%
    pivot_longer(cols = c(`PIC - UIC`, `Zika - UIC`, `Zika - PIC`), names_to = "contrast", values_to = "estimate_draws") %>%
    unnest(cols = c(estimate_draws))
  
  final_summary <- contrasts_df %>%
    group_by(Age_Weeks, contrast) %>%
    summarise(
      median_estimate = median(estimate_draws),
      lower.HPD = hdi(estimate_draws, .width = 0.89)[[1]],
      upper.HPD = hdi(estimate_draws, .width = 0.89)[[2]],
      .groups = "drop"
    )
  
  print(paste("Pairwise Contrasts on Latent Scale for", model_name, ":"))
  #print(as.data.frame(final_summary))
  
  return(final_summary)
}

# --- 2.2. Execute Contrast Analysis ---
# This section checks if the model objects exist in the environment before running.

# Analysis for Orientation Scores
if (exists("model_orientation") && !is.null(model_orientation)) {
  data_orientation_filtered <- data_full %>% filter(!is.na(Orientation_Scores_02))
  orientation_contrasts <- perform_treatment_contrasts(model_orientation, "Orientation Scores", data_orientation_filtered)
}
orientation_contrasts_rounded <- orientation_contrasts %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write_csv(orientation_contrasts_rounded, "orientation_contrasts_results.csv")

# Analysis for Neuromotor Scores
if (exists("model_neuromotor") && !is.null(model_neuromotor)) {
  data_neuromotor_filtered <- data_full %>% filter(!is.na(Neuromotor_Scores_02))
  neuromotor_contrasts <- perform_treatment_contrasts(model_neuromotor, "Neuromotor Scores", data_neuromotor_filtered)
}
neuromotor_contrasts_rounded <- neuromotor_contrasts %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write_csv(neuromotor_contrasts_rounded, "neuromotor_contrasts_results.csv")

# Analysis for Motor Scores
if (exists("model_motor") && !is.null(model_motor)) {
  data_motor_filtered <- data_full %>% filter(!is.na(Motor_Scores_02))
  motor_contrasts <- perform_treatment_contrasts(model_motor, "Motor Scores", data_motor_filtered)
}
motor_contrasts_rounded <- motor_contrasts %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write_csv(motor_contrasts_rounded, "motor_contrasts_results.csv")

# Analysis for Temperament Scores
if (exists("model_temperament") && !is.null(model_temperament)) {
  data_temperament_filtered <- data_full %>% filter(!is.na(Temperament_Scores_02))
  temperament_contrasts <- perform_treatment_contrasts(model_temperament, "Temperament Scores", data_temperament_filtered)
}
temperament_contrasts_rounded <- temperament_contrasts %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write_csv(temperament_contrasts_rounded, "temperament_contrasts_results.csv")
print("--- POST-HOC ANALYSIS FINISHED ---")


##### PLOTS ######

library(ggplot2)
library(readr)
library(dplyr)
library(patchwork) # For combining plots

# --- 2.3. Plotting Contrast Analysis ---

# Function to create a standardized plot for any given contrast result
plot_contrasts <- function(data, title) {
  ggplot(data, aes(x = Age_Weeks, y = median_estimate, ymin = lower.HPD, ymax = upper.HPD)) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
    geom_ribbon(aes(fill = contrast), alpha = 0.2) +
    geom_line(aes(color = contrast), linewidth = 1) +
    facet_wrap(~ contrast, scales = "free_y") +
    labs(
      title = title,
      #subtitle = "Showing the median difference and 89% HPDI. Dashed line at y=0 indicates no difference.",
      x = "Age (Weeks)",
      y = "Estimated Difference on Latent Scale"
    ) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    theme_bw() +
    theme(
      # --- Customizations for Text Size and Alignment ---
      
      # Center-align and increase size of the main plot title
      plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
      
      # From your original code, with added size
      strip.background = element_rect(fill = "grey90", color = "grey90"),
      strip.text = element_text(face = "bold", size = 16),
      
      # Increase the size of the axis titles (e.g., "Age (Weeks)")
      axis.title = element_text(size = 16),
      
      # Increase the size of the axis tick labels (the numbers on the axes)
      axis.text = element_text(size = 14),
      
      # Hide the legend since the facets are labeled
      legend.position = "none"
    )
}

# --- Create and Save Plots ---

# Load the data
orientation_contrasts <- read_csv("orientation_contrasts_results.csv")
neuromotor_contrasts <- read_csv("neuromotor_contrasts_results.csv")
motor_contrasts <- read_csv("motor_contrasts_results.csv")
temperament_contrasts <- read_csv("temperament_contrasts_results.csv")

# Generate the plots
p_orientation <- plot_contrasts(orientation_contrasts, "Orientation Scores")#: Treatment Contrasts Over Age")
p_neuromotor <- plot_contrasts(neuromotor_contrasts, "Neuromotor Scores")#: Treatment Contrasts Over Age")
p_motor <- plot_contrasts(motor_contrasts, "Motor Scores")#: Treatment Contrasts Over Age")
p_temperament <- plot_contrasts(temperament_contrasts, "State Control Scores")#: Treatment Contrasts Over Age")

# Display the plots
p_orientation
p_neuromotor
p_motor
p_temperament

# --- Optional: Save the plots to files ---
ggsave("orientation_contrasts_plot.png", p_orientation, width = 10, height = 4)
ggsave("neuromotor_contrasts_plot.png", p_neuromotor, width = 10, height = 4)
ggsave("motor_contrasts_plot.png", p_motor, width = 10, height = 4)
ggsave("temperament_contrasts_plot.png", p_temperament, width = 10, height = 4)

library(patchwork)

# Arrange the four plots into a single vertical column.
# The `/` operator stacks plots on top of each other.
# `plot_layout(axes = 'collect')` ensures the x-axis label only appears on the bottom plot.
combined_plot_vertical <- p_motor / p_neuromotor / p_orientation / p_temperament +
  plot_layout(axes = 'collect') +
  plot_annotation(title = "Neurobehavioral Assessment: Differences in Scores By Age", 
                  
                  theme = theme(
                    plot.title = element_text(
                      face = "bold",       # bold text
                      size = 18,           # larger font size
                      hjust = 0.5          # center alignment
                    )
                  )
  )



# Display the combined plot
print(combined_plot_vertical)


# Save the final vertical plot to a file.
ggsave("all_contrasts_vertical.png", combined_plot_vertical, width = 8, height = 10, dpi = 300)

print("--- PLOTTING FINISHED ---")

