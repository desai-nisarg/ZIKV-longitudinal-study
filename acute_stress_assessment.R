rm(list = ls())

############ Get weighted priors
library(dplyr)
library(readxl)

# Load data
df <- read_excel("HI_prior.xlsx", skip = 1)

# Rename columns for clarity
colnames(df) <- c("Publication", "Sample_Size", "Females", "Males", "Rearing", "Condition", 
                  "Freeze_Mean", "Freeze_SEM", "Hostile_Mean", "Hostile_SEM", 
                  "Anxious_Mean", "Anxious_SEM", "Coo_Mean", "Coo_SEM", 
                  "SelfDirected_Mean", "SelfDirected_SEM")

# Convert necessary columns to numeric
numeric_cols <- c("Sample_Size", "Freeze_Mean", "Freeze_SEM", "Hostile_Mean", "Hostile_SEM",
                  "Anxious_Mean", "Anxious_SEM", "Coo_Mean", "Coo_SEM", 
                  "SelfDirected_Mean", "SelfDirected_SEM")
df[numeric_cols] <- lapply(df[numeric_cols], as.numeric)

# Function to compute weighted mean
weighted_mean <- function(x, w) {
  sum(x * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

# Function to compute weighted SE (simple approach)
weighted_se <- function(se, w) {
  sum(se * w, na.rm = TRUE) / sum(w, na.rm = TRUE)
}

# Compute weighted means & SEs for each condition
weighted_stats <- df %>%
  group_by(Condition) %>%
  summarise(
    across(
      c(Freeze_Mean, Hostile_Mean, Anxious_Mean, SelfDirected_Mean), 
      list(wmean = ~ weighted_mean(.x, Sample_Size)), 
      .names = "wmean_{.col}"
    ),
    across(
      c(Freeze_SEM, Hostile_SEM, Anxious_SEM, SelfDirected_SEM), 
      list(wse = ~ weighted_se(.x, Sample_Size)), 
      .names = "wse_{.col}"
    )
  )

# View results
print(weighted_stats)




###############

# Load required libraries
library(brms)
library(emmeans)
library(ggplot2)
library(bayestestR)
library(tidybayes)
library(bayesplot)

# Read and prepare the data
data <- read.csv("human_intruder.csv")

# Ensure 'treatment' and 'sex' are factors
data$Group <- factor(data$Treatment.Group..0.UIC..1.PIC..2.Zika., levels = c(0,1,2), labels = c("Control", "Control", "Zika"))
data$Sex <- factor(data$Sex..0.f..1.m., levels = c(0,1), labels = c("Female", "Male"))
data$Self.Directed.Behavior..Dur..sec. <- round(data$Self.Directed.Behavior..Dur..sec.)
data$condition <- factor(data$Condition..0.alone..1.profile..2.stare., levels = c(0,1,2), labels = c("Alone", "Profile", "Stare"))
data$AnimalID <- as.factor(data$Animal.ID)
data$Freezing..Dur..sec. <- round(data$Freezing..Dur..sec.)

# Define behavior-specific information.
behavior_info <- list(
  list(
    var = "Freezing..Dur..sec.",
    control_means = weighted_stats$wmean_Freeze_Mean,
    control_ses   = weighted_stats$wse_Freeze_SEM
  ),
  list(
    var = "All.Hostile.Behavior..Freq.",
    control_means = weighted_stats$wmean_Hostile_Mean,
    control_ses   = weighted_stats$wse_Hostile_SEM
  ),
  list(
    var = "Anxious.Behavior..Freq.",
    control_means = weighted_stats$wmean_Anxious_Mean,
    control_ses   = weighted_stats$wse_Anxious_SEM
  ),
  list(
    var = "Self.Directed.Behavior..Dur..sec.",
    control_means = weighted_stats$wmean_SelfDirected_Mean,
    control_ses   = weighted_stats$wse_SelfDirected_SEM
  )
)

# Create an empty list to store fitted models.
model_list <- list()

# Loop over each behavior-specific info element.
for(info in behavior_info) {
  varname     <- info$var
  ctrl_means  <- info$control_means
  ctrl_ses    <- info$control_ses
  
  # Use Alone as the baseline.
  base_mean <- ctrl_means[1]
  base_se   <- ctrl_ses[1]
  intercept_mean <- log(base_mean)
  intercept_sd   <- base_se / base_mean
  
  # Compute the differences (contrasts) for the other conditions.
  delta2 <- log(ctrl_means[2]) - log(base_mean)
  delta2_sd <- sqrt( (ctrl_ses[2]/ctrl_means[2])^2 + (base_se/base_mean)^2 )
  
  delta3 <- log(ctrl_means[3]) - log(base_mean)
  delta3_sd <- sqrt( (ctrl_ses[3]/ctrl_means[3])^2 + (base_se/base_mean)^2 )
  
  # Define priors.
  priors <- c(
    set_prior(sprintf("normal(%f, %f)", intercept_mean, intercept_sd), class = "Intercept"),
    set_prior(sprintf("normal(%f, %f)", delta2, delta2_sd), class = "b", coef = "conditionProfile"),
    set_prior(sprintf("normal(%f, %f)", delta3, delta3_sd), class = "b", coef = "conditionStare"),
    set_prior("normal(0, 1)", class = "b", coef = "SexMale"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika:SexMale"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika:conditionProfile"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika:conditionStare"),
    set_prior("normal(0, 1)", class = "b", coef = "SexMale:conditionProfile"),
    set_prior("normal(0, 1)", class = "b", coef = "SexMale:conditionStare"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika:SexMale:conditionProfile"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika:SexMale:conditionStare"),
    set_prior("exponential(1)", class = "shape")
  )
  
  family_spec <- negbinomial()
  model_formula <- as.formula(paste(varname, "~ Group * Sex * condition + (1|AnimalID)"))
  
  message(sprintf("\nFitting model for '%s'...", varname))
  
  fit_model <- brm(
    formula = model_formula,
    data = data,
    family = family_spec,
    prior = priors,
    iter = 4000,
    warmup = 1000,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.95)
  )
  
  model_list[[varname]] <- fit_model
  print(summary(fit_model))
  
  # Model diagnostics
  message(sprintf("\n--- Diagnostics for '%s' ---", varname))
  
  # 1. Check Rhat (should all be < 1.01)
  rhat_vals <- brms::rhat(fit_model)
  max_rhat <- max(rhat_vals, na.rm = TRUE)
  message(sprintf("Max Rhat: %.4f (should be < 1.01)", max_rhat))
  if(max_rhat > 1.01) {
    warning(sprintf("Convergence issues detected for %s! Max Rhat = %.4f", varname, max_rhat))
  }
  
  # 2. Check for divergent transitions
  div_trans <- sum(nuts_params(fit_model)$Value[nuts_params(fit_model)$Parameter == "divergent__"])
  message(sprintf("Divergent transitions: %d (should be 0)", div_trans))
  if(div_trans > 0) {
    warning(sprintf("%s has %d divergent transitions!", varname, div_trans))
  }
  
  # 3. Posterior predictive check
  message("Generating posterior predictive check plot...")
  pp <- pp_check(fit_model, ndraws = 100) + 
    ggtitle(paste("Posterior Predictive Check:", varname))
  print(pp)
  
  # 4. Trace plots for key parameters
  #message("Generating trace plots...")
  #trace_plot <- plot(fit_model, pars = c("b_Intercept", "b_GroupZika", "b_SexMale"))
  #print(trace_plot)
}


#### Models without priors
behavior_info_no_prior <- list(
  list(var = "Affiliative.Vocals..Freq.", type = "frequency"),
  list(var = "Screams..Freq.", type = "frequency"),
  list(var = "Threat.Barks..Freq.", type = "frequency"),
  list(var = "Fearful.Behavior..Freq.", type = "frequency"),
  list(var = "Threats.Toward.Intruder..Freq.", type = "frequency"),
  list(var = "Threats.Away.From.Intruder..Freq.", type = "frequency")
)

model_list_no_prior <- list()

for(info in behavior_info_no_prior) {
  varname <- info$var
  
  priors_no_prior <- c(
    set_prior("normal(0, 5)", class = "Intercept"),
    set_prior("normal(0, 1)", class = "b", coef = "conditionProfile"),
    set_prior("normal(0, 1)", class = "b", coef = "conditionStare"),
    set_prior("normal(0, 1)", class = "b", coef = "SexMale"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika:SexMale"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika:conditionProfile"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika:conditionStare"),
    set_prior("normal(0, 1)", class = "b", coef = "SexMale:conditionProfile"),
    set_prior("normal(0, 1)", class = "b", coef = "SexMale:conditionStare"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika:SexMale:conditionProfile"),
    set_prior("normal(0, 1)", class = "b", coef = "GroupZika:SexMale:conditionStare"),
    set_prior("exponential(1)", class = "shape")
  )
  
  family_spec <- negbinomial()
  model_formula <- as.formula(paste(varname, "~ Group * Sex * condition + (1|AnimalID)"))
  
  message(sprintf("\nFitting model for '%s' with generic priors...", varname))
  
  fit_model <- brm(
    formula = model_formula,
    data = data,
    family = family_spec,
    prior = priors_no_prior,
    iter = 4000,
    warmup = 1000,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.95)
  )
  
  model_list_no_prior[[varname]] <- fit_model
  print(summary(fit_model))
  
  # Model diagnostics
  message(sprintf("\n--- Diagnostics for '%s' ---", varname))
  
  # 1. Check Rhat
  rhat_vals <- brms::rhat(fit_model)
  max_rhat <- max(rhat_vals, na.rm = TRUE)
  message(sprintf("Max Rhat: %.4f (should be < 1.01)", max_rhat))
  if(max_rhat > 1.01) {
    warning(sprintf("Convergence issues detected for %s! Max Rhat = %.4f", varname, max_rhat))
  }
  
  # 2. Check for divergent transitions
  div_trans <- sum(nuts_params(fit_model)$Value[nuts_params(fit_model)$Parameter == "divergent__"])
  message(sprintf("Divergent transitions: %d (should be 0)", div_trans))
  if(div_trans > 0) {
    warning(sprintf("%s has %d divergent transitions!", varname, div_trans))
  }
  
  # 3. Posterior predictive check
  message("Generating posterior predictive check plot...")
  pp <- pp_check(fit_model, ndraws = 100) + 
    ggtitle(paste("Posterior Predictive Check:", varname))
  print(pp)
  
  # 4. Trace plots
  #message("Generating trace plots...")
  #trace_plot <- plot(fit_model, pars = c("b_Intercept", "b_GroupZika", "b_SexMale"))
  #print(trace_plot)
}


###############################################################################
# COMPREHENSIVE CONTRAST ANALYSIS
###############################################################################

# Function to perform comprehensive contrasts for each model
perform_contrasts <- function(model, model_name, ci_level = 0.89) {
  
  message(sprintf("\n========================================"))
  message(sprintf("CONTRASTS FOR: %s", model_name))
  message(sprintf("========================================\n"))
  
  # Helper function to back-transform and interpret effect sizes
  interpret_effect <- function(contrast_summary) {
    # Extract the estimate and CI bounds (assuming log scale)
    est <- contrast_summary$estimate
    lower <- contrast_summary$lower.HPD
    upper <- contrast_summary$upper.HPD
    
    # Back-transform to ratio scale
    ratio_est <- exp(est)
    ratio_lower <- exp(lower)
    ratio_upper <- exp(upper)
    
    # Calculate percent change
    pct_change <- (ratio_est - 1) * 100
    pct_lower <- (ratio_lower - 1) * 100
    pct_upper <- (ratio_upper - 1) * 100
    
    # Create interpretation
    if (pct_change > 0) {
      direction <- "increase"
    } else {
      direction <- "decrease"
    }
    
    interp <- sprintf(
      "  → Effect size: %.1f%% %s [%.1f%%, %.1f%%]\n  → Rate ratio: %.2f [%.2f, %.2f]",
      abs(pct_change), direction, abs(pct_lower), abs(pct_upper),
      ratio_est, ratio_lower, ratio_upper
    )
    
    return(interp)
  }
  
  # 1. Main effect of Group (Control vs Zika)
  message("--- Main Group Effect (Control vs Zika) ---")
  emm_group <- emmeans(model, ~ Group)
  contrast_group <- pairs(emm_group, level = ci_level, reverse = T)
  print(contrast_group)
  
  # Add interpretation for Control - Zika contrast
  contrast_df <- summary(contrast_group)
  control_zika <- contrast_df[contrast_df$contrast == "Control - Zika", ]
  if(nrow(control_zika) > 0) {
    message(interpret_effect(control_zika))
  }
  
  # 2. Group effect within each Sex
  message("\n--- Group Effect by Sex (Control vs Zika | Sex) ---")
  emm_group_sex <- emmeans(model, ~ Group | Sex)
  contrast_group_sex <- pairs(emm_group_sex, level = ci_level, reverse = T)
  print(contrast_group_sex)
  
  contrast_df_sex <- summary(contrast_group_sex)
  for(sex_level in unique(contrast_df_sex$Sex)) {
    sex_contrast <- contrast_df_sex[contrast_df_sex$Sex == sex_level & 
                                      contrast_df_sex$contrast == "Zika - Control", ]
    if(nrow(sex_contrast) > 0) {
      message(sprintf("  %s: %s", sex_level, interpret_effect(sex_contrast)))
    }
  }
  
  # 3. Group effect within each Condition
  message("\n--- Group Effect by Condition (Control vs Zika | Condition) ---")
  emm_group_cond <- emmeans(model, ~ Group | condition)
  contrast_group_cond <- pairs(emm_group_cond, level = ci_level, reverse = T)
  print(contrast_group_cond)
  
  contrast_df_cond <- summary(contrast_group_cond)
  for(cond_level in unique(contrast_df_cond$condition)) {
    cond_contrast <- contrast_df_cond[contrast_df_cond$condition == cond_level & 
                                        contrast_df_cond$contrast == "Zika - Control", ]
    if(nrow(cond_contrast) > 0) {
      message(sprintf("  %s: %s", cond_level, interpret_effect(cond_contrast)))
    }
  }
  
  # 4. Group effect within each Sex x Condition combination
  message("\n--- Group Effect by Sex x Condition (Control vs Zika | Sex + Condition) ---")
  emm_group_sex_cond <- emmeans(model, ~ Group | Sex + condition)
  contrast_group_sex_cond <- pairs(emm_group_sex_cond, level = ci_level, reverse = T)
  print(contrast_group_sex_cond)
  
  # 5. Sex effect (overall)
  message("\n--- Main Sex Effect (Female vs Male) ---")
  emm_sex <- emmeans(model, ~ Sex)
  contrast_sex <- pairs(emm_sex, level = ci_level)
  print(contrast_sex)
  
  # 6. Sex effect within each Condition
  message("\n--- Sex Effect by Condition (Female vs Male | Condition) ---")
  emm_sex_cond <- emmeans(model, ~ Sex | condition)
  contrast_sex_cond <- pairs(emm_sex_cond, level = ci_level)
  print(contrast_sex_cond)
  
  # 7. Condition effects (overall)
  message("\n--- Main Condition Effects (Pairwise) ---")
  emm_cond <- emmeans(model, ~ condition)
  contrast_cond <- pairs(emm_cond, level = ci_level)
  print(contrast_cond)
  
  # 8. Condition effects within each Group
  message("\n--- Condition Effects by Group (Pairwise | Group) ---")
  emm_cond_group <- emmeans(model, ~ condition | Group)
  contrast_cond_group <- pairs(emm_cond_group, level = ci_level)
  print(contrast_cond_group)
  
  # Return all contrasts as a list for potential further use
  contrasts_list <- list(
    group = contrast_group,
    group_by_sex = contrast_group_sex,
    group_by_condition = contrast_group_cond,
    group_by_sex_condition = contrast_group_sex_cond,
    sex = contrast_sex,
    sex_by_condition = contrast_sex_cond,
    condition = contrast_cond,
    condition_by_group = contrast_cond_group
  )
  
  return(invisible(contrasts_list))
}

# Apply comprehensive contrasts to all models with priors
message("\n\n###############################################################################")
message("COMPREHENSIVE CONTRASTS: MODELS WITH INFORMATIVE PRIORS")
message("###############################################################################\n")

contrast_results <- list()
for(info in behavior_info) {
  varname <- info$var
  contrast_results[[varname]] <- perform_contrasts(
    model_list[[varname]], 
    varname,
    ci_level = 0.89
  )
}

# Apply comprehensive contrasts to all models without priors
message("\n\n###############################################################################")
message("COMPREHENSIVE CONTRASTS: MODELS WITHOUT INFORMATIVE PRIORS")
message("###############################################################################\n")

contrast_results_no_prior <- list()
for(info in behavior_info_no_prior) {
  varname <- info$var
  contrast_results_no_prior[[varname]] <- perform_contrasts(
    model_list_no_prior[[varname]], 
    varname,
    ci_level = 0.89
  )
}


################################################################################

library(ggplot2)
library(tidybayes)
library(dplyr)
library(forcats)
library(patchwork) 

message("\n\nCreating Refined 2-Panel Visualization...")

# ==============================================================================
# 1. SETUP LISTS & HELPERS (Same as before)
# ==============================================================================

# Panel 1: Overall Group Effect
p1_behaviors <- list(
  "Freezing..Dur..sec."               = "prior",
  "All.Hostile.Behavior..Freq."       = "prior",
  "Anxious.Behavior..Freq."           = "prior",
  "Affiliative.Vocals..Freq."         = "no_prior",
  "Fearful.Behavior..Freq."           = "no_prior",
  "Self.Directed.Behavior..Dur..sec." = "prior",
  "Screams..Freq."      = "no_prior"
)

# Panel 2: Group Effect Split by Sex
p2_behaviors <- list(
  "Threat.Barks..Freq." = "no_prior"
)

# Helper function to clean names
clean_behavior_names <- function(df) {
  df %>%
    mutate(
      behavior_clean = gsub("\\.\\.", " ", behavior),
      behavior_clean = gsub("\\.", " ", behavior_clean),
      behavior_clean = trimws(behavior_clean),
      behavior_clean = gsub(" Freq$", "", behavior_clean),
      behavior_clean = gsub(" Dur sec$", "", behavior_clean),
      behavior_clean = gsub("All ", "", behavior_clean),
      behavior_clean = case_when(
        behavior_clean == "Affiliative Vocals" ~ "Affiliative Vocals",
        TRUE ~ behavior_clean
      )
    )
}

# ==============================================================================
# 2. EXTRACT DATA FOR PANEL 1 (Main Group Effect)
# ==============================================================================
draws_p1 <- list()

for (varname in names(p1_behaviors)) {
  source_type <- p1_behaviors[[varname]]
  
  if (source_type == "prior") {
    res_list <- contrast_results
  } else {
    res_list <- contrast_results_no_prior
  }
  
  # IMPORTANT: Use 'group_by_condition' (Main effect) instead of 'group_by_sex_condition'
  if (!is.null(res_list[[varname]]$group_by_condition)) {
    draws <- res_list[[varname]]$group_by_condition %>%
      gather_emmeans_draws() %>%
      filter(contrast == "Zika - Control") %>%
      mutate(behavior = varname)
    
    draws_p1[[varname]] <- draws
  }
}

df_p1 <- bind_rows(draws_p1) %>%
  clean_behavior_names() %>%
  mutate(behavior_label = factor(behavior_clean, levels = c(
    "Freezing", "Hostile Behavior", "Anxious Behavior", 
    "Affiliative Vocals", "Fearful Behavior", "Self Directed Behavior", "Screams"
  )))

# ==============================================================================
# 3. EXTRACT DATA FOR PANEL 2 (Sex Interaction)
# ==============================================================================
draws_p2 <- list()

for (varname in names(p2_behaviors)) {
  source_type <- p2_behaviors[[varname]]
  
  if (source_type == "prior") {
    res_list <- contrast_results
  } else {
    res_list <- contrast_results_no_prior
  }
  
  # Use 'group_by_sex_condition' here
  if (!is.null(res_list[[varname]]$group_by_sex_condition)) {
    draws <- res_list[[varname]]$group_by_sex_condition %>%
      gather_emmeans_draws() %>%
      filter(contrast == "Zika - Control") %>%
      mutate(behavior = varname)
    
    draws_p2[[varname]] <- draws
  }
}

df_p2 <- bind_rows(draws_p2) %>%
  clean_behavior_names() %>%
  mutate(behavior_label = factor(behavior_clean, levels = c("Threat Barks")))

# ==============================================================================
# 4. PLOTTING
# ==============================================================================

# Shared Theme Elements
common_theme <- theme_classic(base_size = 14) +
  theme(
    axis.text.x = element_text(angle = 35, hjust = 1, size = 11),
    strip.text = element_text(size = 12, face = "bold"),
    strip.background = element_rect(fill = "grey90", color = "grey90"),
    axis.title = element_text(size = 14),
    panel.grid.major.y = element_line(color = "grey95"),
    
    # NEW: Bold, smaller headings for panels
    plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
    plot.tag = element_text(face = "bold") # Ensure (a)/(b) are bold
  )

# --- Left Panel: Overall Effects ---
p_left <- ggplot(df_p1, aes(y = .value, x = behavior_label)) +
  stat_halfeye(
    fill = "#8dd3c7", # NEW: Set color to teal
    alpha = 0.8,
    .width = c(0.66, 0.89),
    interval_color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  facet_wrap(~ condition, nrow = 3, scales = "free_y") +
  labs(
    title = "Overall Effect by Group",
    #subtitle = "Positive values indicate Zika > Control",
    y = "Estimated Difference (log scale)",
    x = NULL
  ) +
  coord_cartesian(ylim = c(-3, 4)) +
  common_theme

# --- Right Panel: Sex-Specific Effects ---
p_right <- ggplot(df_p2, aes(y = .value, x = behavior_label, fill = Sex)) +
  stat_halfeye(
    position = position_dodge(width = 0.8),
    alpha = 0.7,
    .width = c(0.66, 0.89),
    interval_color = "black"
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  facet_wrap(~ condition, nrow = 3, scales = "free_y") +
  scale_fill_manual(values = c("Female" = "#fb8072", "Male" = "#80b1d3")) +
  labs(
    title = "Effect by Sex",
    #subtitle = "By Sex (Zika - Control)",
    y = NULL, 
    x = NULL
  ) +
  coord_cartesian(ylim = c(-3, 4)) +
  common_theme +
  theme(legend.position = "top")

# ==============================================================================
# 5. COMBINE AND SAVE (CORRECTED)
# ==============================================================================

# We combine the plots and define all annotation themes (Title + Tags) 
# inside the plot_annotation() function to avoid the "&" operator error.

final_plot <- (p_left | p_right) + 
  plot_layout(widths = c(3, 1.2)) + 
  plot_annotation(
    title = "Acute Stress Assessment: ZIKV - Control Contrasts",
    subtitle = "Positive values indicate Zika > Control",
    tag_levels = 'A', # Adds (a) and (b)
    theme = theme(
      # Style the Main Title
      plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
      plot.subtitle = element_text(size = 11, hjust = 0.5, color = "gray40"),
      # Style the (a)/(b) Tags
      plot.tag = element_text(size = 14, face = "bold") 
    )
  )

print(final_plot)

ggsave("zika_contrasts_final.png", final_plot, width = 16, height = 10, dpi = 300)
message("Saved zika_contrasts_final.png")





###############################################################################
# VISUALIZATION: ZIKA-CONTROL CONTRASTS BY SEX & CONDITION
###############################################################################
#
# This plot shows the Zika - Control contrast for all behaviors,
# faceted by the condition (Alone, Profile, Stare).
# Within each plot, the x-axis shows the behavior (and its prior type),
# and the fill color shows the sex.
#
###############################################################################

library(ggplot2)
library(tidybayes)
library(dplyr)
library(forcats) # For fct_reorder

message("\n\nCreating visualization: Zika-Control Contrasts by Sex, Condition, and Behavior...")

# 1. Create an empty list to store all contrast draws
all_behavior_contrasts <- list()

# 2. Loop through models WITH priors
message("Extracting contrasts from models with informative priors...")
for (varname in names(contrast_results)) {
  # Get the contrast object for Group by Sex x Condition
  contrast_obj <- contrast_results[[varname]]$group_by_sex_condition
  
  # Gather draws, filter for "Zika - Control", and add metadata
  draws <- contrast_obj %>%
    gather_emmeans_draws() %>%
    filter(contrast == "Zika - Control") %>%
    mutate(behavior = varname, prior_type = "Informative")
  
  all_behavior_contrasts[[varname]] <- draws
}

# 3. Loop through models WITHOUT priors
message("Extracting contrasts from models without informative priors...")
for (varname in names(contrast_results_no_prior)) {
  # Get the contrast object for Group by Sex x Condition
  contrast_obj <- contrast_results_no_prior[[varname]]$group_by_sex_condition
  
  # Gather draws, filter for "Zika - Control", and add metadata
  draws <- contrast_obj %>%
    gather_emmeans_draws() %>%
    filter(contrast == "Zika - Control") %>%
    mutate(behavior = varname, prior_type = "Non-Informative")
  
  all_behavior_contrasts[[paste0(varname, "_no_prior")]] <- draws
}

# 4. Combine into one big data frame
plot_data_behaviors <- bind_rows(all_behavior_contrasts)

# 4.5. Define the desired plotting order
# These are the "clean" names that will be produced in Step 5
behavior_order <- c(
  "Freezing",
  "Hostile Behavior",
  "Anxious Behavior",
  "Threat Barks",
  "Screams",
  "Affiliative Vocals",
  "Fearful Behavior",
  "Self Directed Behavior"
)

# 5. Clean up behavior names for plotting
plot_data_behaviors_cleaned <- plot_data_behaviors %>%
  mutate(
    # Clean up the variable names
    behavior_clean = gsub("\\.\\.", " ", behavior),      # 'All.Hostile.Behavior..Freq.' -> 'All.Hostile.Behavior Freq '
    behavior_clean = gsub("\\.", " ", behavior_clean),  # 'All Hostile Behavior Freq '
    behavior_clean = trimws(behavior_clean),            # 'All Hostile Behavior Freq' (removes trailing space)
    
    # Now, remove the suffixes from the end of the string
    behavior_clean = gsub(" Freq$", "", behavior_clean),
    behavior_clean = gsub(" Dur sec$", "", behavior_clean),
    
    # Clean up prefixes and specific terms
    behavior_clean = gsub("All ", "", behavior_clean)
  ) %>%
  
  # NEW: Filter out the behaviors you no longer want
  filter(behavior_clean %in% behavior_order) %>%
  
  # NEW: Create the final x-axis label as an ordered factor
  mutate(
    behavior_label = factor(behavior_clean, levels = behavior_order)
  )

# 6. Create the plot
message("Generating plot...")

p_zika_by_sex_cond <- ggplot(plot_data_behaviors_cleaned, aes(y = .value, x = behavior_label, fill = Sex)) +
  stat_halfeye(
    position = position_dodge(width = 0.8), # Dodge bars for each sex
    alpha = 0.7,
    .width = c(0.66, 0.89), # 66% and 89% CIs
    interval_color = "black",
    slab_color = NA
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  
  # Facet by condition in 3 rows, as requested
  facet_wrap(~ condition, nrow = 3, scales = "free_y") + 
  
  scale_fill_manual(values = c("Female" = "#fb8072", "Male" = "#80b1d3"), name = "Sex") +
  labs(
    title = "Acute Stress Assessment: ZIKV - Control Contrasts by Sex and Condition",
    subtitle = "Positive values indicate Zika > Control",
    y = "Estimated Difference (log scale)",
    x = "Behavior" # <-- MODIFIED
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 35, hjust = 1, size = 12), # <-- MODIFIED
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    axis.text.y = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40"),
    strip.text = element_text(size = 14, face = "bold"), # Make facet labels clearer
    strip.background = element_rect(fill = "grey90", color = "grey90") # Darker background for labels
  ) +
  coord_cartesian(ylim = c(-3, 4)) # Set a reasonable y-limit for log-scale contrasts

print(p_zika_by_sex_cond)

# 7. Save the plot
ggsave("zika_contrasts_by_sex_condition_all_behaviors.png", p_zika_by_sex_cond, width = 16, height = 12, dpi = 300)
message("Saved zika_contrasts_by_sex_condition_all_behaviors.png")

