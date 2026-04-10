rm(list = ls())

# Load required libraries
library(brms)
library(emmeans)
library(bayesplot)
library(tidybayes)
library(ggplot2)
library(ggridges)
library(tidyverse)
library(parameters)
library(patchwork)
library(cowplot)

options(mc.cores = parallel::detectCores())
set.seed(12)


# =============================================================================
# Fit Models for Each Behavior with Behavior-Specific Priors
# =============================================================================

# Read and prepare the data
data <- read.csv("attachment.csv")

# Ensure 'treatment' and 'sex' are factors
data$Group <- factor(data$Treatment.Group..0.UIC..1.PIC..2.Zika., levels = c(0,1,2), labels = c("Control", "Control", "Zika"))
data$Sex <- factor(data$Sex..0.f..1.m., levels = c(0,1), labels = c("Female", "Male"))
data$Self.Directed.Behaviors..Dur..sec. <- round(data$Self.Directed.Behaviors..Dur..sec.)

# We'll store the models in a list for easy management
models <- list()
behavior_vars <- c(names(data[7:14]))
### ---- A. Frequency Models (Negative Binomial) ---- ###

# For count data with a log link, we need to convert control means and SEs:
# Delta method: log-transform approximate SE = SE_raw / mean_raw

# ----Scream----
control_mean_scream <- 32       # Control group raw mean for behavior 1
control_se_scream   <- 22       # Control standard error

mu_log_scream <- log(control_mean_scream)            
se_log_scream <- control_se_scream / control_mean_scream   

cat("Scream - Control (log-scale): mean =", round(mu_log_scream, 3),
    ", SE =", round(se_log_scream, 3), "\n")


prior_scream <- c(
  set_prior(sprintf("normal(%f, %f)", mu_log_scream, se_log_scream), class = "Intercept"),
  set_prior("normal(0, 1)", class = "b", coef = "GroupZika"),
  set_prior("normal(0, 1)", class = "b", coef = "SexMale"),
  set_prior("normal(0, 1)", class = "b", coef = "GroupZika:SexMale"),
  set_prior("exponential(1)", class = "shape")
)

model_freq_scream <- brm(
  formula = Scream..Freq. ~ Group * Sex,
  data = data,
  family = negbinomial(),    
  prior = prior_scream,
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  control = list(adapt_delta = 0.95)
)
models$freq_scream <- model_freq_scream
pp_check(model_freq_scream) + scale_x_log10()
summary(model_freq_scream, prob = 0.89)


# ----Remaining behaviors
behavior_vars <- c(names(data[7:14]))

for (behavior in behavior_vars) {
  formula <- as.formula(
    paste0(behavior, " ~ Group * Sex")
  )
  
  models[[behavior]] <- brm(
    formula = formula,
    data = data,
    family = negbinomial(),  
    iter = 4000,
    warmup = 1000,
    chains = 4,
    cores = 4,
    control = list(adapt_delta = 0.95),
    seed = 111
  )
}

posterior_summary(models[[4]], probs = c(0.055, 0.945))
posterior_summary(models[[5]], probs = c(0.055, 0.945))
posterior_summary(models[[7]], probs = c(0.055, 0.945))

for (i in 1:9){
  print(summary(models[[i]], prob = 0.89))
}

# For Fearful Behaviors
fearful_emms <- emmeans(models[[4]], ~ Group | Sex)
fearful_contrasts <- contrast(fearful_emms, method = "pairwise", level = 0.89)
print(fearful_contrasts)

# For Hostile/Defensive Behaviors
hostile_emms <- emmeans(models[[5]], ~ Group | Sex)
hostile_contrasts <- contrast(hostile_emms, method = "pairwise", level = 0.89)
print(hostile_contrasts)

# For Self-Directed Behaviors
selfdirected_emms <- emmeans(models[[7]], ~ Group | Sex)
selfdirected_contrasts <- contrast(selfdirected_emms, method = "pairwise", level = 0.89)
print(selfdirected_contrasts)

# For Screams
screams_emms <- emmeans(models[[1]], ~ Group | Sex)
screams_contrasts <- contrast(screams_emms, method = "pairwise", level = 0.89)
print(screams_contrasts)


### ---- B. Index of Preference Model (Beta Regression with logit link) ---- ###

# Transform the index from [-1,1] to (0,1)
data$pref_adj <- (data$Index.of.Preference..Scale...1.to.1. + 1) / 2
data$pref_adj <- pmin(pmax(data$pref_adj, 0.001), 0.999)  # trim limits to avoid 0/1

# Suppose that for the control group the transformed mean is 0.5 with SE = 0.1.
pref_mean_raw <- 0.9
pref_se_raw   <- 0.532

# Define the logit function and transform:
logit <- function(p) log(p / (1 - p))
pref_mu_logit <- logit(pref_mean_raw)  
# Delta method: derivative of logit() at 0.9 is 1/(0.9 * 0.1) = 11.111.
pref_se_logit <- pref_se_raw * 11.111          

cat("Index of Preference - Control (logit-scale): mean =", round(pref_mu_logit, 3),
    ", SE =", round(pref_se_logit, 3), "\n")

prior_index <- c(
  set_prior(sprintf("normal(%f, %f)", pref_mu_logit, pref_se_logit), class = "Intercept"),
  set_prior("normal(0, 1)", class = "b", coef = "GroupZika"),
  set_prior("normal(0, 1)", class = "b", coef = "SexMale"),
  set_prior("normal(0, 1)", class = "b", coef = "GroupZika:SexMale")
)

model_index <- brm(
  formula = pref_adj ~ Group * Sex,
  data = data,
  family = Beta(link = "logit"),
  prior = prior_index,
  iter = 4000,
  warmup = 1000,
  chains = 4,
  cores = 4,
  control = list(adapt_delta = 0.95)
)
models$index <- model_index
summary(model_index)
posterior_interval(model_index, prob = 0.89)
parameters::model_parameters(model_index, ci = 0.89) 
#The model has a log- or logit-link. Consider using `exponentiate = TRUE` to interpret coefficients as ratios.

coefs <- posterior_summary(model_index, probs = c(0.055, 0.945))
coefs <- round(coefs, 3)
write.csv(coefs, "index_model.csv", row.names = TRUE)

# =============================================================================
# Model Summaries, Pairwise Contrasts, and Posterior Plots
# =============================================================================

# Print summaries for each model
for(name in names(models)) {
  cat("\n------ Summary for model:", name, "------\n")
  print(model_parameters(models[[name]], ci = 0.89, ci_method = "HDI"))
}

# Perform pairwise contrasts (by sex) for each frequency model using emmeans

contrast_results <- list()

for(name in names(models)) {
  emm_obj <- emmeans(models[[name]], ~ Group | Sex)
  contrast_results[[name]] <- pairs(emm_obj, reverse = T, level = 0.89)
  cat("\nPairwise Contrasts for", name, ":\n")
  print(contrast_results[[name]])
}

##### Create table of significant effects.

# Load necessary libraries
library(dplyr)
library(purrr)

# 1. Combine all your contrast results into a single data frame
all_contrasts <- map_dfr(names(contrast_results), function(behavior) {
  # Convert the emmeans object to a standard data frame
  df <- as.data.frame(contrast_results[[behavior]])
  # Add the behavior name as a column
  df$Behavior <- behavior
  return(df)
})

# 2. Filter for "significant" effects and clean up the table
significant_table <- all_contrasts %>%
  # Filter to keep only rows where the 89% HPD interval does NOT cross 0
  filter(lower.HPD > 0 | upper.HPD < 0) %>%
  # Select and rename the columns to match your desired output
  select(
    Behavior,
    Contrast = contrast,
    Sex,
    `Posterior median` = estimate,
    `Lower 89% HPD` = lower.HPD,
    `Upper 89% HPD` = upper.HPD
  )

# View the result in your R console
print(significant_table)

# 3. Export to CSV so you can paste it into your manuscript/report
write.csv(significant_table, "attachment_behavior_contrasts.csv", row.names = FALSE)



# Load libraries if you are running this in a new session
library(emmeans)
library(tidybayes)
library(ggplot2)
library(dplyr)

# --- 1. Prepare Data from Multiple Models ---

# List the specific behaviors and their model names you want to plot
behaviors_to_plot <- list(
  "Index of Preference" = models[[10]],
  "Screams" = models[[1]],
  "Attention Seeking" = models[[2]],
  "Affiliative" = models[[3]],
  "Fearful" = models[[4]],
  "Hostile" = models[[5]],
  "Anxious" = models[[6]],
  "Self-Directed" = models[[7]]
)

# Use purrr::map_dfr to loop, calculate contrasts, and combine into one data frame
all_contrasts <- purrr::map_dfr(behaviors_to_plot, ~{
  emmeans(.x, pairwise ~ Group | Sex) %>%
    pluck("contrasts") %>% # Extract just the contrast results
    gather_emmeans_draws() %>%
    filter(contrast == "Control - Zika")
}, .id = "behavior") # .id creates the 'behavior' column from the list names

all_contrasts$Z_C = -all_contrasts$.value

# --- Alternative desired_order vector (if you want to include "Anxious") ---

desired_order_with_anxious <- c(
  "Index of Preference",
  "Attention Seeking",  
  "Affiliative",
  "Screams",
  "Hostile",
  "Self-Directed",
  "Fearful",
  "Anxious"
)

# Then you would use this new vector in the mutate() call:
all_contrasts <- all_contrasts %>%
  mutate(behavior = factor(behavior, levels = desired_order_with_anxious))

# --- 2. Create the Single Combined Plot ---

p_combined_behaviors <- ggplot(all_contrasts, aes(y = Z_C, x = behavior, fill = Sex)) +
  # Use position_dodge() to place the half-eye plots for each sex side-by-side
  stat_halfeye(
    position = position_dodge(width = 0.75),
    alpha = 0.7,
    .width = c(0.66, 0.89),
    interval_color = "black",
    slab_color = NA
  ) + ylim(-8, 4) +
  # Add a reference line at zero
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray30") +
  
  # Set custom colors and a title for the legend
  scale_fill_manual(values = c("Female" = "#fb8072", "Male" = "#80b1d3"), name = "Sex") +
  
  labs(
    title = "Attachment Assessment: ZIKV - Control Contrast Across Behaviors",
    subtitle = "Positive values indicate Zika > Control",
    y = "Estimated Difference (log scale)",
    x = "Behavior Category"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.position = "top",
    axis.text.x = element_text(angle = 15, hjust = 0.7), # Angle text slightly
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 14),
    plot.title = element_text(size = 18, face = "bold", hjust = 0.5),
    plot.subtitle = element_text(size = 12, hjust = 0.5, color = "gray40")
  )

print(p_combined_behaviors)