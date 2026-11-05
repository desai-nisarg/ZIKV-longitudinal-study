rm(list=ls())

# Load required packages
library(brms)
library(dplyr)
library(ggplot2)
library(emmeans)
library(tidybayes)
library(tidyverse)

# Read in the data
data <- read.csv("visual_acuity.csv", stringsAsFactors = FALSE)


# Convert numeric factors to factors
data <- data %>%
  mutate(
    #Group = factor(Treatment.Group..0.UIC..1.PIC..2.Zika., levels = c(0,1,2), labels = c("UIC", "PIC", "Zika")),
    Group = factor(Treatment.Group..0.UIC..1.PIC..2.Zika., levels = c(0,1,2), labels = c("Control", "Control", "Zika")), 
    Contrast_level = factor(Contrast..1.low..2.medium..3.high., levels = c(1,2,3), labels = c("Low", "Medium", "High")),
    Age = factor(Age..months., levels = c(4,6), labels = c("4 Months", "6 Months")),
    Sex = factor(Sex..0.f..1.m., levels = c(0,1), labels = c("Female", "Male")),
    left_right_acuity = factor(Left.or.Right.Acuity..0.left..1.right., levels = c(0,1), labels = c("Left", "Right")),
    acuity_number = factor(Acuity..),
    AnimalID = factor(Animal.ID)
  )


# Number of observations
n <- nrow(data)

# Adjust Best % looking
data$Best_looking_adj <- ((data$Best...Looking...acuity / 100) * (n - 1) + 0.5) / n

# Check for missing values
sum(is.na(data$Best_looking_adj))

# Model formula without interactions
formula_global <- bf(
  Best_looking_adj | mi() ~ Group +
    (1 | left_right_acuity) +
    (1 | acuity_number) +
    (1 | AnimalID)
)

# Define informative priors for the regression coefficients
priors <- c(
  prior(normal(0, 1), class = "b"),        # Fixed effects coefficients
  prior(exponential(1), class = "phi"),    # Precision parameter
  prior(normal(0, 1), class = "Intercept") # Intercept
)

# Refit the model with priors
model_global <- brm(
  formula = formula_global,
  data = data,
  family = Beta(),
  prior = priors,
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 123,
  control = list(adapt_delta = 0.99)
)

summary(model_global, prob = 0.89)

# Obtain estimated marginal means
emm_global <- emmeans(model_global, ~ Group, epred = TRUE)

# Perform pairwise contrasts between treatments 0 and 1
contrast_0_vs_1 <- contrast(emm_global, method = "pairwise", adjust = "none", level = 0.89)

# View the contrast
print(contrast_0_vs_1)


# Model formula including interactions
formula_interaction <- bf(
  Best_looking_adj | mi() ~ Group * Contrast_level * Age +
    (1 | left_right_acuity) +
    (1 | acuity_number) +
    (1 | AnimalID)
)


# Set weakly informative priors
priors <- c(
  prior(normal(0, 2.5), class = "b"),       # Coefficients
  prior(exponential(1), class = "phi")      # Precision parameter
)


# Fit the Beta regression model
model_interaction <- brm(
  formula = formula_interaction,
  data = data,
  family = Beta(),
  prior = priors,          
  chains = 4,
  cores = 4,
  iter = 4000,
  seed = 148,
  control = list(adapt_delta = 0.95)
)

# Summarize the model
summary(model_interaction, prob = 0.89)
coefs <- posterior_summary(model_interaction, probs = c(0.055, 0.945))
coefs <- round(coefs, 3)
write.csv(coefs, "vis_act_interaction_model.csv", row.names = TRUE)

# Plot trace plots
plot(model_interaction)


# Posterior predictive check
pp_check(model_interaction) + theme_minimal()


# Obtain estimated marginal means
emm_interaction <- emmeans(
  model_interaction,
  ~ Group | Contrast_level + Age,
  epred = TRUE
)


# Perform pairwise contrasts between treatments
contrasts_interaction <- contrast(emm_interaction, method = "pairwise", level = 0.89)

# Summarize contrasts
summary_contrasts <- summary(contrasts_interaction, infer = TRUE)
summary_contrasts <- summary_contrasts %>%
  mutate(across(where(is.numeric), ~ round(.x, 3)))
write_csv(summary_contrasts, "viz_act_contrasts_results.csv")

# View the contrasts
print(summary_contrasts)

# Extract posterior samples for contrasts
contrast_samples <- as.data.frame(contrasts_interaction)

# Extract estimated marginal means for treatments by age and contrast_level
emm_interaction <- emmeans(
  model_interaction,
  pairwise ~ Group | Age + Contrast_level,
  epred = TRUE
)

# Extract the posterior samples of the contrasts
contrast_samples <- emm_interaction$contrasts %>%
  gather_emmeans_draws() %>%
  ungroup()


# Ensure 'age' and 'contrast_level' are treated as factors with appropriate labels
contrast_samples <- contrast_samples %>%
  mutate(
    age = factor(Age, levels = unique(Age), labels = paste("Age", unique(Age))),
    contrast_level = factor(Contrast_level, levels = unique(Contrast_level), labels = paste("Contrast", unique(Contrast_level))),
    contrast = as.factor(contrast)
  )


# Create the plot
ggplot(contrast_samples, aes(x = Contrast_level, y = .value)) + 
  stat_halfeye(
    fill = "#8dd3c7",
    alpha = 0.8,
    .width = c(0.66, 0.89),
    slab_color = NA,
    interval_color = "black",
    position = position_dodge(width = 0.8)
    
  ) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "gray40") +
  facet_wrap(~ Age, ncol = 1) +
  scale_x_discrete(limits = c("High", "Medium", "Low")) + 
  labs(
    #title = "Visual Acuity: ZIKV - Control Contrasts by Age Group",
    y = "Estimated Difference in % Looking",
    x = "Contrast"
  ) +
  theme_classic(base_size = 16) +
  theme(
    legend.title = element_text(size = 14),
    legend.text = element_text(size = 12),
    axis.title = element_text(size = 16),
    axis.text = element_text(size = 12),
    plot.title = element_text(size = 18, face = "bold")
  )



###############################################################################
# SENSITIVITY ANALYSIS: CONTROL GROUPS ANALYSED SEPARATELY (UIC, PIC, Zika)
#
# The primary analysis above pools UIC + PIC into "Control" (data$Group,
# line 19). Because the acuity result is a NULL (no ZIKV-Control difference),
# this section tests if (a) the null holds when ZIKV is compared against
# EACH control separately (Zika - UIC, Zika - PIC) and (b) the two controls 
# differ (PIC - UIC).
#
# Both primary models use a blanket prior(normal(0, .), class = "b"),
# so moving from 2-level Group to 3-level Group3 needs NO prior change: the new
# Group3PIC coefficients are covered automatically. Contrasts use epred = TRUE
# (proportion scale, as in the primary analysis) and are summarized as
# median + 89% HDI.
#
# CAUTION: per-group cells are small; treat separate-control contrasts as a
# consistency check with wide intervals.
# NOTE: models are cached to brms_3group_va_*.rds; delete those to refit.
###############################################################################

message("\n\n=== VISUAL ACUITY SENSITIVITY: UIC, PIC, Zika kept separate ===")

# 3-level treatment factor (does NOT overwrite data$Group)
data$Group3 <- factor(
  data$Treatment.Group..0.UIC..1.PIC..2.Zika.,
  levels = c(0, 1, 2), labels = c("UIC", "PIC", "Zika")
)
cat("Group3 x Sex counts:\n"); print(table(data$Group3, data$Sex))

## --- Refit the two primary models with Group3 (priors unchanged) ---
priors_global3 <- c(
  prior(normal(0, 1), class = "b"),          # blanket: Group3PIC + Group3Zika
  prior(exponential(1), class = "phi"),
  prior(normal(0, 1), class = "Intercept")
)
model_global_3 <- brm(
  formula = bf(Best_looking_adj | mi() ~ Group3 +
                 (1 | left_right_acuity) + (1 | acuity_number) + (1 | AnimalID)),
  data = data, family = Beta(), prior = priors_global3,
  chains = 4, cores = 4, iter = 4000, seed = 123,
  control = list(adapt_delta = 0.99),
  file = "brms_3group_va_global"
)

priors_inter3 <- c(
  prior(normal(0, 2.5), class = "b"),        # blanket for all Group3 * Contrast * Age terms
  prior(exponential(1), class = "phi")
)
model_interaction_3 <- brm(
  formula = bf(Best_looking_adj | mi() ~ Group3 * Contrast_level * Age +
                 (1 | left_right_acuity) + (1 | acuity_number) + (1 | AnimalID)),
  data = data, family = Beta(), prior = priors_inter3,
  chains = 4, cores = 4, iter = 4000, seed = 148,
  control = list(adapt_delta = 0.95),
  file = "brms_3group_va_interaction"
)

## --- Extract UIC/PIC/Zika contrasts (epred/proportion scale, median + 89% HDI) ---
# revpairwise on (UIC, PIC, Zika) -> "PIC - UIC", "Zika - UIC", "Zika - PIC".
extract_va <- function(model, spec, model_label, ci = 0.89) {
  emm <- emmeans(model, spec, epred = TRUE)
  ctr <- contrast(emm, method = "revpairwise")
  dr  <- gather_emmeans_draws(ctr)
  grp <- intersect(c("contrast", "Contrast_level", "Age"), names(dr))
  dr %>%
    group_by(across(all_of(grp))) %>%
    median_hdci(.value, .width = ci) %>%
    ungroup() %>%
    mutate(model = model_label,
           credible = (.lower > 0) | (.upper < 0))
}

va_sens <- bind_rows(
  extract_va(model_global_3,      ~ Group3,                          "Global (overall)"),
  extract_va(model_interaction_3, ~ Group3 | Contrast_level + Age,   "Interaction (contrast x age)")
)
write.csv(va_sens, "va_sensitivity_3group_all_contrasts.csv", row.names = FALSE)
message("Saved va_sensitivity_3group_all_contrasts.csv")

## --- PIC vs UIC : pooling justification ---
va_pic_uic <- va_sens %>%
  filter(contrast == "PIC - UIC") %>%
  select(model, Contrast_level, Age, .value, .lower, .upper, credible) %>%
  arrange(model, Age, Contrast_level)
cat("\n=== PIC vs UIC (proportion scale; credible == TRUE flags a difference) ===\n")
print(as.data.frame(va_pic_uic))
write.csv(va_pic_uic, "va_control_PIC_vs_UIC.csv", row.names = FALSE)
cat(sprintf("\nPIC vs UIC: %d of %d contrasts credible.\n",
            sum(va_pic_uic$credible, na.rm = TRUE), nrow(va_pic_uic)))

## --- ZIKV vs each control : confirm the null holds against both ---
va_zika <- va_sens %>%
  filter(contrast %in% c("Zika - UIC", "Zika - PIC")) %>%
  select(model, Contrast_level, Age, contrast, .value, .lower, .upper, credible) %>%
  arrange(model, Age, Contrast_level, contrast)
cat("\n=== ZIKV vs UIC and ZIKV vs PIC (expect null / intervals spanning 0) ===\n")
print(as.data.frame(va_zika))
write.csv(va_zika, "va_sensitivity_zika_vs_each_control.csv", row.names = FALSE)
cat(sprintf("\nZika vs each control: %d of %d contrasts credible.\n",
            sum(va_zika$credible, na.rm = TRUE), nrow(va_zika)))

message("\n=== VISUAL ACUITY SENSITIVITY ANALYSIS COMPLETE ===")
