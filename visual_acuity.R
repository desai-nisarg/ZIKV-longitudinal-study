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
    Group = factor(Treatment.Group..0.UIC..1.PIC..2.Zika., levels = c(0,1,2), labels = c("Control", "Control", "Zika")), 
    Contrast_level = factor(Contrast..1.low..2.medium..3.high., levels = c(1,2,3), labels = c("Low", "Medium", "High")),
    Age = factor(Age..months., levels = c(4,6), labels = c("4 Months", "6 Months")),
    left_right_acuity = factor(Left.or.Right.Acuity..0.left..1.right., levels = c(0,1), labels = c("Left", "Right")),    
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
