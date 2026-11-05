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

# Regularizing priors for the remaining count behaviors (revision).
# These previously used brms flat defaults, which are unstable for the
# low-frequency behaviors (near-separation with small cells). Here we apply
# the SAME weakly-informative group/sex priors already used for Screams and
# the Index of Preference -- normal(0,1) on class="b" (covers GroupZika,
# SexMale, and their interaction) -- plus a weakly-informative intercept and
# shape. This matches the acute-stress no-prior models as well.
prior_behaviors <- c(
  set_prior("normal(0, 5)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "b"),
  set_prior("exponential(1)", class = "shape")
)

for (behavior in behavior_vars) {
  formula <- as.formula(
    paste0(behavior, " ~ Group * Sex")
  )
  
  models[[behavior]] <- brm(
    formula = formula,
    data = data,
    prior = prior_behaviors,
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
write.csv(all_contrasts, "attachment_all_contrasts.csv")


###############################################################################
# ADD-ON: Threat Barks and Hostile as SEPARATE outcomes
###############################################################################
library(brms); library(emmeans); library(tidybayes); library(dplyr)

# 1. Pull the two new columns into the loaded `data` (the .RData `data` predates them).
#    read.csv renames "Treat Barks (Freq)" -> Treat.Barks..Freq. and
#    "Hostile Behavior (Freq)" -> Hostile.Behavior..Freq. 
.nd <- read.csv("attachment.csv")
data$ThreatBarks  <- .nd$Treat.Barks..Freq.
data$HostileAlone <- .nd$Hostile.Behavior..Freq.

# Re-assert the factors exactly as the script defines them (idempotent).
data$Group  <- factor(data$Treatment.Group..0.UIC..1.PIC..2.Zika., levels=c(0,1,2), labels=c("Control","Control","Zika"))
data$Group3 <- factor(data$Treatment.Group..0.UIC..1.PIC..2.Zika., levels=c(0,1,2), labels=c("UIC","PIC","Zika"))
data$Sex    <- factor(data$Sex..0.f..1.m.,                         levels=c(0,1),   labels=c("Female","Male"))

# 2. Same weakly-informative scheme used for the other (regularized) count behaviors.
prior_new <- c(
  set_prior("normal(0, 5)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "b"),
  set_prior("exponential(1)", class = "shape")
)

# 3. Primary (pooled Control vs Zika) models -- these feed the plot.
model_freq_threatbark <- brm(ThreatBarks  ~ Group * Sex, data=data, family=negbinomial(),
                             prior=prior_new, iter=4000, warmup=1000, chains=4, cores=4,
                             control=list(adapt_delta=0.95), seed=111, file="brms_attach_threatbark")
model_freq_hostile    <- brm(HostileAlone ~ Group * Sex, data=data, family=negbinomial(),
                             prior=prior_new, iter=4000, warmup=1000, chains=4, cores=4,
                             control=list(adapt_delta=0.95), seed=111, file="brms_attach_hostile")

# 4. Three-group (UIC, PIC, Zika) models -- for the vs-each-control sensitivity.
model3_threatbark <- brm(ThreatBarks  ~ Group3 * Sex, data=data, family=negbinomial(),
                         prior=prior_new, iter=4000, warmup=1000, chains=4, cores=4,
                         control=list(adapt_delta=0.95), seed=111, file="brms_3group_attach_reg_ThreatBarks")
model3_hostile    <- brm(HostileAlone ~ Group3 * Sex, data=data, family=negbinomial(),
                         prior=prior_new, iter=4000, warmup=1000, chains=4, cores=4,
                         control=list(adapt_delta=0.95), seed=111, file="brms_3group_attach_reg_HostileAlone")

# 5. Contrast extractor: by Sex + overall, median + 89% HDI + pd 
extract_new <- function(model, label, group_var) {
  do.call(bind_rows, lapply(c("By sex","Overall"), function(st) {
    sp  <- if (st == "By sex") as.formula(paste("~", group_var, "| Sex")) else as.formula(paste("~", group_var))
    dr  <- emmeans(model, sp) %>% contrast(method = "revpairwise") %>% gather_emmeans_draws()
    grp <- intersect(c("contrast","Sex"), names(dr))
    hdi <- dr %>% group_by(across(all_of(grp))) %>% median_hdci(.value, .width = 0.89) %>% ungroup()
    pd  <- dr %>% group_by(across(all_of(grp))) %>%
      summarise(pd = max(mean(.value > 0), mean(.value < 0)), .groups = "drop")
    left_join(hdi, pd, by = grp) %>%
      mutate(behavior = label, stratum = st, exp_effect = exp(.value),
             credible = (.lower > 0) | (.upper < 0))
  }))
}

# Primary (Zika - Control) and separate-control (PIC-UIC, Zika-UIC, Zika-PIC).
attach_new_primary <- bind_rows(
  extract_new(model_freq_threatbark, "Threat Barks",     "Group"),
  extract_new(model_freq_hostile,    "Hostile Behavior", "Group"))
attach_new_sens <- bind_rows(
  extract_new(model3_threatbark, "Threat Barks",     "Group3"),
  extract_new(model3_hostile,    "Hostile Behavior", "Group3"))

write.csv(attach_new_primary, "attach_NEWVARS_primary_contrasts.csv", row.names = FALSE)
write.csv(attach_new_sens,    "attach_NEWVARS_3group_contrasts.csv",  row.names = FALSE)
print(attach_new_sens)







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
  "Hostile" = model_freq_hostile,
  "Threat Barks" = model_freq_threatbark,
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

# --- Alternative desired_order vector ---

desired_order_with_anxious <- c(
  "Index of Preference",
  "Attention Seeking",  
  "Affiliative",
  "Screams",
  "Hostile",
  "Threat Barks",
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



###############################################################################
# SENSITIVITY ANALYSIS: CONTROL GROUPS ANALYSED SEPARATELY (UIC, PIC, Zika)
#
#
# The primary analysis above pools UIC + PIC into "Control" (data$Group,
# line 27). Here we refit every attachment-assessment model with the three
# groups kept SEPARATE (data$Group3: UIC, PIC, Zika) so we can (a) confirm the
# ZIKV effects hold against EACH control (Zika - UIC and Zika - PIC) and
# (b) show whether the two controls differ (PIC - UIC), which justifies pooling.
#
# Faithful to the primary analysis, model by model:
#   * Scream..Freq. : informative intercept prior (control mean). Every Group3
#     term gets the same normal(0,1) used for GroupZika, applied as a blanket
#     class = "b" prior so the new Group3PIC coefficients are covered.
#   * The count behaviors in cols 7:14 : brms DEFAULT priors, exactly as the
#     primary loop (no prior argument).
#   * Index of Preference : Beta(logit) with informative intercept prior; Group3
#     terms via blanket normal(0,1).
#
# Contrasts are stratified by Sex (the stratification used throughout this
# script) plus an overall marginal, summarized as median + 89% HDI. NOTE on
# exp(estimate): it is a RATE ratio for the negbinomial behaviors but an ODDS
# ratio for the Beta Index model -- see the 'family' column. Behaviors are
# keyed by their raw column name (no positional relabelling).
#
# CAUTION: per-sex control subgroups are very small (n ~ 3), so the
# separate-control contrasts have wide intervals and are a consistency check.
#
# NOTE: models are cached to brms_3group_attach_*.rds; delete those to refit.
###############################################################################

message("\n\n=== ATTACHMENT SENSITIVITY: UIC, PIC, Zika kept separate ===")

# 3-level treatment factor (does NOT overwrite data$Group)
data$Group3 <- factor(
  data$Treatment.Group..0.UIC..1.PIC..2.Zika.,
  levels = c(0, 1, 2), labels = c("UIC", "PIC", "Zika")
)
cat("Group3 x Sex counts:\n"); print(table(data$Group3, data$Sex))

models3       <- list()
model_family  <- c()   # track family per behaviour for correct exp() reading

## --- Scream (informative intercept + blanket b prior) ---
prior_scream3 <- c(
  set_prior(sprintf("normal(%f, %f)", mu_log_scream, se_log_scream), class = "Intercept"),
  set_prior("normal(0, 1)", class = "b"),          # covers Group3PIC, Group3Zika + interactions
  set_prior("exponential(1)", class = "shape")
)
models3[["Scream..Freq."]] <- brm(
  formula = Scream..Freq. ~ Group3 * Sex, data = data, family = negbinomial(),
  prior = prior_scream3, iter = 4000, warmup = 1000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.95), seed = 111,
  file = "brms_3group_attach_Scream"
)
model_family["Scream..Freq."] <- "negbinomial"

## --- Count behaviours (cols 7:14): regularized priors ---
## Same weakly-informative scheme as the (now regularized) primary loop:
## normal(0,5) intercept, normal(0,1) on all group/sex terms (covers
## Group3PIC, Group3Zika + interactions), exponential(1) shape. This tames
## the separation that flat defaults produce when the small control group is
## split three ways (e.g. Fearful, Attention-Seeking, Displacement).
## Cache filenames carry a 'reg' tag so brms refits rather than reloading the
## earlier default-prior fits.
prior_behaviors3 <- c(
  set_prior("normal(0, 5)", class = "Intercept"),
  set_prior("normal(0, 1)", class = "b"),
  set_prior("exponential(1)", class = "shape")
)
for (behavior in behavior_vars) {
  f <- as.formula(paste0(behavior, " ~ Group3 * Sex"))
  models3[[behavior]] <- brm(
    formula = f, data = data, family = negbinomial(),
    prior = prior_behaviors3,
    iter = 4000, warmup = 1000, chains = 4, cores = 4,
    control = list(adapt_delta = 0.95), seed = 111,
    file = paste0("brms_3group_attach_reg_", make.names(behavior))
  )
  model_family[behavior] <- "negbinomial"
}

## --- Index of Preference (Beta logit, informative intercept + blanket b) ---
prior_index3 <- c(
  set_prior(sprintf("normal(%f, %f)", pref_mu_logit, pref_se_logit), class = "Intercept"),
  set_prior("normal(0, 1)", class = "b")
)
models3[["Index of Preference"]] <- brm(
  formula = pref_adj ~ Group3 * Sex, data = data, family = Beta(link = "logit"),
  prior = prior_index3, iter = 4000, warmup = 1000, chains = 4, cores = 4,
  control = list(adapt_delta = 0.95), seed = 111,
  file = "brms_3group_attach_Index"
)
model_family["Index of Preference"] <- "Beta"

## --- Extract UIC/PIC/Zika contrasts (median + 89% HDI): by Sex + overall ---
# revpairwise on (UIC, PIC, Zika) -> "PIC - UIC", "Zika - UIC", "Zika - PIC".
extract_attach <- function(model, label, fam, ci = 0.89) {
  out <- list()
  for (stratum in c("By sex", "Overall")) {
    sp  <- if (stratum == "By sex") (~ Group3 | Sex) else (~ Group3)
    emm <- emmeans(model, sp)
    ctr <- contrast(emm, method = "revpairwise")
    dr  <- gather_emmeans_draws(ctr)
    grp <- intersect(c("contrast", "Sex"), names(dr))
    out[[stratum]] <- dr %>%
      group_by(across(all_of(grp))) %>%
      median_hdci(.value, .width = ci) %>%
      ungroup() %>%
      mutate(behavior = label, family = fam, stratum = stratum,
             exp_effect = exp(.value),
             credible   = (.lower > 0) | (.upper < 0))
  }
  bind_rows(out)
}

attach_sens <- bind_rows(lapply(names(models3), function(nm) {
  message(sprintf("Extracting 3-group contrasts for '%s'...", nm))
  extract_attach(models3[[nm]], nm, model_family[[nm]])
}))

write.csv(attach_sens, "attach_sensitivity_3group_all_contrasts.csv", row.names = FALSE)
message("Saved attach_sensitivity_3group_all_contrasts.csv")

## --- PIC vs UIC : pooling justification (mostly non-credible) ---
attach_pic_uic <- attach_sens %>%
  filter(contrast == "PIC - UIC") %>%
  select(behavior, family, stratum, Sex, .value, .lower, .upper, exp_effect, credible) %>%
  arrange(behavior, stratum)
cat("\n=== PIC vs UIC (credible == TRUE flags a control-group difference) ===\n")
print(as.data.frame(attach_pic_uic))
write.csv(attach_pic_uic, "attach_control_PIC_vs_UIC.csv", row.names = FALSE)
cat(sprintf("\nPIC vs UIC: %d of %d contrasts credible (89%% HDI excludes 0).\n",
            sum(attach_pic_uic$credible, na.rm = TRUE), nrow(attach_pic_uic)))

## --- ZIKV vs each control, by sex (for the sex-specific attachment findings) ---
attach_zika <- attach_sens %>%
  filter(stratum == "By sex", contrast %in% c("Zika - UIC", "Zika - PIC")) %>%
  select(behavior, family, Sex, contrast, .value, .lower, .upper, exp_effect, credible) %>%
  arrange(behavior, Sex, contrast)
cat("\n=== ZIKV vs UIC and ZIKV vs PIC, by sex (all behaviours) ===\n")
print(as.data.frame(attach_zika))
write.csv(attach_zika, "attach_sensitivity_zika_vs_each_control.csv", row.names = FALSE)

message("\n=== ATTACHMENT SENSITIVITY ANALYSIS COMPLETE ===")