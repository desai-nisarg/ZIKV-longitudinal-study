################################################################################
# Supplementary individual-data figures
#
# One or more figures per assessment showing the underlying per-animal data as
# simple box-and-whisker plots with every individual animal overlaid as a point
# (colored/shaped by sex). Boxes are deliberately wide so the individual
# points are easy to see. All measured variables are shown for each assessment.
#
# To keep each figure legible, the larger assessments are split into panels of
# roughly 6-9 per figure (letter-suffixed: S3a-S3e, S4a-S4b, S5a-S5c).
#
# Outputs (PNG, 300 dpi) to the working directory:
#   FigS2_INAS_individual.png              (neurobehavioral trajectories, 4 domains)
#   FigS3a_acutestress_hostility.png       (Threat barks, All hostile)
#   FigS3b_acutestress_fear.png            (Fearful, Freezing)
#   FigS3c_acutestress_anxiety.png         (Anxious, Self-directed)
#   FigS3d_acutestress_affiliation.png     (Affiliative vocals, Screams)
#   FigS4a_attachment_affiliation.png      (Index, Screams, Attention-seeking, Affiliative, Fearful, Anxious)
#   FigS4b_attachment_hostility.png        (Hostile-defensive, Threat barks, Hostile, Self-directed, Displacement, Falls)
#   FigS5a_volumes_global.png              (Intracranial, Total brain, Total CSF, Lateral ventricle)
#   FigS5b_volumes_limbic.png              (Hippocampus, Amygdala, L/R temporal-limbic)
#   FigS5c_volumes_sensory.png             (L/R temporal-auditory, visual, occipital)
#   FigS6_cytokine_individual.png          (all analytes, UIC vs PIC)
#   FigS7_visualacuity_individual.png      (all ages x contrasts)
#
# Run from the folder that contains the data files.
################################################################################

library(tidyverse)
library(readxl)

sex_cols   <- c(Female = "#fb8072", Male = "#80b1d3")
grp_levels <- c("UIC", "PIC", "ZIKV")

prep <- function(d, group_col, sex_col) {
  d %>%
    mutate(Group = factor(.data[[group_col]], levels = c(0,1,2), labels = grp_levels),
           Sex   = factor(.data[[sex_col]],   levels = c(0,1),   labels = c("Female","Male")))
}

# ---- simple wide box-and-whisker + individual points (facet_wrap version) ----
boxfig <- function(df, outcomes, title, ylab = "Value", ncol = 3,
                   box_width = 0.6, jitter_w = 0.14, pt_size = 2.3) {
  d <- df %>% filter(outcome %in% outcomes) %>%
    mutate(outcome = factor(outcome, levels = outcomes))
  ggplot(d, aes(x = Group, y = value)) +
    geom_boxplot(width = box_width, outlier.shape = NA, alpha = 0.35,
                 fill = "grey85", colour = "grey30", linewidth = 0.4, na.rm = TRUE) +
    geom_point(aes(colour = Sex, shape = Sex),
               position = position_jitter(width = jitter_w, height = 0, seed = 1),
               size = pt_size, alpha = 0.9, stroke = 0.4, na.rm = TRUE) +
    scale_colour_manual(values = sex_cols) +
    facet_wrap(~ outcome, scales = "free_y", ncol = ncol) +
    labs(title = title, x = "Treatment group", y = ylab) +
    theme_bw(base_size = 13) +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold", size = 10),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
}

## --- Figure S2: Neurobehavioral (INAS) individual trajectories ---------------
## Longitudinal outcome: kept as per-animal trajectories (all 4 domains, every
## animal) because a box-and-whisker cannot show the course over age.
inas <- read_csv("SNIP_cleaned.csv", show_col_types = FALSE) %>%
  prep("Treatment_Group_0UIC_1PIC_2Zika", "Sex_0f_1m")
inas_long <- inas %>%
  pivot_longer(c(Orientation_Scores_02, Neuromotor_Scores_02,
                 Motor_Scores_02, Temperament_Scores_02),
               names_to = "domain", values_to = "score") %>%
  mutate(domain = recode(domain,
           Orientation_Scores_02 = "Orientation", Neuromotor_Scores_02 = "Neuromotor",
           Motor_Scores_02 = "Motor", Temperament_Scores_02 = "State Control")) %>%
  group_by(Animal_ID, Age_Weeks, Group, Sex, domain) %>%
  summarise(score = mean(score, na.rm = TRUE), .groups = "drop")

pS2 <- ggplot(inas_long, aes(Age_Weeks, score, colour = Group)) +
  geom_vline(xintercept = 5, linetype = "dotted", colour = "grey40") +
  geom_line(aes(group = Animal_ID), alpha = 0.35, linewidth = 0.5) +
  stat_summary(aes(group = Group), fun = mean, geom = "line", linewidth = 1.3) +
  facet_wrap(~ domain, scales = "free_y") +
  scale_colour_brewer(palette = "Dark2") +
  labs(title = "Neurobehavioral (INAS): individual trajectories",
       x = "Age (weeks)", y = "Mean item score",
       subtitle = "Thin lines = individual animals; bold lines = group means; dotted line = challenge (~5 weeks)") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("FigS2_INAS_individual.png", pS2, width = 10, height = 7, dpi = 300)

## --- Acute stress (Human Intruder): behavior x condition, split into groups --
hi <- read_csv("human_intruder.csv", show_col_types = FALSE) %>%
  prep("Treatment Group (0=UIC, 1=PIC, 2=Zika)", "Sex (0=f, 1=m)") %>%
  mutate(Condition = factor(`Condition (0=alone, 1=profile, 2=stare)`,
                            levels = c(0,1,2), labels = c("Alone","Profile","Stare")))
# The 8 behaviors analyzed for acute stress. Displacement, threats-toward, and
# threats-away are NOT shown separately: they are subsumed within All Hostile
# Behavior and were not modeled as separate outcomes.
hi_behaviours <- c(
  `Affiliative Vocals (Freq)`         = "Affiliative vocals",
  `Screams (Freq)`                    = "Screams",
  `Threat Barks (Freq)`               = "Threat barks",
  `Anxious Behavior (Freq)`           = "Anxious",
  `Self-Directed Behavior (Dur, sec)` = "Self-directed (s)",
  `Fearful Behavior (Freq)`           = "Fearful",
  `All Hostile Behavior (Freq)`       = "All hostile",
  `Freezing (Dur, sec)`               = "Freezing (s)")
hi_long <- hi %>%
  pivot_longer(all_of(names(hi_behaviours)), names_to = "outcome", values_to = "value") %>%
  mutate(outcome = recode(outcome, !!!hi_behaviours))

hi_fig <- function(behavs, file, title, w, h) {
  d <- hi_long %>% filter(outcome %in% behavs) %>%
    mutate(outcome = factor(outcome, levels = behavs))
  p <- ggplot(d, aes(Group, value)) +
    geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.35,
                 fill = "grey85", colour = "grey30", linewidth = 0.4, na.rm = TRUE) +
    geom_point(aes(colour = Sex, shape = Sex),
               position = position_jitter(width = 0.14, height = 0, seed = 1),
               size = 1.9, alpha = 0.9, stroke = 0.4, na.rm = TRUE) +
    scale_colour_manual(values = sex_cols) +
    facet_grid(outcome ~ Condition, scales = "free_y") +
    labs(title = title, x = "Treatment group", y = "Frequency (or duration, s)") +
    theme_bw(base_size = 12) +
    theme(legend.position = "top",
          strip.text = element_text(face = "bold", size = 10),
          plot.title = element_text(face = "bold", hjust = 0.5, size = 15))
  ggsave(file, p, width = w, height = h, dpi = 300)
}
hi_fig(c("Threat barks","All hostile"),
       "FigS3a_acutestress_hostility.png",
       "Acute stress: hostility-related behaviors", w = 9, h = 6.5)  # 6 panels
hi_fig(c("Fearful","Freezing (s)"),
       "FigS3b_acutestress_fear.png",
       "Acute stress: fear and freezing", w = 9, h = 6.5)
hi_fig(c("Anxious","Self-directed (s)"),
       "FigS3c_acutestress_anxiety.png",
       "Acute stress: anxiety-related behaviors", w = 9, h = 6.5)
hi_fig(c("Affiliative vocals","Screams"),
       "FigS3d_acutestress_affiliation.png",
       "Acute stress: affiliation and distress", w = 9, h = 6.5)

## --- Attachment: all behaviors + index, split into two figures of 6 ---------
att <- read_csv("attachment.csv", show_col_types = FALSE) %>%
  prep("Treatment Group (0=UIC, 1=PIC, 2=Zika)", "Sex (0=f, 1=m)")
att_outcomes <- c(
  `Index of Preference (Scale, -1 to 1)` = "Index of preference",
  `Scream (Freq)`                        = "Screams",
  `Attention Seeking Behaviors (Freq)`   = "Attention-seeking",
  `Affiliative Behaviors (Freq)`         = "Affiliative",
  `Fearful Behaviors (Freq)`             = "Fearful",
  `Hostile Defensive Behaviors (Freq)`   = "Hostile-defensive",
  `Anxious Behaviors (Freq)`             = "Anxious",
  `Self-Directed Behaviors (Dur, sec)`   = "Self-directed (s)",
  `Displacement Behaviors (Freq)`        = "Displacement",
  `Falls (Freq)`                         = "Falls",
  `Treat Barks (Freq)`                   = "Threat barks",
  `Hostile Behavior (Freq)`              = "Hostile behaviour")
att_long <- att %>%
  pivot_longer(all_of(names(att_outcomes)), names_to = "outcome", values_to = "value") %>%
  mutate(outcome = recode(outcome, !!!att_outcomes))
pS4a <- boxfig(att_long,
  c("Index of preference","Screams","Attention-seeking","Affiliative","Fearful","Anxious"),
  "Attachment: affiliation, distress, and preference",
  ylab = "Frequency (or duration, s, or index)", ncol = 3)
ggsave("FigS4a_attachment_affiliation.png", pS4a, width = 11, height = 7, dpi = 300)
pS4b <- boxfig(att_long,
  c("Hostile-defensive","Threat barks","Hostile behaviour","Self-directed (s)","Displacement","Falls"),
  "Attachment: hostility and agitation",
  ylab = "Frequency (or duration, s)", ncol = 3)
ggsave("FigS4b_attachment_hostility.png", pS4b, width = 11, height = 7, dpi = 300)

## --- Neuroimaging: all volumes, split into global / limbic / sensory ---------
mri <- read_csv("ZIKV_3moMRI.csv", show_col_types = FALSE) %>%
  prep("Treatment Group (0=UIC, 1=PIC, 2=Zika)", "Sex (0=f, 1=m)")
mri_outcomes <- c(
  `Total Hippocampus (mm3)`     = "Hippocampus",
  `Total Amygdala (mm3)`        = "Amygdala",
  `Total Intracranial Vol (mm3)`= "Intracranial",
  `Total Brain Vol (mm3)`       = "Total brain",
  `Total CSF volume (mm3)`      = "Total CSF",
  `Lateral Ventricle (mm3)`     = "Lateral ventricle",
  `RT Temporal Auditory (GM+WM)`= "R temporal-auditory",
  `RT Temporal Visual  (GM+WM)` = "R temporal-visual",
  `RT Temporal Limbic (GM+WM)`  = "R temporal-limbic",
  `RT Occipital (GM+WM)`        = "R occipital",
  `LT Temporal Auditory (GM+WM)`= "L temporal-auditory",
  `LT Temporal Visual  (GM+WM)` = "L temporal-visual",
  `LT Temporal Limbic (GM+WM)`  = "L temporal-limbic",
  `LT Occipital (GM+WM)`        = "L occipital")
mri_long <- mri %>%
  pivot_longer(all_of(names(mri_outcomes)), names_to = "outcome", values_to = "value") %>%
  mutate(outcome = recode(outcome, !!!mri_outcomes))
vol.lab <- expression(Volume~(mm^3))
pS5a <- boxfig(mri_long,
  c("Intracranial","Total brain","Total CSF","Lateral ventricle"),
  "Global brain and CSF volumes", ylab = vol.lab, ncol = 2)
ggsave("FigS5a_volumes_global.png", pS5a, width = 8, height = 7, dpi = 300)
pS5b <- boxfig(mri_long,
  c("Hippocampus","Amygdala","R temporal-limbic","L temporal-limbic"),
  "Subcortical and limbic volumes", ylab = vol.lab, ncol = 2)
ggsave("FigS5b_volumes_limbic.png", pS5b, width = 8, height = 7, dpi = 300)
pS5c <- boxfig(mri_long,
  c("R temporal-auditory","L temporal-auditory","R temporal-visual",
    "L temporal-visual","R occipital","L occipital"),
  "Sensory cortical volumes (auditory, visual, occipital)", ylab = vol.lab, ncol = 3)
ggsave("FigS5c_volumes_sensory.png", pS5c, width = 11, height = 7, dpi = 300)

## --- Figure S6: Cytokines (UIC vs PIC; ZIKV not measured) ---------------------
cyt <- read_excel("ZIKV-R01_LongitudinalStudy_animals&sex_8-14-26.xlsx",
                  sheet = "All_subjects")
names(cyt) <- c("pub","cohort","sex","group","vrna","IFNa2a","IFNg","IL1Ra","IP10")
cyt <- cyt %>%
  filter(group %in% c(0,1)) %>%
  mutate(Group = factor(group, levels = c(0,1), labels = c("UIC","PIC")),
         Sex   = factor(sex, levels = c(0,1), labels = c("Female","Male")))
cyt_long <- cyt %>%
  pivot_longer(c(IFNa2a, IFNg, IL1Ra, IP10), names_to = "outcome", values_to = "value") %>%
  mutate(outcome = recode(outcome, IFNa2a = "IFN-alpha2a", IFNg = "IFN-gamma",
                          IL1Ra = "IL-1Ra", IP10 = "IP-10"))
pS6 <- boxfig(cyt_long, c("IFN-alpha2a","IFN-gamma","IL-1Ra","IP-10"),
              "Cytokine responses (Poly-IC vs uninfected control)",
              ylab = "log10 fold change", ncol = 2)
ggsave("FigS6_cytokine_individual.png", pS6, width = 8, height = 6, dpi = 300)

## --- Figure S7: Visual acuity -- ALL ages x contrasts ------------------------
va <- read_csv("visual_acuity.csv", show_col_types = FALSE) %>%
  prep("Treatment Group (0=UIC, 1=PIC, 2=Zika)", "Sex (0=f, 1=m)") %>%
  mutate(Age = paste0(`Age (months)`, " mo"),
         Contrast = factor(`Contrast (1=low, 2=medium, 3=high)`,
                           levels = c(1,2,3), labels = c("Low","Medium","High")),
         value = `Best % Looking @ acuity`)
pS7 <- ggplot(va, aes(Group, value)) +
  geom_boxplot(width = 0.6, outlier.shape = NA, alpha = 0.35,
               fill = "grey85", colour = "grey30", linewidth = 0.4, na.rm = TRUE) +
  geom_point(aes(colour = Sex, shape = Sex),
             position = position_jitter(width = 0.14, height = 0, seed = 1),
             size = 1.9, alpha = 0.9, stroke = 0.4, na.rm = TRUE) +
  scale_colour_manual(values = sex_cols) +
  facet_grid(Age ~ Contrast) +
  labs(title = "Visual acuity: individual data",
       x = "Treatment group", y = "Percent looking") +
  theme_bw(base_size = 13) +
  theme(legend.position = "top", strip.text = element_text(face = "bold"),
        plot.title = element_text(face = "bold", hjust = 0.5))
ggsave("FigS7_visualacuity_individual.png", pS7, width = 10, height = 7, dpi = 300)

cat("--- INDIVIDUAL-DATA FIGURES WRITTEN: FigS2, S3a-e, S4a-b, S5a-c, S6, S7 ---\n")
