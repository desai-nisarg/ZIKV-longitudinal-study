## Altered neurodevelopment, CSF volume, and emotional reactivity following postnatal Zika virus infection in infant macaques

Manuscript authors: Nisarg Desai, Kaitlyn Love, Alex van Schoor, Sienna Freeman, Muskan Ali, Rebecca Richardson, Zsofia A Kovacs-Balint, Ruy Amaru Tobar Mosqueira, Rachel Lebovic, Jose Acevedo-Polo, Roza Vlasova, Martin Styner, Mar M Sanchez, Kathryn Moore, Nils Schoof, Patrick Whang, Vidisha Singh, Venkata-Viswanadh Edara, Mehul S. Suthar, Ann Chahroudi, Jessica Raper

Last update: November 5, 2026

bioRxiv [link](https://www.biorxiv.org/content/10.64898/2026.03.27.714817v1.abstract)

This repository contains the **R scripts** and input data to reproduce the analyses reported in the paper: viral-load, cytokine, neurobehavioral, visual-acuity, acute-stress, attachment, and neuroimaging analyses in infant rhesus macaques comparing **Zika-infected** and **control** conditions (with sex and, where relevant, other factors). Models are fit with **`brms`** (Stan) using weakly-informative or study-specific priors (synthesized from prior-study summaries where noted), and effects are summarized by the posterior median with an 89% highest-density interval (HDI) and the probability of direction (pd).

The primary analyses compare Zika-infected animals against a pooled control group (uninfected controls, UIC, and Poly-IC-treated controls, PIC). Where noted below, scripts also include sensitivity analyses that keep the three treatment groups separate (UIC, PIC, Zika), to check that the Zika effects hold against each control group and to compare the two control groups directly.

### Analysis scripts

Scripts are listed in the order the corresponding analyses appear in the manuscript.

| Script | Focus |
|--------|--------|
| `viral_load_AUC_analysis.R` | Viral-load burden (area under the viremia curve, AUC) within the Zika group: sex differences in burden (Gaussian model on log AUC), and associations of burden with the neuroimaging and behavioral outcomes (each outcome in its own model family). Every model is also refit excluding one extreme low-AUC animal as a sensitivity check. Reads the outcome datasets listed under Data. |
| `cytokine_analysis.R` | Innate-immune activation: cytokine/chemokine responses in Poly-IC versus uninfected controls (per-analyte Gaussian models on the log10 fold-change scale). |
| `neurobehavioral_assessment.R` | Infant Neurobehavioral Assessment Scale (INAS) ordinal outcomes (cumulative-logit models). Three-part workflow: (1) model fitting with age as a continuous predictor, (2) contrasts and plots, (3) sensitivity analysis treating age as a factor, with per-week contrasts including the weeks before infection. |
| `visual_acuity.R` | Visual acuity (beta regression with hierarchical structure), plus the separate-control sensitivity analysis. |
| `acute_stress_assessment.R` | Acute-stress (Human Intruder) behaviors (negative binomial models); weighted priors from `HI_prior.xlsx` plus `human_intruder.csv`. Includes the separate-control sensitivity analysis. |
| `attachment_assessment.R` | Attachment behaviors (negative binomial models for counts; beta regression for the index of preference); behavior-specific priors. Models threat barks and hostile behavior as separate outcomes, and includes the separate-control sensitivity analysis. |
| `neuroimaging_assessment.R` | Brain MRI volumetric outcomes; optional volume adjustment (ICV/TBV); Zika-vs-control contrasts by sex and the group × sex interaction; priors from `prior.csv`. Section 10 contains the separate-control sensitivity analysis, including regional UIC-vs-PIC comparisons. |
| `compute_interaction_pd.R` | Holistic Bayesian summaries (pd, 89% HDI/ETI, and P(\|d\| > 0.2)) for the neuroimaging group × sex interactions, computed from the saved regional model objects. |
| `individual_data_figures.R` | Per-animal data figures (box-and-whisker plots with individual points, distinguished by sex) across all assessments. |

### Data

Place the following input files in the working directory before running:

- `ZIKV-R01_LongitudinalStudy_animals&sex_8-14-26.xlsx` — animal metadata (publication order, cohort, sex, treatment group), viral-load burden (AUC; Zika-infected animals only), and cytokine/chemokine log10 fold-change values (IFN-α2a, IFN-γ, IL-1Ra, IP-10); used by the viral-load and cytokine analyses
- `SNIP_cleaned.csv`, `SNIP_prior_data.csv` — INAS scores and prior summaries
- `visual_acuity.csv` — visual-acuity looking data
- `human_intruder.csv`, `HI_prior.xlsx` — acute-stress behaviors and priors
- `attachment.csv` — attachment behaviors
- `ZIKV_3moMRI.csv` — 3-month MRI regional volumes
- `prior.csv` — neuroimaging priors

### Requirements

- R (analyses in the paper used R 4.4.0 and `brms` 2.21.0) with a working **Stan** toolchain for `brms`
- R packages: `brms`, `tidyverse` (`dplyr`, `tidyr`, `ggplot2`, `readr`, `stringr`, `purrr`, `forcats`, `tibble`), `emmeans`, `tidybayes`, `posterior`, `bayesplot`, `bayestestR`, `parameters`, `readxl`, `janitor`, `patchwork`, `cowplot`, `ggridges`, `pacman`. Install any that are missing with `install.packages("package-name")`.

**Note:** Model fitting can be slow. Some scripts save fitted model objects to disk and reload them on later runs instead of refitting. `compute_interaction_pd.R` loads models fit by `neuroimaging_assessment.R`, so run the model-fitting scripts before the scripts that summarize their outputs. Scripts write their result tables and figures to the working directory.
