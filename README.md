## Increased CSF volume, altered brain development and emotional reactivity after postnatal Zika virus infection in infant rhesus macaques

Manuscript authors: Nisarg Desai, Kaitlyn Love, Alex van Schoor, Sienna Freeman, Muskan Ali, Rebecca Richardson, Zsofia A Kovacs-Balint, Ruy Amaru Tobar Mosqueira, Rachel Lebovic, Jose Acevedo-Polo, Roza Vlasova, Martin Styner, Mar M Sanchez, Kathryn Moore, Nils Schoof, Patrick Whang, Vidisha Singh, Venkata-Viswanadh Edara, Mehul S. Suthar, Ann Chahroudi, Jessica Raper

Code author: Nisarg P. Desai, nisarg.parimal.desai[at]emory[dot]edu

Last update: April 10, 2026

bioRxiv [link](https://www.biorxiv.org/content/10.64898/2026.03.27.714817v1.abstract)

Bayesian analyses in R (`brms`) for postnatal Zika-related neuroimaging and neurobehavioral outcomes, using informative priors from prior literature where noted.
This repository contains **R scripts** for **Bayesian models** fit with **`brms`** (Stan) on experimental data comparing **control and Zika** conditions (with **sex** and, where relevant, other factors). Analyses use **study-specific priors** synthesized from control-like or prior-study summaries when those inputs are provided.
| Script | Focus |
|--------|--------|
| `neuroimaging_assessment.R` | Brain MRI outcomes (ZIKV cohort); optional volume adjustment (e.g. ICV/TBV); priors from `prior.csv`. |
| `neurobehavioral_assessment.R` | **SNIP**-style ordinal outcomes; two-part workflow (model fitting, then contrasts/plots). |
| `acute_stress_assessment.R` | **Human intruder** / acute stress–related behaviors; weighted priors from `HI_prior.xlsx` plus `human_intruder.csv`. |
| `attachment_assessment.R` | **Attachment**-related behaviors (counts/durations); behavior-specific priors. |
| `visual_acuity.R` | **Visual acuity** (beta regression with hierarchical structure). |
## Requirements
- R and a working **Stan** toolchain for `brms`
- Packages each script loads (e.g. `tidyverse`, `readr`, `emmeans`, `tidybayes`; `readxl` for the human-intruder workflow)
Place the expected **CSV** files (and `HI_prior.xlsx` where needed) in the working directory before running.
**Note:** Model fitting can be slow; some scripts are split so you can save/load fitted objects between sessions.
