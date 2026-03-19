# Do Voice Identity Ratings of AI Cloned Voices Affect Their Intelligibility in Noise?

Division of Psychology and Language Sciences, University College London (UCL)

March 2026

## Project Summary
This project investigates whether voice identity ratings of AI cloned voices are associated with differences in intelligibility and perceived social traits under noisy listening conditions.

The analysis pipeline focuses on:
- deriving latent social-perception dimensions using PCA
- comparing human vs clone voices on these dimensions
- examining voice-specific effects
- describing acoustic differences between human and cloned speech
- producing publication figures

## Ethics
Ethical approval was granted by the UCL Research Ethics Committee.

Approval Number: `0599.005`

## Repository Structure
- `Code/01_generate_publication_plots.R`: main script to generate publication figures
- `Code/VoiceIdentity_PCA_analysis.r`: full analysis script (legacy/original workflow)
- `Code/assumption_checks_axis_only.R`: assumption-check script
- `Code/Figures/`: output folder for publication-ready PNG figures
- `Code/Outputs/`: output folder for generated tables
- `Data/VoiceIdentity_forAnalysis (1).csv`: primary ratings dataset
- `Data/Acoustics_short.csv`: acoustic feature dataset

## Figure Outputs
Running the plotting script generates the following files in `Code/Figures/`:
- `Figure_01_H1_PCA_Scree_and_Loadings.png`
- `Figure_02_H2_Clone_vs_Human_PCA_Axes.png`
- `Figure_03_H3_Valence_By_Voice.png`
- `Figure_04_H3_Dominance_By_Voice.png`
- `Figure_05_H3_Acoustic_Differences_By_Voice.png`

## How to Run
From the project root:

```bash
Rscript Code/01_generate_publication_plots.R
```

## R Requirements
The plotting script uses these R packages:
- `dplyr`
- `tidyr`
- `ggplot2`
- `grid`

Additional scripts in this repository may require:
- `lme4`
- `lmerTest`
- `psych`
- `FactoMineR`
- `corrplot`
- `broom`
- `purrr`
- `forcats`

## Reproducibility Notes
- Run analyses from the project root so relative paths resolve correctly.
- Keep files in `Data/` with their current names unless script paths are updated.
- For submission, include an R session record (`sessionInfo()`) and package versions.

## Suggested Citation / Identification Text
Do Voice Identity Ratings of AI Cloned Voices Affect Their Intelligibility in Noise?
Division of Psychology and Language Sciences, University College London.
March 2026.
Ethical approval: UCL Research Ethics Committee (`0599.005`).
