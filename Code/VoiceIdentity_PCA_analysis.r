#######
# Statistical processign data Chang, voice identity processing experiment 
# Author: Patti Adank
# Date 12 Feb 2026
#######

#######################################################
# Preparatory steps 
######################################################

# Packages

library(dplyr)          # data wrangling
library(tidyr)          # data wrangling
library(FactoMineR)     # PCA
library(psych)          # rotation
library(corrplot)       # visualization
library(ggplot2)        # plotting

library(broom)
library(purrr)

library(lme4)
library(lmerTest)
library(forcats)

# Project paths
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- "--file="
  idx <- grep(file_arg, args)
  if (length(idx) > 0) {
    return(dirname(normalizePath(sub(file_arg, "", args[idx[1]]))))
  }
  normalizePath(getwd())
}

script_dir <- get_script_dir()
project_dir <- normalizePath(file.path(script_dir, ".."))
fig_dir <- file.path(project_dir, "Code", "Figures")
out_dir <- file.path(project_dir, "Code", "Outputs")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# Load in data file
Dat <- read.csv(file.path(project_dir, "Data", "VoiceIdentity_forAnalysis (1).csv"))

Dat <- Dat |>
  dplyr::rename(Stimulus = Spreadsheet..Audio) |>
  dplyr::mutate(Stimulus = gsub("\\.mp3$", "", Stimulus))

# restructure Dat so it is in the correct format for PCA 

voice_identity <- Dat |>
  tidyr::pivot_wider(
    id_cols = c(Participant.Public.ID,
                Stimulus,
                Voice,
                Type),
    names_from = Trait,
    values_from = Response
  )

trait_cols <- c(  "Trustworthiness",  "Dominance",  "Attractiveness",  "Competence",  "Warmth",  "Aggression", "Confidence", "Friendliness")

# PCA

## RUN PCA:
pca <- prcomp(voice_identity[trait_cols], scale = TRUE)

### principal components:
pca$rotation

# Extract PCs
pc_scores <- as.data.frame(pca$x)

# Keep scores for PC1 (valence) and PC2 (dominance) only. I checked the results and it looked
# like PC1 only had positive contributions from mostly positive/friendly traits
# and PC2 had Dominance loading highest

pc12 <- pc_scores[, c("PC1", "PC2")]

# Add them to your originial data frame and rename them and place them at the end for tinyness
voice_identity$Valence <- pc12$PC1
voice_identity$Dominance <- pc12$PC2
voice_identity <- voice_identity |>
  relocate(Valence, Dominance, .after = Friendliness)

# AI SECTION:
# Simple hypothesis plots (H1 and H2), saved to Code/Figures

# H1: Review-focused PCA figure with two panels (scree + loading map)
scree_data <- data.frame(
  PC = seq_along(pca$sdev),
  Variance = (pca$sdev^2) / sum(pca$sdev^2)
) %>%
  mutate(
    Cumulative = cumsum(Variance),
    Highlight = ifelse(PC %in% c(1, 2), "PC1_PC2", "Other"),
    VarianceLabel = ifelse(PC %in% c(1, 2), paste0("PC", PC, ": ", sprintf("%.1f%%", 100 * Variance)), NA_character_)
  )

h1_scree <- ggplot(scree_data, aes(x = PC, y = Variance, fill = Highlight)) +
  geom_col(width = 0.65) +
  geom_line(aes(y = Cumulative), colour = "#1F5AA6", linewidth = 0.9, group = 1) +
  geom_point(aes(y = Cumulative), colour = "#1F5AA6", size = 1.8) +
  geom_text(
    data = subset(scree_data, PC %in% c(1, 2)),
    aes(label = VarianceLabel),
    vjust = -0.6,
    size = 3.1,
    colour = "#1F5AA6",
    fontface = "bold",
    show.legend = FALSE
  ) +
  scale_x_continuous(breaks = scree_data$PC) +
  scale_fill_manual(values = c("PC1_PC2" = "#1F5AA6", "Other" = "#7A7A7A"), guide = "none") +
  expand_limits(y = max(scree_data$Variance) + 0.06) +
  labs(
    title = "H1 Panel A: Scree plot",
    x = "Principal component",
    y = "Proportion of variance"
  ) +
  theme_minimal(base_size = 11)

h1_load_data <- as.data.frame(pca$rotation[, c("PC1", "PC2")])
h1_load_data$Trait <- rownames(h1_load_data)

h1_loadings <- ggplot(h1_load_data, aes(x = PC1, y = PC2, label = Trait)) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey70", linetype = 2) +
  geom_point(colour = "#1F5AA6", size = 2.2) +
  geom_text(nudge_y = 0.03, size = 3.1) +
  coord_cartesian(
    xlim = range(h1_load_data$PC1) + c(-0.08, 0.08),
    ylim = range(h1_load_data$PC2) + c(-0.08, 0.08),
    clip = "off"
  ) +
  labs(
    title = "H1 Panel B: Trait loadings on PC1 and PC2",
    x = "PC1 loading",
    y = "PC2 loading"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.margin = margin(8, 20, 8, 20))

png(file.path(fig_dir, "figure_h1_pca_variance.png"), width = 3000, height = 1400, res = 300)
grid::grid.newpage()
grid::pushViewport(grid::viewport(layout = grid::grid.layout(nrow = 1, ncol = 2)))
print(h1_scree, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 1))
print(h1_loadings, vp = grid::viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# H2: Paired estimation plot for clone vs human on Valence and Dominance
h2_pair_data <- voice_identity %>%
  group_by(Participant.Public.ID, Type) %>%
  summarise(
    Valence = mean(Valence, na.rm = TRUE),
    Dominance = mean(Dominance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(Type = factor(Type, levels = c("human", "clone"))) %>%
  pivot_longer(
    cols = c(Valence, Dominance),
    names_to = "Dimension",
    values_to = "Score"
  )

star_from_p <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

h2_sig <- h2_pair_data %>%
  pivot_wider(names_from = Type, values_from = Score) %>%
  group_by(Dimension) %>%
  summarise(
    p = t.test(clone, human, paired = TRUE)$p.value,
    y = max(c(clone, human), na.rm = TRUE) + 0.20,
    stars = star_from_p(p),
    .groups = "drop"
  ) %>%
  filter(stars != "")

p_h2 <- ggplot(
  h2_pair_data,
  aes(x = Type, y = Score)
) +
  geom_violin(aes(fill = Type), alpha = 0.55, width = 0.85, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.35, fill = "white", alpha = 0.8) +
  geom_jitter(aes(color = Type), width = 0.12, height = 0, alpha = 0.35, size = 1.1) +
  geom_text(
    data = h2_sig,
    aes(x = 1.5, y = y, label = stars),
    inherit.aes = FALSE,
    size = 5
  ) +
  facet_wrap(~ Dimension, ncol = 2, scales = "fixed") +
  scale_fill_manual(values = c("human" = "#8A8A8A", "clone" = "#4D6FA8")) +
  scale_color_manual(values = c("human" = "#5F5F5F", "clone" = "#1E4A8A")) +
  labs(
    title = "H2: Clone vs human distributions",
    x = "",
    y = "PCA axis score"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(file.path(fig_dir, "figure_h2_average_type_valence_dominance.png"), p_h2, width = 8, height = 4, dpi = 300)
# END OF AI SECTION

# Some boxplots per voice of both PCs

voices_to_plot <- c(1,3,4,6,7,8,9,10,11,12)

ggplot(
  voice_identity |> dplyr::filter(Voice %in% voices_to_plot),
  aes(x = Type, y = Valence, fill = Type)
) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(~ Voice, ncol = 5) +
  labs(
    title = "Valence by Voice Pair (Human vs Clone)",
    x = "",
    y = "Valence (PC1)"
  ) +
  theme_minimal()

# Some basic stats
# Very basic t-tests to see whether the PC scores for both types differ and how

data_part <- voice_identity %>%
  group_by(Participant.Public.ID, Type) %>%
  summarise(
    Valence   = mean(Valence, na.rm = TRUE),
    Dominance = mean(Dominance, na.rm = TRUE),
    .groups   = "drop"
  ) %>%
  pivot_wider(
    names_from  = Type,
    values_from = c(Valence, Dominance)
  )

# Valence t-test
t.test(
  data_part$Valence_clone,
  data_part$Valence_human,
  paired = TRUE
)

# Dominance t-test
t.test(
  data_part$Dominance_clone,
  data_part$Dominance_human,
  paired = TRUE
)

# LMEs, the same as above but more sophisticated

# Ensure Type factor (human baseline)
voice_identity <- voice_identity %>%
  mutate(Type = factor(Type, levels = c("human", "clone")))

# Valence model
m_val <- lmer(
  Valence ~ Type + (1 | Participant.Public.ID) + (1 | Voice),
  data = voice_identity
)
summary(m_val)

# Dominance model
m_dom <- lmer(
  Dominance ~ Type + (1 | Participant.Public.ID) + (1 | Voice),
  data = voice_identity
)
summary(m_dom)

# Now we start to investigate which voices *drive* these differences
# First run some per‑voice paired t‑tests for Valence and Dominance
# This code collapses to participant × voice × type means for Valence and Dominance and 
# reshapes so each row has Human and Clone columns per voice (e.g., Valence_human, Valence_clone).
# Next, for each voice, it runs a paired t‑test across participants: Clone vs Human.

by_voice <- voice_identity %>%
  group_by(Participant.Public.ID, Voice, Type) %>%
  summarise(Valence = mean(Valence), Dominance = mean(Dominance), .groups = "drop") %>%
  pivot_wider(names_from = Type, values_from = c(Valence, Dominance))

valence_voice_tests <- by_voice %>%
  group_by(Voice) %>%
  summarise(
    diff = mean(Valence_clone - Valence_human, na.rm = TRUE),
    p    = t.test(Valence_clone, Valence_human, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_fdr = p.adjust(p, method = "BH"))

dominance_voice_tests <- by_voice %>%
  group_by(Voice) %>%
  summarise(
    diff = mean(Dominance_clone - Dominance_human, na.rm = TRUE),
    p    = t.test(Dominance_clone, Dominance_human, paired = TRUE)$p.value,
    .groups = "drop"
  ) %>%
  mutate(p_fdr = p.adjust(p, method = "BH"))

# Output per voice: mean difference (Clone − Human) and p‑value.
# Read it as: big |difference| + small p → that voice shows a strong Clone–Human effect.
valence_voice_tests
dominance_voice_tests

# Significance-starred per-voice plots (no star for non-significant voices)
sig_stars <- function(p) {
  ifelse(p < 0.001, "***",
         ifelse(p < 0.01, "**",
                ifelse(p < 0.05, "*", "")))
}

voice_order <- c(1, 3, 4, 6, 7, 8, 9, 10, 11, 12)
voice_map <- tibble(
  Voice = voice_order,
  VoiceLabel = seq_along(voice_order)
)

valence_sig_voices <- valence_voice_tests %>%
  filter(p_fdr < 0.05) %>%
  pull(Voice)

dominance_sig_voices <- dominance_voice_tests %>%
  filter(p_fdr < 0.05) %>%
  pull(Voice)

valence_plot_data <- voice_identity %>%
  filter(Voice %in% valence_sig_voices) %>%
  left_join(voice_map, by = "Voice") %>%
  mutate(VoiceLabel = factor(VoiceLabel, levels = seq_along(voice_order)))

dominance_plot_data <- voice_identity %>%
  filter(Voice %in% dominance_sig_voices) %>%
  left_join(voice_map, by = "Voice") %>%
  mutate(VoiceLabel = factor(VoiceLabel, levels = seq_along(voice_order)))

valence_ann <- valence_plot_data %>%
  group_by(Voice, VoiceLabel) %>%
  summarise(y = max(Valence, na.rm = TRUE) + 0.30, .groups = "drop") %>%
  left_join(valence_voice_tests, by = "Voice") %>%
  mutate(stars = sig_stars(p_fdr)) %>%
  filter(stars != "")

dominance_ann <- dominance_plot_data %>%
  group_by(Voice, VoiceLabel) %>%
  summarise(y = max(Dominance, na.rm = TRUE) + 0.30, .groups = "drop") %>%
  left_join(dominance_voice_tests, by = "Voice") %>%
  mutate(stars = sig_stars(p_fdr)) %>%
  filter(stars != "")

p_valence_stars <- ggplot(
  valence_plot_data,
  aes(x = Type, y = Valence, fill = Type)
) +
  geom_violin(alpha = 0.55, width = 0.85, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.35, fill = "white", alpha = 0.8) +
  geom_jitter(aes(color = Type), width = 0.12, height = 0, alpha = 0.35, size = 1.1) +
  geom_text(
    data = valence_ann,
    aes(x = 1.5, y = y, label = stars),
    inherit.aes = FALSE,
    size = 5
  ) +
  facet_wrap(~ VoiceLabel, ncol = 5) +
  scale_fill_manual(values = c("human" = "#8A8A8A", "clone" = "#4D6FA8")) +
  scale_color_manual(values = c("human" = "#5F5F5F", "clone" = "#1E4A8A")) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "",
    y = "Valence (PC1)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 15, margin = margin(t = 4, unit = "pt"))
  )

ggsave(file.path(fig_dir, "figure1_valence_all_voices.png"), p_valence_stars, width = 10.5, height = 5.5, dpi = 300)

p_dominance_stars <- ggplot(
  dominance_plot_data,
  aes(x = Type, y = Dominance, fill = Type)
) +
  geom_violin(alpha = 0.55, width = 0.85, trim = FALSE, colour = NA) +
  geom_boxplot(width = 0.15, outlier.shape = NA, linewidth = 0.35, fill = "white", alpha = 0.8) +
  geom_jitter(aes(color = Type), width = 0.12, height = 0, alpha = 0.35, size = 1.1) +
  geom_text(
    data = dominance_ann,
    aes(x = 1.5, y = y, label = stars),
    inherit.aes = FALSE,
    size = 5
  ) +
  facet_wrap(~ VoiceLabel, ncol = 5) +
  scale_fill_manual(values = c("human" = "#8A8A8A", "clone" = "#4D6FA8")) +
  scale_color_manual(values = c("human" = "#5F5F5F", "clone" = "#1E4A8A")) +
  labs(
    title = NULL,
    subtitle = NULL,
    x = "",
    y = "Dominance (PC2)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 15, margin = margin(t = 4, unit = "pt"))
  )

ggsave(file.path(fig_dir, "figure3_dominance_all_voices.png"), p_dominance_stars, width = 10.5, height = 5.5, dpi = 300)

# Leave One Out code LMEs (remove one vocie to see if significance of model changes)
# Fits  a baseline mixed model on all voices: Outcome ~ Type + (1|Participant) + (1|Voice).
# Records the Typeclone estimate (Clone − Human).
# Then loops: remove one voice, refit, record the new Typeclone estimate.
# Compares each new estimate to the baseline; the change (delta) shows that voice’s influence.
# Large absolute delta = that voice is driving the global Type effect.
# Sign of delta tells whether removing that voice makes the Clone–Human gap bigger or smaller.

voices <- unique(voice_identity$Voice)

m_val_all <- lmer(Valence ~ Type + (1 | Participant.Public.ID) + (1 | Voice), data = voice_identity)
m_dom_all <- lmer(Dominance ~ Type + (1 | Participant.Public.ID) + (1 | Voice), data = voice_identity)

base_val <- fixef(m_val_all)["Typeclone"]
base_dom <- fixef(m_dom_all)["Typeclone"]

loo_val <- sapply(voices, function(v){
  fit <- lmer(Valence ~ Type + (1 | Participant.Public.ID) + (1 | Voice),
              data = subset(voice_identity, Voice != v))
  fixef(fit)["Typeclone"]
})
loo_dom <- sapply(voices, function(v){
  fit <- lmer(Dominance ~ Type + (1 | Participant.Public.ID) + (1 | Voice),
              data = subset(voice_identity, Voice != v))
  fixef(fit)["Typeclone"]
})

loo_val <- data.frame(Voice = voices, estimate = as.numeric(loo_val),
                      delta = as.numeric(loo_val) - base_val)
loo_dom <- data.frame(Voice = voices, estimate = as.numeric(loo_dom),
                      delta = as.numeric(loo_dom) - base_dom)


voices <- unique(voice_identity$Voice)

m_val_all <- lmer(Valence ~ Type + (1 | Participant.Public.ID) + (1 | Voice), data = voice_identity)
m_dom_all <- lmer(Dominance ~ Type + (1 | Participant.Public.ID) + (1 | Voice), data = voice_identity)

base_val <- fixef(m_val_all)["Typeclone"]
base_dom <- fixef(m_dom_all)["Typeclone"]

loo_val <- sapply(voices, function(v){
  fit <- lmer(Valence ~ Type + (1 | Participant.Public.ID) + (1 | Voice),
              data = subset(voice_identity, Voice != v))
  fixef(fit)["Typeclone"]
})
loo_dom <- sapply(voices, function(v){
  fit <- lmer(Dominance ~ Type + (1 | Participant.Public.ID) + (1 | Voice),
              data = subset(voice_identity, Voice != v))
  fixef(fit)["Typeclone"]
})

loo_val <- data.frame(Voice = voices, estimate = as.numeric(loo_val),
                      delta = as.numeric(loo_val) - base_val)
loo_dom <- data.frame(Voice = voices, estimate = as.numeric(loo_dom),
                      delta = as.numeric(loo_dom) - base_dom)

# Large absolute delta = that voice is driving the global Type effect.
# Sign of delta tells whether removing that voice makes the Clone–Human gap bigger or smaller.

loo_val
loo_dom

# If you ask GenAI to interpret the result, it appears that it is voices 8 and 7 
# and to some degree 6 and 9 that drive the effects
# Here is a plot that shows the effects on the PCs again for four voices

voices_to_plot <- c(6,7,8,9)

ggplot(
  voice_identity |> dplyr::filter(Voice %in% voices_to_plot),
  aes(x = Type, y = Valence, fill = Type)
) +
  geom_violin(trim = FALSE, alpha = 0.5) +
  geom_boxplot(width = 0.15, outlier.shape = NA) +
  facet_wrap(~ Voice, ncol = 5) +
  labs(
    title = "Valence by Voice Pair (Human vs Clone)",
    x = "",
    y = "Valence (PC1)"
  ) +
  theme_minimal()

# Now we move onto the acoustic analysis, I chose five factors that 
# might give you the best results, pitch, pitch variation, intensity variation
# shimmer (voice roughness), and duration (speaking rate).

acoustic_data <- read.csv(file.path(project_dir, "Data", "Acoustics_short.csv"))  # read in Acoustics_short.csv

acoustic_data$shimmerLocal <- acoustic_data$shimmerLocal. # it had a stupid name

# reformat the data so you can measure differences per voice for 80 sentences
acoustic_diff <- acoustic_data %>%
  # set factors
  mutate(
    type  = factor(type, levels = c("human","clone")),
    voice = as.factor(voice)
  ) %>%
  # per-voice × type means over all sentences
  group_by(voice, type) %>%
  summarise(
    pitch         = mean(pitch, na.rm = TRUE),
    pitchdev      = mean(pitchdev, na.rm = TRUE),
    durationMsec  = mean(durationMsec, na.rm = TRUE),
    shimmerLocal  = mean(shimmerLocal, na.rm = TRUE),
    intensitydev  = mean(intensitydev, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  # wide by type
  pivot_wider(
    names_from  = type,
    values_from = c(pitch, pitchdev, durationMsec, shimmerLocal, intensitydev)
  ) %>%
  # clone − human differences
  mutate(
    d_pitch         = pitch_clone        - pitch_human,
    d_pitchdev      = pitchdev_clone     - pitchdev_human,
    d_durationMsec  = durationMsec_clone - durationMsec_human,
    d_shimmerLocal  = shimmerLocal_clone - shimmerLocal_human,
    d_intensitydev  = intensitydev_clone - intensitydev_human
  ) %>%
  select(voice, starts_with("d_")) %>%
  arrange(voice)

acoustic_diff

# Overall results indicate that:
# d_pitch: + = clone has higher pitch
# d_pitchdev: + = clone has more pitch variability (rougher / more dominant)
# d_durationMsec: + = clone speaks slower / longer
# d_shimmerLocal: + = clone has more shimmer (breathier / rougher)
# d_intensitydev: + = clone has more intensity fluctuation (more dynamic)

# Here are some voice-specific results that I got after dumping the table in copilot:

# VALENCE
# Voice 8
# 
# pitchdev: –25.4 (clone is much more unstable)
# shimmer: –0.0184 (clone more rough)
# duration: +88.1 (clone slower)
# → Strong driver of lower Valence
# 
# Voice 7
# 
# pitchdev: –16.8 (clone more unstable)
# shimmer: –0.00745 (clone rougher)
# duration: +72 (clone slower)
# → Also contributes to lower Valence
# 
# Voice 3
# 
# pitchdev: –15.1
# duration: –248 (clone faster)
# → Mixed, but pitchdev points to lower Valence
# 
# Voice 4
# 
# pitchdev −16, shimmer negative
# → Lower Valence trend

# DOMINANCE
# Voices contributing to higher Dominance
# (Dominance = clone sounds more forceful/assertive)
# Typically driven by:
# • higher pitch
# • more pitchdev
# • more intensitydev

# Voice 8
# 
# pitch: +4.74
# pitchdev: –25.4 (more variable)
# intensitydev: +0.267 (more dynamic)
# → Strong driver of higher Dominance
# 
# Voice 6
# 
# pitch: –9.45 (clone lower pitch)
# pitchdev: –4.18 (more unstable)
# duration: +177
# → Less clear, but duration + instability can push dominance perception
# 
# Voice 9
# 
# pitch: +3.98
# intensitydev: –0.0825 (slightly more stable)
# → Mildly dominant based on pitch

# ------------------------------
# 5) Plot A: raw-score lollipop by voice (shared x-axis across traits)
# ------------------------------
acoustic_long_all <- acoustic_diff %>%
  mutate(Voice = as.integer(as.character(voice))) %>%
  pivot_longer(
    cols = starts_with("d_"),
    names_to = "feature",
    values_to = "diff"
  ) %>%
  mutate(
    voice_f = factor(Voice, levels = sort(unique(Voice))),
    feature = factor(
      feature,
      levels = c("d_pitch", "d_pitchdev", "d_durationMsec", "d_shimmerLocal", "d_intensitydev"),
      labels = c("Delta pitch", "Delta pitch variability", "Delta duration (ms)", "Delta shimmer", "Delta intensity variability")
    )
  )

p_acoustic_lollipop_raw <- ggplot(acoustic_long_all, aes(y = feature, x = diff, color = feature)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_segment(aes(x = 0, xend = diff, y = feature, yend = feature), linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ voice_f, ncol = 5, scales = "fixed") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Raw acoustic clone-human differences by voice",
    subtitle = "Shared x-axis across traits; each point is a trait difference within a voice",
    x = "Difference (clone minus human, raw units)",
    y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

print(p_acoustic_lollipop_raw)
ggsave(file.path(fig_dir, "figure5_acoustic_raw_lollipop_by_voice.png"), p_acoustic_lollipop_raw, width = 13, height = 7, dpi = 300)

# ------------------------------
# 6) Table B: acoustic clone-human differences by voice
# ------------------------------
acoustic_table <- acoustic_diff %>%
  mutate(Voice = as.integer(as.character(voice))) %>%
  transmute(
    Voice,
    d_pitch,
    d_pitchdev,
    d_durationMsec,
    d_shimmerLocal,
    d_intensitydev
  ) %>%
  arrange(Voice)

print(acoustic_table)
write.csv(acoustic_table, file.path(out_dir, "table1_acoustic_differences_by_voice.csv"), row.names = FALSE)
