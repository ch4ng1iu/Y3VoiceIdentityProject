# Publication plots for the Voice Identity project.
# Run from project root: Rscript Code/01_generate_publication_plots.R

suppressPackageStartupMessages({
  library(dplyr)
  library(tidyr)
  library(ggplot2)
  library(grid)
})

# Resolve paths relative to this script for reproducibility.
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
out_dir <- file.path(project_dir, "Code", "Figures")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

ratings_path <- file.path(project_dir, "Data", "VoiceIdentity_forAnalysis (1).csv")
acoustics_path <- file.path(project_dir, "Data", "Acoustics_short.csv")

if (!file.exists(ratings_path)) stop("Missing data file: VoiceIdentity_forAnalysis (1).csv")
if (!file.exists(acoustics_path)) stop("Missing data file: Acoustics_short.csv")

Dat <- read.csv(ratings_path)
Dat <- Dat %>%
  rename(Stimulus = Spreadsheet..Audio) %>%
  mutate(Stimulus = gsub("\\.mp3$", "", Stimulus))

voice_identity <- Dat %>%
  pivot_wider(
    id_cols = c(Participant.Public.ID, Stimulus, Voice, Type),
    names_from = Trait,
    values_from = Response
  )

trait_cols <- c(
  "Trustworthiness", "Dominance", "Attractiveness", "Competence",
  "Warmth", "Aggression", "Confidence", "Friendliness"
)

if (!all(trait_cols %in% names(voice_identity))) {
  missing_traits <- paste(setdiff(trait_cols, names(voice_identity)), collapse = ", ")
  stop(paste("Missing trait columns:", missing_traits))
}

pca <- prcomp(voice_identity[trait_cols], scale. = TRUE)
voice_identity$Valence <- as.data.frame(pca$x)$PC1
voice_identity$Dominance <- as.data.frame(pca$x)$PC2
voice_identity$Type <- factor(voice_identity$Type, levels = c("human", "clone"))

star_from_p <- function(p) {
  ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
}

paired_p <- function(x, y) {
  out <- tryCatch(t.test(x, y, paired = TRUE)$p.value, error = function(e) NA_real_)
  out
}

# Figure 1: PCA scree + loadings.
scree_data <- data.frame(
  PC = seq_along(pca$sdev),
  Variance = (pca$sdev^2) / sum(pca$sdev^2)
) %>%
  mutate(
    Cumulative = cumsum(Variance),
    Highlight = ifelse(PC %in% c(1, 2), "PC1_PC2", "Other"),
    VarianceLabel = ifelse(PC %in% c(1, 2), paste0("PC", PC, ": ", sprintf("%.1f%%", 100 * Variance)), NA_character_)
  )

p1_scree <- ggplot(scree_data, aes(x = PC, y = Variance, fill = Highlight)) +
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
    title = "Panel A: PCA scree plot",
    x = "Principal component",
    y = "Proportion of variance"
  ) +
  theme_minimal(base_size = 11)

load_data <- as.data.frame(pca$rotation[, c("PC1", "PC2")])
load_data$Trait <- rownames(load_data)

p1_loadings <- ggplot(load_data, aes(x = PC1, y = PC2, label = Trait)) +
  geom_hline(yintercept = 0, colour = "grey70", linetype = 2) +
  geom_vline(xintercept = 0, colour = "grey70", linetype = 2) +
  geom_point(colour = "#1F5AA6", size = 2.2) +
  geom_text(nudge_y = 0.03, size = 3.1) +
  coord_cartesian(
    xlim = range(load_data$PC1) + c(-0.08, 0.08),
    ylim = range(load_data$PC2) + c(-0.08, 0.08),
    clip = "off"
  ) +
  labs(
    title = "Panel B: Trait loadings",
    x = "PC1 loading",
    y = "PC2 loading"
  ) +
  theme_minimal(base_size = 11) +
  theme(plot.margin = margin(8, 20, 8, 20))

png(file.path(out_dir, "Figure_01_H1_PCA_Scree_and_Loadings.png"), width = 3000, height = 1400, res = 300)
grid.newpage()
pushViewport(viewport(layout = grid.layout(nrow = 1, ncol = 2)))
print(p1_scree, vp = viewport(layout.pos.row = 1, layout.pos.col = 1))
print(p1_loadings, vp = viewport(layout.pos.row = 1, layout.pos.col = 2))
dev.off()

# Figure 2: clone vs human across PCA axes.
h2_data <- voice_identity %>%
  group_by(Participant.Public.ID, Type) %>%
  summarise(
    Valence = mean(Valence, na.rm = TRUE),
    Dominance = mean(Dominance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_longer(
    cols = c(Valence, Dominance),
    names_to = "Dimension",
    values_to = "Score"
  )

h2_sig <- h2_data %>%
  pivot_wider(names_from = Type, values_from = Score) %>%
  group_by(Dimension) %>%
  summarise(
    p = paired_p(clone, human),
    y = max(c(clone, human), na.rm = TRUE) + 0.20,
    stars = star_from_p(p),
    .groups = "drop"
  ) %>%
  filter(stars != "")

p2 <- ggplot(h2_data, aes(x = Type, y = Score)) +
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
    title = "Clone vs human distributions",
    x = "",
    y = "PCA axis score"
  ) +
  theme_minimal() +
  theme(legend.position = "none")

ggsave(
  filename = file.path(out_dir, "Figure_02_H2_Clone_vs_Human_PCA_Axes.png"),
  plot = p2,
  width = 8,
  height = 4,
  dpi = 300
)

# Figure 3 and 4: significant voices only.
by_voice <- voice_identity %>%
  group_by(Participant.Public.ID, Voice, Type) %>%
  summarise(
    Valence = mean(Valence, na.rm = TRUE),
    Dominance = mean(Dominance, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = Type, values_from = c(Valence, Dominance))

valence_voice_tests <- by_voice %>%
  group_by(Voice) %>%
  summarise(
    diff = mean(Valence_clone - Valence_human, na.rm = TRUE),
    p = paired_p(Valence_clone, Valence_human),
    .groups = "drop"
  ) %>%
  mutate(p_fdr = p.adjust(p, method = "BH"))

dominance_voice_tests <- by_voice %>%
  group_by(Voice) %>%
  summarise(
    diff = mean(Dominance_clone - Dominance_human, na.rm = TRUE),
    p = paired_p(Dominance_clone, Dominance_human),
    .groups = "drop"
  ) %>%
  mutate(p_fdr = p.adjust(p, method = "BH"))

voice_order <- sort(unique(voice_identity$Voice))
voice_map <- tibble(Voice = voice_order, VoiceLabel = seq_along(voice_order))

valence_sig_voices <- valence_voice_tests %>% filter(p_fdr < 0.05) %>% pull(Voice)
dominance_sig_voices <- dominance_voice_tests %>% filter(p_fdr < 0.05) %>% pull(Voice)

# Fallback: use all voices if none are significant after FDR.
if (length(valence_sig_voices) == 0) valence_sig_voices <- voice_order
if (length(dominance_sig_voices) == 0) dominance_sig_voices <- voice_order

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
  mutate(stars = star_from_p(p_fdr)) %>%
  filter(stars != "")

dominance_ann <- dominance_plot_data %>%
  group_by(Voice, VoiceLabel) %>%
  summarise(y = max(Dominance, na.rm = TRUE) + 0.30, .groups = "drop") %>%
  left_join(dominance_voice_tests, by = "Voice") %>%
  mutate(stars = star_from_p(p_fdr)) %>%
  filter(stars != "")

p3 <- ggplot(valence_plot_data, aes(x = Type, y = Valence, fill = Type)) +
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
    title = "Valence by voice",
    x = "",
    y = "Valence (PC1)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 13, margin = margin(t = 3, unit = "pt"))
  )

ggsave(
  filename = file.path(out_dir, "Figure_03_H3_Valence_By_Voice.png"),
  plot = p3,
  width = 10.5,
  height = 5.5,
  dpi = 300
)

p4 <- ggplot(dominance_plot_data, aes(x = Type, y = Dominance, fill = Type)) +
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
    title = "Dominance by voice",
    x = "",
    y = "Dominance (PC2)"
  ) +
  theme_minimal() +
  theme(
    legend.position = "none",
    strip.text = element_text(size = 13, margin = margin(t = 3, unit = "pt"))
  )

ggsave(
  filename = file.path(out_dir, "Figure_04_H3_Dominance_By_Voice.png"),
  plot = p4,
  width = 10.5,
  height = 5.5,
  dpi = 300
)

# Figure 5: acoustic clone-human differences by voice.
acoustic_data <- read.csv(acoustics_path)
if ("shimmerLocal." %in% names(acoustic_data)) {
  acoustic_data$shimmerLocal <- acoustic_data$shimmerLocal.
}

acoustic_diff <- acoustic_data %>%
  mutate(
    type = factor(type, levels = c("human", "clone")),
    voice = as.factor(voice)
  ) %>%
  group_by(voice, type) %>%
  summarise(
    pitch = mean(pitch, na.rm = TRUE),
    pitchdev = mean(pitchdev, na.rm = TRUE),
    durationMsec = mean(durationMsec, na.rm = TRUE),
    shimmerLocal = mean(shimmerLocal, na.rm = TRUE),
    intensitydev = mean(intensitydev, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  pivot_wider(names_from = type, values_from = c(pitch, pitchdev, durationMsec, shimmerLocal, intensitydev)) %>%
  mutate(
    Voice = as.integer(as.character(voice)),
    d_pitch = pitch_clone - pitch_human,
    d_pitchdev = pitchdev_clone - pitchdev_human,
    d_durationMsec = durationMsec_clone - durationMsec_human,
    d_shimmerLocal = shimmerLocal_clone - shimmerLocal_human,
    d_intensitydev = intensitydev_clone - intensitydev_human
  ) %>%
  select(Voice, d_pitch, d_pitchdev, d_durationMsec, d_shimmerLocal, d_intensitydev)

acoustic_long <- acoustic_diff %>%
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

p5 <- ggplot(acoustic_long, aes(y = feature, x = diff, color = feature)) +
  geom_vline(xintercept = 0, colour = "grey55", linetype = 2) +
  geom_segment(aes(x = 0, xend = diff, y = feature, yend = feature), linewidth = 0.8) +
  geom_point(size = 2.2) +
  facet_wrap(~ voice_f, ncol = 5, scales = "fixed") +
  scale_color_brewer(palette = "Set2") +
  labs(
    title = "Acoustic clone-human differences by voice",
    subtitle = "Shared x-axis across acoustic features",
    x = "Difference (clone minus human; raw units)",
    y = ""
  ) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom")

ggsave(
  filename = file.path(out_dir, "Figure_05_H3_Acoustic_Differences_By_Voice.png"),
  plot = p5,
  width = 13,
  height = 7,
  dpi = 300
)

message("Saved publication figures to: ", out_dir)
