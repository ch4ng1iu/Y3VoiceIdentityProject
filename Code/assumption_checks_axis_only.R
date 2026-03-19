library(dplyr)
library(tidyr)
library(lme4)
library(lmerTest)
library(psych)

set.seed(20260225)

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

Dat <- read.csv(file.path(project_dir, "Data", "VoiceIdentity_forAnalysis (1).csv"))
Dat <- Dat %>%
  rename(Stimulus = Spreadsheet..Audio) %>%
  mutate(Stimulus = gsub('\\.mp3$', '', Stimulus))

voice_identity <- Dat %>%
  pivot_wider(
    id_cols = c(Participant.Public.ID, Stimulus, Voice, Type),
    names_from = Trait,
    values_from = Response
  )

trait_cols <- c('Trustworthiness', 'Dominance', 'Attractiveness', 'Competence',
                'Warmth', 'Aggression', 'Confidence', 'Friendliness')

# -----------------------------
# H1 assumption checks (PCA)
# -----------------------------
R <- cor(voice_identity[trait_cols], use = 'pairwise.complete.obs')
kmo <- KMO(R)
bart <- cortest.bartlett(R, n = nrow(voice_identity))
fa_par <- fa.parallel(voice_identity[trait_cols], fa = 'pc', n.iter = 100, plot = FALSE)

# PCA scores for H2/H3
pca <- prcomp(voice_identity[trait_cols], scale = TRUE)
voice_identity$Valence <- as.data.frame(pca$x)$PC1
voice_identity$DominancePC <- as.data.frame(pca$x)$PC2
voice_identity$Type <- factor(voice_identity$Type, levels = c('human', 'clone'))

# -----------------------------
# H2 assumption checks
# -----------------------------
data_part <- voice_identity %>%
  group_by(Participant.Public.ID, Type) %>%
  summarise(
    Valence = mean(Valence, na.rm = TRUE),
    Dominance = mean(DominancePC, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_wider(names_from = Type, values_from = c(Valence, Dominance))

val_diff <- data_part$Valence_clone - data_part$Valence_human
dom_diff <- data_part$Dominance_clone - data_part$Dominance_human

sw_val <- shapiro.test(val_diff)
sw_dom <- shapiro.test(dom_diff)

outlier_n <- function(x) {
  q <- quantile(x, probs = c(0.25, 0.75), na.rm = TRUE)
  iqr <- q[2] - q[1]
  lo <- q[1] - 1.5 * iqr
  hi <- q[2] + 1.5 * iqr
  sum(x < lo | x > hi, na.rm = TRUE)
}

m_val <- lmer(Valence ~ Type + (1 | Participant.Public.ID) + (1 | Voice), data = voice_identity)
m_dom <- lmer(DominancePC ~ Type + (1 | Participant.Public.ID) + (1 | Voice), data = voice_identity)

# -----------------------------
# H3 assumption checks (axis-only)
# -----------------------------
by_voice <- voice_identity %>%
  group_by(Participant.Public.ID, Voice, Type) %>%
  summarise(
    Valence = mean(Valence, na.rm = TRUE),
    Dominance = mean(DominancePC, na.rm = TRUE),
    .groups = 'drop'
  ) %>%
  pivot_wider(names_from = Type, values_from = c(Valence, Dominance))

valence_voice_tests <- by_voice %>%
  group_by(Voice) %>%
  summarise(
    diff = mean(Valence_clone - Valence_human, na.rm = TRUE),
    p = t.test(Valence_clone, Valence_human, paired = TRUE)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(p_fdr = p.adjust(p, method = 'BH'))

dominance_voice_tests <- by_voice %>%
  group_by(Voice) %>%
  summarise(
    diff = mean(Dominance_clone - Dominance_human, na.rm = TRUE),
    p = t.test(Dominance_clone, Dominance_human, paired = TRUE)$p.value,
    .groups = 'drop'
  ) %>%
  mutate(p_fdr = p.adjust(p, method = 'BH'))

acoustic_data <- read.csv(file.path(project_dir, "Data", "Acoustics_short.csv"))
acoustic_data$shimmerLocal <- acoustic_data$shimmerLocal.
acoustic_diff <- acoustic_data %>%
  mutate(type = factor(type, levels = c('human', 'clone')), voice = as.factor(voice)) %>%
  group_by(voice, type) %>%
  summarise(
    pitch = mean(pitch, na.rm = TRUE),
    pitchdev = mean(pitchdev, na.rm = TRUE),
    durationMsec = mean(durationMsec, na.rm = TRUE),
    shimmerLocal = mean(shimmerLocal, na.rm = TRUE),
    intensitydev = mean(intensitydev, na.rm = TRUE),
    .groups = 'drop'
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

h3 <- left_join(
  valence_voice_tests %>% select(Voice, d_val = diff),
  dominance_voice_tests %>% select(Voice, d_dom = diff),
  by = 'Voice'
) %>%
  left_join(acoustic_diff, by = 'Voice')

features <- c('d_pitch', 'd_pitchdev', 'd_durationMsec', 'd_shimmerLocal', 'd_intensitydev')
outcomes <- c('d_val', 'd_dom')

res <- do.call(rbind, lapply(features, function(f) {
  do.call(rbind, lapply(outcomes, function(o) {
    ct <- suppressWarnings(cor.test(h3[[f]], h3[[o]], method = 'spearman', exact = FALSE))
    data.frame(
      feature = f,
      outcome = o,
      rho = unname(ct$estimate),
      p_asym = ct$p.value
    )
  }))
})) %>%
  mutate(
    p_asym_fdr = p.adjust(p_asym, method = 'BH')
  )

# -----------------------------
# Output summary
# -----------------------------
cat('H1_CHECKS\n')
cat(sprintf('KMO_overall=%.6f\n', kmo$MSA))
cat(sprintf('Bartlett_chisq=%.6f df=%d p=%.12f\n', unname(bart$chisq), unname(bart$df), unname(bart$p.value)))
cat(sprintf('Parallel_ncomp=%d\n', fa_par$ncomp))

cat('H2_CHECKS\n')
cat(sprintf('Shapiro_ValenceDiff_W=%.6f p=%.12f\n', unname(sw_val$statistic), sw_val$p.value))
cat(sprintf('Shapiro_DominanceDiff_W=%.6f p=%.12f\n', unname(sw_dom$statistic), sw_dom$p.value))
cat(sprintf('Outliers_ValenceDiff=%d\n', outlier_n(val_diff)))
cat(sprintf('Outliers_DominanceDiff=%d\n', outlier_n(dom_diff)))
cat(sprintf('LME_Valence_singular=%s converged=%s\n', as.character(isSingular(m_val)), as.character(is.null(m_val@optinfo$conv$lme4$messages))))
cat(sprintf('LME_Dominance_singular=%s converged=%s\n', as.character(isSingular(m_dom)), as.character(is.null(m_dom@optinfo$conv$lme4$messages))))

cat('H3_CHECKS\n')
cat(sprintf('Valence_voice_tests_sig_FDR=%d_of_%d\n', sum(valence_voice_tests$p_fdr < 0.05), nrow(valence_voice_tests)))
cat(sprintf('Dominance_voice_tests_sig_FDR=%d_of_%d\n', sum(dominance_voice_tests$p_fdr < 0.05), nrow(dominance_voice_tests)))
cat(sprintf('Spearman_total=%d asym_sig=%d asym_sig_FDR=%d\n',
            nrow(res),
            sum(res$p_asym < 0.05),
            sum(res$p_asym_fdr < 0.05)))

best <- res[order(res$p_asym), ][1, ]
cat(sprintf('Top_acoustic_axis_assoc=%s_with_%s rho=%.6f p_asym=%.6f\n',
            best$feature, best$outcome, best$rho, best$p_asym))
