# ==========================================================
# Pupillometry: 2×2 models, Harmony only model, and the plotting stuff
# Ali Çağan Kaya, 2025
# alicagankaya.com
# ==========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(broom)
  library(broom.mixed)
  library(lmtest)
  library(sandwich)
  library(emmeans)
  library(ggplot2)
  library(tidyr)
})

# ---------- PATHS ----------
IN_FILE  <- "/Users/path/file.csv"
OUT_DIR  <- "/Users/path/analysis_out"
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)

# ---------- UTILS ----------
save_txt <- function(lines, file) {
  con <- file(file, open = "wt", encoding = "UTF-8")
  on.exit(close(con), add = TRUE)
  writeLines(lines, con)
}

ggsave_safe <- function(plot_obj, filename, width = 8, height = 5, dpi = 300) {
  # only save if it's a ggplot object
  if (inherits(plot_obj, "ggplot")) {
    ggsave(filename = filename, plot = plot_obj, width = width, height = height, dpi = dpi, limitsize = FALSE)
  }
}

# ---------- READ & CLEAN ----------
cat("\n=== Loading data ===\n")
DT <- fread(IN_FILE)

need_cols <- c("Subject","Trial","PupilBc","TimeRel","HarmonyType","SuffixType")
stopifnot(all(need_cols %in% names(DT)))
DT <- DT[, ..need_cols]

DT[, HarmonyType := tolower(trimws(as.character(HarmonyType)))]
DT[, SuffixType   := tolower(trimws(as.character(SuffixType)))]

DT <- DT[ HarmonyType %in% c("harmonic","disharmonic") &
            SuffixType %in% c("inf","der") ]

DT <- DT[!is.na(PupilBc) & !is.na(TimeRel)]
DT[, Subject := as.factor(Subject)]
DT[, Trial   := as.factor(Trial)]

# ---------- ANALYSIS WINDOWS ----------
windows <- list(
  Overall_700_2000 = c(700, 2000),
  Early_700_1200   = c(700, 1200),
  Late_1200_2000   = c(1200, 2000)
)

# ---------- HELPERS ----------
aggregate_window <- function(DT, tmin, tmax, min_samples = 5) {
  sub <- DT[TimeRel >= tmin & TimeRel <= tmax]
  agg <- sub[, .(
    pupil_mean = mean(PupilBc),
    pupil_med  = median(PupilBc),
    n_samples  = .N
  ), by = .(Subject, Trial, HarmonyType, SuffixType)]
  agg <- agg[n_samples >= min_samples]
  agg[, HarmonyType := factor(HarmonyType, levels = c("harmonic","disharmonic"))]
  agg[, SuffixType   := factor(SuffixType,   levels = c("inf","der"))]
  agg[]
}

print_counts <- function(agg, title) {
  cat("\n---", title, ": trial counts by condition ---\n")
  cnt <- agg[, .N, by = .(HarmonyType, SuffixType)][order(HarmonyType, SuffixType)]
  print(cnt)
  invisible(cnt)
}

fit_lmm_2x2 <- function(agg, win_name) {
  cat("\n================ LMM (lmerTest) =================\n")
  cat("Window:", win_name, "\n")
  m <- lmer(pupil_mean ~ HarmonyType * SuffixType + (1 | Subject), data = agg, REML = TRUE)
  print(summary(m), corr = FALSE)
  cat("\nFixed effects (tidy):\n")
  print(broom.mixed::tidy(m, effects = "fixed"))
  invisible(m)
}

fit_ols_cluster_2x2 <- function(agg, win_name) {
  cat("\n================ OLS with cluster-robust SEs =================\n")
  cat("Window:", win_name, "\n")
  lm0 <- lm(pupil_mean ~ HarmonyType * SuffixType, data = agg)
  Vc  <- sandwich::vcovCL(lm0, cluster = agg$Subject)
  ct  <- lmtest::coeftest(lm0, vcov. = Vc)
  
  cat("\nOLS (non-robust) summary — coefficients only:\n")
  print(summary(lm0))
  cat("\nOLS with cluster-robust SE by Subject:\n")
  print(ct)
  cat("\nCluster-robust 95% CIs:\n")
  est <- coef(lm0); se <- sqrt(diag(Vc))
  ci  <- cbind(
    estimate = est,
    ci_low   = est + qnorm(0.025) * se,
    ci_high  = est + qnorm(0.975) * se
  )
  print(ci)
  invisible(list(lm = lm0, vc = Vc, ct = ct))
}

fit_harmony_only <- function(agg, win_name) {
  cat("\n================ Harmony-only models =================\n")
  cat("Window:", win_name, "\n")
  m_lmm <- lmer(pupil_mean ~ HarmonyType + (1 | Subject), data = agg, REML = TRUE)
  cat("\nLMM summary (Harmony only):\n")
  print(summary(m_lmm), corr = FALSE)
  
  em <- emmeans(m_lmm, ~ HarmonyType)
  cat("\nEstimated marginal means (LMM):\n")
  print(em)
  cat("\nContrast (disharmonic - harmonic) — LMM:\n")
  print(contrast(em, "revpairwise"))
  
  lm0 <- lm(pupil_mean ~ HarmonyType, data = agg)
  Vc  <- sandwich::vcovCL(lm0, cluster = agg$Subject)
  ct  <- lmtest::coeftest(lm0, vcov. = Vc)
  cat("\nOLS with cluster-robust SEs (Harmony only):\n")
  print(ct)
  
  est <- coef(lm0); se <- sqrt(diag(Vc))
  ci  <- cbind(estimate = est,
               ci_low   = est + qnorm(0.025) * se,
               ci_high  = est + qnorm(0.975) * se)
  cat("\nCluster-robust 95% CIs (Harmony only):\n")
  print(ci)
  
  invisible(list(lmm = m_lmm, lm = lm0, vc = Vc, ct = ct))
}

# ---------- PLOTTING HELPERS (now RETURN plots) ----------
plot_timecourse <- function(DT, tmin = 700, tmax = 2000, bin_ms = 20) {
  lab_suffix  <- c(inf = "Inflectional", der = "Derivational")
  lab_harmony <- c(harmonic = "Harmonic", disharmonic = "Disharmonic")
  
  D <- DT %>%
    filter(TimeRel >= tmin, TimeRel <= tmax) %>%
    transmute(
      Subject = Subject, Trial = Trial,
      Pupil = PupilBc,
      TimeRel = TimeRel,
      Suffix = SuffixType,
      Harmony = HarmonyType
    ) %>%
    mutate(
      Suffix = factor(Suffix, levels = c("inf","der")),
      Harmony = factor(Harmony, levels = c("harmonic","disharmonic")),
      TimeBin = floor(TimeRel / bin_ms) * bin_ms
    )
  
  subj_time <- D %>%
    group_by(Subject, Suffix, Harmony, TimeBin) %>%
    summarise(pupil = mean(Pupil), .groups = "drop")
  
  time_summ <- subj_time %>%
    group_by(Suffix, Harmony, TimeBin) %>%
    summarise(
      mean = mean(pupil),
      sd   = sd(pupil),
      n    = dplyr::n(),
      se   = sd / sqrt(n),
      ci   = qt(0.975, df = pmax(n - 1, 1)) * se,
      lo   = mean - ci,
      hi   = mean + ci,
      .groups = "drop"
    )
  
  p1 <- ggplot(time_summ, aes(TimeBin, mean, color = Harmony, fill = Harmony)) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, linewidth = 0) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ Suffix, labeller = as_labeller(lab_suffix), nrow = 1) +
    labs(
      title = "Pupil time-course by Harmony × Suffix",
      subtitle = sprintf("Mean ± 95%% CI across subjects, %d–%d ms (bin = %d ms)", tmin, tmax, bin_ms),
      x = "Time relative to event (ms)",
      y = "Baseline-corrected pupil (a.u.)",
      color = "Harmony", fill = "Harmony"
    ) +
    scale_color_discrete(labels = lab_harmony) +
    scale_fill_discrete(labels = lab_harmony) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
  
  early_end <- 1200
  p2 <- ggplot(time_summ, aes(TimeBin, mean, color = Harmony, fill = Harmony)) +
    annotate("rect", xmin = tmin, xmax = early_end, ymin = -Inf, ymax = Inf, alpha = 0.05) +
    annotate("rect", xmin = early_end, xmax = tmax, ymin = -Inf, ymax = Inf, alpha = 0.08) +
    geom_ribbon(aes(ymin = lo, ymax = hi), alpha = 0.15, linewidth = 0) +
    geom_line(linewidth = 0.9) +
    facet_wrap(~ Suffix, labeller = as_labeller(lab_suffix), nrow = 1) +
    labs(
      title = "Pupil time-course with Early (700–1200) and Late (1200–2000) windows",
      x = "Time (ms)", y = "Baseline-corrected pupil (a.u.)"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
  
  return(list(timecourse = p1, timecourse_shaded = p2))
}

plot_windowed_interaction <- function(DT) {
  lab_suffix  <- c(inf = "Inflectional", der = "Derivational")
  
  W_OVERALL <- DT %>% filter(TimeRel >= 700,  TimeRel <= 2000) %>% mutate(Window = "Overall")
  W_EARLY   <- DT %>% filter(TimeRel >= 700,  TimeRel <= 1200) %>% mutate(Window = "Early")
  W_LATE    <- DT %>% filter(TimeRel >= 1200, TimeRel <= 2000) %>% mutate(Window = "Late")
  DW <- bind_rows(W_OVERALL, W_EARLY, W_LATE)
  
  trial_means <- DW %>%
    group_by(Window, Subject, Trial, Suffix = SuffixType, Harmony = HarmonyType) %>%
    summarise(pupil = mean(PupilBc), .groups = "drop") %>%
    mutate(Suffix = factor(Suffix, levels = c("inf","der")),
           Harmony = factor(Harmony, levels = c("harmonic","disharmonic")))
  
  subj_means <- trial_means %>%
    group_by(Window, Subject, Suffix, Harmony) %>%
    summarise(pupil = mean(pupil), .groups = "drop")
  
  win_summ <- subj_means %>%
    group_by(Window, Suffix, Harmony) %>%
    summarise(
      mean = mean(pupil),
      sd   = sd(pupil),
      n    = dplyr::n(),
      se   = sd / sqrt(n),
      ci   = qt(0.975, df = pmax(n - 1, 1)) * se,
      lo   = mean - ci,
      hi   = mean + ci,
      .groups = "drop"
    )
  
  p_int <- ggplot(win_summ, aes(Harmony, mean, group = Harmony)) +
    geom_point(position = position_dodge(width = 0.3), size = 2.2) +
    geom_errorbar(aes(ymin = lo, ymax = hi), width = 0.12,
                  position = position_dodge(width = 0.3)) +
    geom_line(aes(x = as.numeric(Harmony)), alpha = 0.35) +
    facet_grid(Window ~ Suffix, labeller = labeller(Suffix = lab_suffix)) +
    labs(
      title = "Windowed means by Harmony × Suffix",
      y = "Mean pupil (a.u.)", x = "Harmony"
    ) +
    theme_minimal(base_size = 13)
  
  p_dist <- ggplot(subj_means, aes(Harmony, pupil, fill = Harmony)) +
    geom_violin(width = 0.9, alpha = 0.2, color = NA) +
    geom_boxplot(width = 0.25, outlier.shape = 16, outlier.size = 1.5) +
    facet_grid(Window ~ Suffix, labeller = labeller(Suffix = lab_suffix)) +
    labs(
      title = "Subject-level distribution by window",
      y = "Mean pupil per subject (a.u.)", x = "Harmony"
    ) +
    theme_minimal(base_size = 13) +
    theme(panel.grid.minor = element_blank())
  
  p_spag <- ggplot(subj_means, aes(x = Harmony, y = pupil, group = Subject)) +
    geom_line(alpha = 0.25) +
    geom_point(alpha = 0.6, position = position_jitter(width = 0.03, height = 0)) +
    facet_grid(Window ~ Suffix, labeller = labeller(Suffix = lab_suffix)) +
    labs(
      title = "Subject trajectories: harmonic → disharmonic (per window × suffix)",
      y = "Mean pupil per subject (a.u.)", x = "Harmony"
    ) +
    theme_minimal(base_size = 13)
  
  p_bar <- ggplot(win_summ, aes(x = Suffix, y = mean, fill = Harmony)) +
    geom_col(position = position_dodge(width = 0.8), width = 0.7) +
    geom_errorbar(aes(ymin = lo, ymax = hi),
                  position = position_dodge(width = 0.8), width = 0.2) +
    facet_wrap(~ Window, nrow = 1) +
    labs(
      title = "Grouped means (±95% CI) by Suffix and Harmony across windows",
      x = "Suffix", y = "Mean pupil (a.u.)", fill = "Harmony"
    ) +
    scale_x_discrete(labels = c(inf = "Inflectional", der = "Derivational")) +
    theme_minimal(base_size = 13)
  
  return(list(window_means = p_int, subject_dist = p_dist, spaghetti = p_spag, grouped_bar = p_bar))
}

plot_simple_effects <- function(DT, which_window = "Overall") {
  lab_suffix  <- c(inf = "Inflectional", der = "Derivational")
  
  if (which_window == "Overall") {
    D <- DT %>% filter(TimeRel >= 700, TimeRel <= 2000) %>% mutate(Window = "Overall")
  } else if (which_window == "Early") {
    D <- DT %>% filter(TimeRel >= 700, TimeRel <= 1200) %>% mutate(Window = "Early")
  } else if (which_window == "Late") {
    D <- DT %>% filter(TimeRel >= 1200, TimeRel <= 2000) %>% mutate(Window = "Late")
  } else stop("which_window must be one of: Overall, Early, Late")
  
  trial_means <- D %>%
    group_by(Subject, Trial, Suffix = SuffixType, Harmony = HarmonyType) %>%
    summarise(pupil = mean(PupilBc), .groups = "drop") %>%
    mutate(Suffix = factor(Suffix, levels = c("inf","der")),
           Harmony = factor(Harmony, levels = c("harmonic","disharmonic")))
  
  subj_means <- trial_means %>%
    group_by(Subject, Suffix, Harmony) %>%
    summarise(pupil = mean(pupil), .groups = "drop")
  
  diff_df <- subj_means %>%
    select(Subject, Suffix, Harmony, pupil) %>%
    tidyr::pivot_wider(names_from = Harmony, values_from = pupil) %>%
    mutate(diff = disharmonic - harmonic)
  
  diff_summ <- diff_df %>%
    group_by(Suffix) %>%
    summarise(
      mean = mean(diff, na.rm = TRUE),
      sd   = sd(diff, na.rm = TRUE),
      n    = dplyr::n(),
      se   = sd / sqrt(n),
      ci   = qt(0.975, df = pmax(n - 1, 1)) * se,
      lo   = mean - ci,
      hi   = mean + ci,
      .groups = "drop"
    )
  
  p <- ggplot(diff_df, aes(Suffix, diff)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_jitter(width = 0.08, height = 0, alpha = 0.5, size = 1.8) +
    geom_point(data = diff_summ,
               aes(x = Suffix, y = mean),
               size = 3, inherit.aes = FALSE) +
    geom_errorbar(data = diff_summ,
                  aes(x = Suffix, ymin = lo, ymax = hi),
                  width = 0.15, inherit.aes = FALSE) +
    scale_x_discrete(labels = lab_suffix) +
    labs(
      title = paste0("Simple effects: (disharmonic − harmonic) by Suffix (", which_window, " window)"),
      y = "Difference in pupil (a.u.)", x = NULL
    ) +
    theme_minimal(base_size = 13)
  return(p)
}

fit_suffix_within_disharmonic <- function(agg, win_name) {
  cat("\n================ Disharmonic only: Suffix comparison =================\n")
  cat("Window:", win_name, "\n")
  
  A <- agg %>% filter(HarmonyType == "disharmonic")
  if (nrow(A) == 0) { cat("No disharmonic rows.\n"); return(invisible(NULL)) }
  
  m_lmm <- lmer(pupil_mean ~ SuffixType + (1 | Subject), data = A, REML = TRUE)
  cat("\nLMM summary (Disharmonic; Suffix only):\n")
  print(summary(m_lmm), corr = FALSE)
  
  em <- emmeans(m_lmm, ~ SuffixType)
  cat("\nEstimated marginal means (Disharmonic; LMM):\n")
  print(em)
  cat("\nContrast (der - inf) — Disharmonic; LMM:\n")
  print(contrast(em, "revpairwise"))
  
  lm0 <- lm(pupil_mean ~ SuffixType, data = A)
  Vc  <- sandwich::vcovCL(lm0, cluster = A$Subject)
  ct  <- lmtest::coeftest(lm0, vcov. = Vc)
  cat("\nOLS with cluster-robust SEs (Disharmonic; Suffix only):\n")
  print(ct)
  
  est <- coef(lm0); se <- sqrt(diag(Vc))
  ci  <- cbind(estimate = est,
               ci_low   = est + qnorm(0.025) * se,
               ci_high  = est + qnorm(0.975) * se)
  cat("\nCluster-robust 95% CIs (Disharmonic; Suffix only):\n")
  print(ci)
  
  invisible(list(lmm = m_lmm, lm = lm0))
}

plot_disharmonic_suffix <- function(DT, which_window = "Overall") {
  if (which_window == "Overall") {
    D <- DT %>% filter(TimeRel >= 700, TimeRel <= 2000) %>% mutate(Window = "Overall")
  } else if (which_window == "Early") {
    D <- DT %>% filter(TimeRel >= 700, TimeRel <= 1200) %>% mutate(Window = "Early")
  } else if (which_window == "Late") {
    D <- DT %>% filter(TimeRel >= 1200, TimeRel <= 2000) %>% mutate(Window = "Late")
  } else stop("which_window must be one of: Overall, Early, Late")
  
  trial_means <- D %>%
    group_by(Subject, Trial, Suffix = SuffixType, Harmony = HarmonyType) %>%
    summarise(pupil = mean(PupilBc), .groups = "drop")
  
  dish <- trial_means %>%
    filter(Harmony == "disharmonic") %>%
    group_by(Subject, Suffix) %>%
    summarise(pupil = mean(pupil), .groups = "drop") %>%
    mutate(Suffix = factor(Suffix, levels = c("inf","der")))
  
  summ <- dish %>%
    group_by(Suffix) %>%
    summarise(
      mean = mean(pupil), sd = sd(pupil), n = dplyr::n(),
      se = sd/sqrt(n), ci = qt(0.975, df = pmax(n - 1, 1)) * se,
      lo = mean - ci, hi = mean + ci, .groups = "drop"
    )
  
  p_pairs <- ggplot(dish, aes(x = Suffix, y = pupil, group = Subject)) +
    geom_line(alpha = 0.25) +
    geom_point(alpha = 0.6, position = position_jitter(width = 0.03, height = 0)) +
    geom_point(data = summ, aes(x = Suffix, y = mean), size = 3, inherit.aes = FALSE) +
    geom_errorbar(data = summ, aes(x = Suffix, ymin = lo, ymax = hi),
                  width = 0.15, inherit.aes = FALSE) +
    scale_x_discrete(labels = c(inf = "Inflectional", der = "Derivational")) +
    labs(
      title = paste0("Disharmonic only: Suffix comparison (", which_window, " window)"),
      x = "Suffix", y = "Mean pupil (a.u.)"
    ) +
    theme_minimal(base_size = 13)
  
  wide <- dish %>%
    tidyr::pivot_wider(names_from = Suffix, values_from = pupil) %>%
    mutate(diff = der - inf)
  
  summ_diff <- wide %>%
    summarise(
      mean = mean(diff, na.rm = TRUE),
      sd   = sd(diff, na.rm = TRUE),
      n    = dplyr::n(),
      se   = sd / sqrt(n),
      ci   = qt(0.975, df = pmax(n - 1, 1)) * se,
      lo   = mean - ci,
      hi   = mean + ci
    )
  
  p_diff <- ggplot(wide, aes(x = 1, y = diff)) +
    geom_hline(yintercept = 0, linetype = 2) +
    geom_jitter(width = 0.08, height = 0, alpha = 0.5, size = 1.8) +
    geom_point(data = summ_diff, aes(x = 1, y = mean), size = 3, inherit.aes = FALSE) +
    geom_errorbar(data = summ_diff, aes(x = 1, ymin = lo, ymax = hi),
                  width = 0.15, inherit.aes = FALSE) +
    scale_x_continuous(breaks = 1, labels = "Disharmonic\n(der − inf)") +
    labs(
      title = paste0("Disharmonic only: difference (der − inf) per subject (", which_window, " window)"),
      x = NULL, y = "Difference in pupil (a.u.)"
    ) +
    theme_minimal(base_size = 13)
  
  return(list(dish_suffix_pairs = p_pairs, dish_suffix_diff = p_diff))
}

# ---------- RUN MODELS + SAVE OUTPUTS ----------
cat("\n=== Modeling start ===\n")
for (nm in names(windows)) {
  rng <- windows[[nm]]
  agg <- aggregate_window(DT, tmin = rng[1], tmax = rng[2], min_samples = 5)
  
  # counts (print + save CSV)
  cnt <- print_counts(agg, nm)
  fwrite(cnt, file.path(OUT_DIR, sprintf("counts_%s.csv", nm)))
  
  # Save model outputs to txt per window
  m2_lmm_log <- capture.output({ m2_lmm <- fit_lmm_2x2(agg, nm) })
  save_txt(m2_lmm_log, file.path(OUT_DIR, sprintf("models_LMM_2x2_%s.txt", nm)))
  
  m2_cr_log  <- capture.output({ m2_cr <- fit_ols_cluster_2x2(agg, nm) })
  save_txt(m2_cr_log,  file.path(OUT_DIR, sprintf("models_OLScluster_2x2_%s.txt", nm)))
  
  mh_log <- capture.output({
    cat("\n--- Harmony-only analysis:", nm, "---\n")
    mh <- fit_harmony_only(agg, nm)
  })
  save_txt(mh_log, file.path(OUT_DIR, sprintf("models_HarmonyOnly_%s.txt", nm)))
  
  dish_log <- capture.output({ fit_suffix_within_disharmonic(agg, nm) })
  save_txt(dish_log, file.path(OUT_DIR, sprintf("models_Disharmonic_SuffixOnly_%s.txt", nm)))
}

cat("\n=== Plots start & savving ===\n")

# Time-course (2 plots)
tc <- plot_timecourse(DT, tmin = 700, tmax = 2000, bin_ms = 20)
ggsave_safe(tc$timecourse,        file.path(OUT_DIR, "plot_timecourse_700_2000.png"), width = 9.5, height = 4.8)
ggsave_safe(tc$timecourse_shaded, file.path(OUT_DIR, "plot_timecourse_700_2000_shaded.png"), width = 9.5, height = 4.8)

# Windowed interaction (4 plots)
wi <- plot_windowed_interaction(DT)
ggsave_safe(wi$window_means, file.path(OUT_DIR, "plot_windowed_means.png"), width = 9.5, height = 7)
ggsave_safe(wi$subject_dist, file.path(OUT_DIR, "plot_subject_distribution.png"), width = 9.5, height = 7)
ggsave_safe(wi$spaghetti,    file.path(OUT_DIR, "plot_spaghetti_subjects.png"), width = 9.5, height = 7)
ggsave_safe(wi$grouped_bar,  file.path(OUT_DIR, "plot_grouped_means_bar.png"), width = 9.5, height = 4.8)

# Simple-effects per window (3 plots)
se_overall <- plot_simple_effects(DT, which_window = "Overall")
se_early   <- plot_simple_effects(DT, which_window = "Early")
se_late    <- plot_simple_effects(DT, which_window = "Late")
ggsave_safe(se_overall, file.path(OUT_DIR, "plot_simple_effects_Overall.png"))
ggsave_safe(se_early,   file.path(OUT_DIR, "plot_simple_effects_Early.png"))
ggsave_safe(se_late,    file.path(OUT_DIR, "plot_simple_effects_Late.png"))

# Disharmonic-only suffix comparisons per window (2 x 3 plots)
ds_overall <- plot_disharmonic_suffix(DT, which_window = "Overall")
ds_early   <- plot_disharmonic_suffix(DT, which_window = "Early")
ds_late    <- plot_disharmonic_suffix(DT, which_window = "Late")
ggsave_safe(ds_overall$dish_suffix_pairs, file.path(OUT_DIR, "plot_dish_suffix_pairs_Overall.png"))
ggsave_safe(ds_overall$dish_suffix_diff,  file.path(OUT_DIR, "plot_dish_suffix_diff_Overall.png"))
ggsave_safe(ds_early$dish_suffix_pairs,   file.path(OUT_DIR, "plot_dish_suffix_pairs_Early.png"))
ggsave_safe(ds_early$dish_suffix_diff,    file.path(OUT_DIR, "plot_dish_suffix_diff_Early.png"))
ggsave_safe(ds_late$dish_suffix_pairs,    file.path(OUT_DIR, "plot_dish_suffix_pairs_Late.png"))
ggsave_safe(ds_late$dish_suffix_diff,     file.path(OUT_DIR, "plot_dish_suffix_diff_Late.png"))

# Session info to file
sess_log <- capture.output({ cat("\n=== Session info ===\n"); print(sessionInfo()) })
save_txt(sess_log, file.path(OUT_DIR, "session_info.txt"))

cat(sprintf("\nAll outputs saved to: %s\n", OUT_DIR))

