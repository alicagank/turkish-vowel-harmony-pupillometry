# ============================
# accuracy_rt_2x2.R
# Ali Çağan Kaya, 2025
# alicagankaya.com
# ============================
suppressPackageStartupMessages({
  library(tidyverse)
  library(stringr)
  library(lme4)
  library(lmerTest)
  library(emmeans)
  library(broom)
  library(broom.mixed)
  library(readr)
  library(tools)
})

# --------- OPTIONS ----------
# Allow df adjustments in emmeans (may slow down on very large n)
emmeans::emm_options(pbkrtest.limit = 1e6, lmerTest.limit = 1e6)

# --------- PATHS (EDIT ME) ----------
IN_DIR  <- "/Users/path"  # folder with participant CSVs
OUT_DIR <- file.path(IN_DIR, "analysisR_focus") # same fil path
PLOTS_DIR <- file.path(OUT_DIR, "plots")
dir.create(OUT_DIR, showWarnings = FALSE, recursive = TRUE)
dir.create(PLOTS_DIR, showWarnings = FALSE, recursive = TRUE)

# Optional: test ONE file only (else NULL to load all CSVs in IN_DIR)
ONE_FILE <- NULL  # e.g., ONE_FILE <- file.path(IN_DIR, "P01_clean_mapped.csv")

# --------- PARAMETERS ----------
RT_MIN <- 200   # ms (trim for RT analyses)
RT_MAX <- 5000  # ms

# --------- HELPERS ----------
norm_keys_to_12 <- function(pressed_key_vec) {
  # Accepts 0/1 or 1/2; returns 1 = harmonic, 2 = disharmonic; else NA
  vals <- sort(unique(na.omit(pressed_key_vec)))
  if (identical(vals, c(0,1))) {
    if_else(pressed_key_vec == 1, 1L, if_else(pressed_key_vec == 0, 2L, NA_integer_))
  } else if (all(vals %in% c(1,2))) {
    if_else(pressed_key_vec %in% c(1,2), pressed_key_vec, NA_integer_)
  } else {
    warning("PressedKey contains unexpected values; setting non {0,1,2} to NA.")
    if_else(pressed_key_vec %in% c(1,2), pressed_key_vec, NA_integer_)
  }
}

parse_design <- function(item) {
  base <- str_remove(item, "^filler_")
  tibble(
    Harmony    = if_else(str_detect(base, "disharmonic"), "disharmonic",
                         if_else(str_detect(base, "harmonic"), "harmonic", NA_character_)),
    SuffixType = case_when(
      str_detect(base, "_der$") ~ "der",
      str_detect(base, "_inf$") ~ "inf",
      TRUE ~ NA_character_
    )
  )
}

load_one <- function(fp) {
  raw <- read_csv(fp, show_col_types = FALSE)
  stopifnot(all(c("Item","PressedKey","RT(ms)") %in% names(raw)))
  if (!"OriginalItem" %in% names(raw)) raw$OriginalItem <- raw$Item
  pid <- str_extract(basename(fp), "^[^_]+")
  design <- parse_design(raw$Item)
  
  df <- bind_cols(raw, design) %>%
    mutate(
      Participant    = pid,
      is_filler      = str_detect(Item, "^filler_"),
      PressedKeyNorm = norm_keys_to_12(PressedKey),
      # 1/0 elements as requested:
      PressedHarmonic = if_else(PressedKeyNorm == 1L, 1L,
                                if_else(PressedKeyNorm == 2L, 0L, NA_integer_)),  # 1 = pressed 'harmonic', 0 = 'disharmonic'
      Response       = if_else(PressedKeyNorm == 1L, "harmonic",
                               if_else(PressedKeyNorm == 2L, "disharmonic", NA_character_)),
      Correct        = as.integer(Response == Harmony),   # 1 = correct, 0 = incorrect
      ValidRT        = dplyr::between(`RT(ms)`, RT_MIN, RT_MAX)
    ) %>%
    filter(!is_filler, !is.na(Harmony), !is.na(SuffixType))
  df
}

# --------- LOAD DATA ----------
files <- if (!is.null(ONE_FILE)) ONE_FILE else list.files(IN_DIR, pattern = "\\.csv$", full.names = TRUE)
stopifnot(length(files) > 0)
dat <- purrr::map_dfr(files, load_one)

message("Loaded N = ", nrow(dat), " trials from ", dplyr::n_distinct(dat$Participant), " participant(s).")

# --------- SAVE BASIC QA ----------
dat %>% count(Harmony, SuffixType, Response) %>% arrange(Harmony, SuffixType, Response) %>%
  write_csv(file.path(OUT_DIR, "qa_confusion_counts.csv"))
dat %>% count(Participant, Harmony, SuffixType) %>%
  pivot_wider(names_from = c(Harmony, SuffixType), values_from = n, values_fill = 0) %>%
  write_csv(file.path(OUT_DIR, "qa_cell_sizes_by_participant.csv"))

# ============================================================
# PART 1 — ACCURACY (0/1) : Understanding harmony vs disharmony,
# and disharmonic_der vs disharmonic_inf
# ============================================================

# ---- GLMM: Correct ~ Harmony * SuffixType ----
glmm_acc_full <- glmer(
  Correct ~ Harmony * SuffixType + (1 | Participant) + (1 | OriginalItem),
  data = dat, family = binomial(),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
saveRDS(glmm_acc_full, file.path(OUT_DIR, "glmm_accuracy_full.rds"))
# If item variance is ~0, refit without item RE
# Burası önemli
if (isSingular(glmm_acc_full, tol = 1e-4)) {
  message("Accuracy GLMM is singular; refitting without item RE...")
  glmm_acc <- glmer(
    Correct ~ Harmony * SuffixType + (1 | Participant),
    data = dat, family = binomial(),
    control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
  )
  saveRDS(glmm_acc, file.path(OUT_DIR, "glmm_accuracy_no_itemRE.rds"))
} else {
  glmm_acc <- glmm_acc_full
}

# EMMs on probability (understanding harmonic vs disharmonic, and the 2x2 cell means)
emm_acc_cells <- emmeans(glmm_acc, ~ Harmony * SuffixType, type = "response")  # prob column
acc_cells_df <- as.data.frame(emm_acc_cells)
write_csv(acc_cells_df, file.path(OUT_DIR, "acc_emmeans_cells_prob.csv"))

# Targeted contrast: Disharmonic only, derivatives vs inflectional
emm_acc_by_h <- emmeans(glmm_acc, ~ SuffixType | Harmony, type = "response")
acc_dish_contrast <- contrast(emm_acc_by_h, method = "revpairwise") %>%
  as.data.frame() %>% filter(Harmony == "disharmonic")
write_csv(acc_dish_contrast, file.path(OUT_DIR, "acc_contrast_disharmonic_inf_vs_der.csv"))

# Optional: direct "pressed 1 vs 0" understanding check
# (probability of pressing '1' differs by Harmony)
glmm_press <- glmer(
  PressedHarmonic ~ Harmony + (1 | Participant) + (1 | OriginalItem),
  data = dat, family = binomial(),
  control = glmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 2e5))
)
saveRDS(glmm_press, file.path(OUT_DIR, "glmm_pressed1_by_harmony.rds"))
emm_press <- emmeans(glmm_press, ~ Harmony, type = "response")
write_csv(as.data.frame(emm_press), file.path(OUT_DIR, "pressed1_emmeans_by_harmony.csv"))

# ---- Accuracy summary for plots ----
acc_by_cond <- dat %>%
  group_by(Participant, Harmony, SuffixType) %>%
  summarise(acc = mean(Correct, na.rm = TRUE), n = dplyr::n(), .groups = "drop")
write_csv(acc_by_cond, file.path(OUT_DIR, "acc_by_condition_participant.csv"))

# ============================================================
# PART 2 — REACTION TIME (ms) : Harmony cognitive load? and disharmonic der vs inf
# ============================================================

# Use correct trials within RT window
rt_df <- dat %>% filter(Correct == 1, ValidRT) %>%
  mutate(logRT = log(`RT(ms)`))

# ---- LMM: logRT ~ Harmony * SuffixType ----
lmm_rt <- lmer(
  logRT ~ Harmony * SuffixType + (1 | Participant) + (1 | OriginalItem),
  data = rt_df,
  control = lmerControl(check.nobs.vs.nRE = "ignore")
)
saveRDS(lmm_rt, file.path(OUT_DIR, "lmm_rt_full.rds"))

# Cell means (EMMs) back-transformed to ms
emm_rt_cells <- emmeans(lmm_rt, ~ Harmony * SuffixType)
emm_rt_cells_df <- broom::tidy(emm_rt_cells) %>%
  mutate(emmean_ms = exp(estimate),
         lower_ms  = exp(estimate - 1.96*std.error),
         upper_ms  = exp(estimate + 1.96*std.error)) %>%
  select(Harmony, SuffixType, emmean_ms, lower_ms, upper_ms)
write_csv(emm_rt_cells_df, file.path(OUT_DIR, "rt_emmeans_cells_ms.csv"))

# Targeted contrast: within disharmony, der vs inf (back-transformed)
emm_rt_by_h <- emmeans(lmm_rt, ~ SuffixType | Harmony)
rt_dish_contrast <- contrast(emm_rt_by_h, method = "revpairwise") %>%
  summary(type = "response") %>%  # back-transform from log
  as.data.frame() %>%
  filter(Harmony == "disharmonic")
write_csv(rt_dish_contrast, file.path(OUT_DIR, "rt_contrast_disharmonic_inf_vs_der_ms.csv"))

# Descriptive participant-level effect sizes (paired tests; optional)
# Accuracy: participant-wise mean accuracy by Harmony (collapsing SuffixType)
acc_part_h <- dat %>%
  group_by(Participant, Harmony) %>%
  summarise(acc = mean(Correct, na.rm = TRUE), .groups = "drop") %>%
  pivot_wider(names_from = Harmony, values_from = acc)
if (all(c("harmonic","disharmonic") %in% names(acc_part_h))) {
  acc_t <- t.test(acc_part_h$harmonic, acc_part_h$disharmonic, paired = TRUE)
  capture.output(acc_t, file = file.path(OUT_DIR, "acc_paired_t_harmonic_vs_disharmonic.txt"))
}

# RT: participant-wise mean RT by Harmony (correct trials)
rt_part_h <- rt_df %>%
  group_by(Participant, Harmony) %>%
  summarise(mean_rt = mean(`RT(ms)`), .groups = "drop") %>%
  pivot_wider(names_from = Harmony, values_from = mean_rt)
if (all(c("harmonic","disharmonic") %in% names(rt_part_h))) {
  rt_t <- t.test(rt_part_h$harmonic, rt_part_h$disharmonic, paired = TRUE)
  capture.output(rt_t, file = file.path(OUT_DIR, "rt_paired_t_harmonic_vs_disharmonic.txt"))
}

# ============================================================
# PLOTS
# ============================================================

# ---- Accuracy plot: Harmony × SuffixType (per-participant dots + mean ± SE)
p_acc <- acc_by_cond %>%
  ggplot(aes(x = Harmony, y = acc, colour = SuffixType, group = SuffixType)) +
  geom_point(position = position_jitter(width = .06, height = 0), alpha = .5) +
  stat_summary(fun = mean, geom = "point", size = 3,
               position = position_dodge(width = .4)) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .2,
               position = position_dodge(width = .4)) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(title = "Accuracy by Harmony × SuffixType",
       x = "Harmony", y = "Accuracy", colour = "Suffix type") +
  theme_minimal(base_size = 13)
ggsave(file.path(PLOTS_DIR, "accuracy_by_condition.png"), p_acc, width = 7, height = 4.5, dpi = 300)

# ---- Accuracy (disharmonic only): der vs inf
acc_dish_part <- dat %>%
  filter(Harmony == "disharmonic") %>%
  group_by(Participant, SuffixType) %>%
  summarise(acc = mean(Correct), .groups = "drop")
p_acc_dish <- acc_dish_part %>%
  ggplot(aes(x = SuffixType, y = acc)) +
  geom_point(position = position_jitter(width = .05), alpha = .5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .15) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1), limits = c(0,1)) +
  labs(title = "Accuracy within DISHARMONIC trials: der vs inf",
       x = "Suffix type", y = "Accuracy") +
  theme_minimal(base_size = 13)
ggsave(file.path(PLOTS_DIR, "accuracy_disharmonic_der_vs_inf.png"),
       p_acc_dish, width = 6, height = 4.2, dpi = 300)

# ---- RT plot: Harmony × SuffixType (correct trials)
p_rt <- rt_df %>%
  ggplot(aes(x = interaction(Harmony, SuffixType, sep = " × "), y = `RT(ms)`)) +
  geom_violin(trim = TRUE, alpha = .4) +
  geom_boxplot(width = .15, outlier.shape = NA) +
  geom_jitter(width = .05, alpha = .2, size = .8) +
  labs(title = "RT (ms) for correct trials by Harmony × SuffixType",
       x = "Condition", y = "RT (ms)") +
  theme_minimal(base_size = 13)
ggsave(file.path(PLOTS_DIR, "rt_by_condition.png"), p_rt, width = 7.5, height = 4.5, dpi = 300)

# ---- RT (disharmonic only): der vs inf
rt_dish_part <- rt_df %>%
  filter(Harmony == "disharmonic") %>%
  group_by(Participant, SuffixType) %>%
  summarise(mean_rt = mean(`RT(ms)`), .groups = "drop")
p_rt_dish <- rt_dish_part %>%
  ggplot(aes(x = SuffixType, y = mean_rt)) +
  geom_point(position = position_jitter(width = .05), alpha = .5) +
  stat_summary(fun = mean, geom = "point", size = 3) +
  stat_summary(fun.data = mean_se, geom = "errorbar", width = .15) +
  labs(title = "RT within DISHARMONIC trials: der vs inf (participant means)",
       x = "Suffix type", y = "Mean RT (ms)") +
  theme_minimal(base_size = 13)
ggsave(file.path(PLOTS_DIR, "rt_disharmonic_der_vs_inf.png"),
       p_rt_dish, width = 6, height = 4.2, dpi = 300)

# ============================================================
# CONSOLE SUMMARY (quick read-out)
# ============================================================
message("\n=== GLMM Accuracy (Correct ~ Harmony * SuffixType) ===")
print(summary(glmm_acc))
message("\nEMMs (probability) by Harmony × SuffixType:")
print(emm_acc_cells)

message("\nDisharmonic contrast (Accuracy, inf vs der):")
print(acc_dish_contrast)

message("\n=== LMM logRT (correct, trimmed) ===")
print(summary(lmm_rt))
message("\nEMMs (ms) by Harmony × SuffixType (back-transformed):")
print(emm_rt_cells_df)

message("\nDisharmonic contrast (RT, inf vs der; ms, back-transformed):")
print(rt_dish_contrast)

message("\nPressed '1' (harmonic) probability by Harmony (sanity check):")
print(emm_press)

