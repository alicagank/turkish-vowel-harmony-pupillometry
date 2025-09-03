# ==========================================================
# merge_pupil_clean_csvs_keep_all.R
# - Merge *_clean.csv
# - Keep ALL existing columns
# - Add: Participant, HarmonyType, SuffixType
# - Drop fillers (Audio startswith "filler_" or Condition contains "filler")
# - Ensure canonical columns exist by COPYING (no renames): Pupil, Time
# Ali Çağan Kaya, 2025
# alicagankaya.com
# ==========================================================

suppressPackageStartupMessages({
  library(data.table)
  library(stringr)
  library(tools)
})

# ---- paths ----
IN_DIR   <- "/Users/path/"
OUT_FILE <- "/Users/path/file.csv"
dir.create(dirname(OUT_FILE), showWarnings = FALSE, recursive = TRUE)

# ---- options ----
DROP_FILLERS <- TRUE
# Preferred sources to COPY into canonical columns (originals kept):
PUPIL_PREF <- c("Pupil", "PupInterp", "PupilBc")
TIME_PREF  <- c("Time", "Timestamp", "TimeRel")

# ---- helpers ----
rename_if_present_ci <- function(DT, old_ci, new) {
  idx <- which(tolower(names(DT)) == tolower(old_ci))
  if (length(idx)) setnames(DT, idx, new)
}

first_present <- function(cols, choices) {
  hit <- intersect(choices, cols)
  if (length(hit)) hit[1] else NA_character_
}

add_canonical_copy <- function(DT, target, candidates) {
  # If target column is missing, create it by copying the first candidate present
  if (!(target %in% names(DT))) {
    src <- first_present(names(DT), candidates)
    if (!is.na(src)) DT[, (target) := DT[[src]]]
  }
  invisible(NULL)
}

clean_one <- function(fp) {
  dt <- fread(fp, showProgress = FALSE)
  
  # Normalize a few names (case-insensitive) but KEEP all columns
  rename_if_present_ci(dt, "subject",   "Subject")
  rename_if_present_ci(dt, "trial",     "Trial")
  rename_if_present_ci(dt, "condition", "Condition")
  rename_if_present_ci(dt, "audio",     "Audio")
  rename_if_present_ci(dt, "timestamp", "Timestamp")
  rename_if_present_ci(dt, "timerel",   "TimeRel")
  rename_if_present_ci(dt, "pupil",     "Pupil")
  rename_if_present_ci(dt, "pupinterp", "PupInterp")
  rename_if_present_ci(dt, "pupilbc",   "PupilBc")
  
  # Ensure canonical columns exist by COPYING from preferred sources
  add_canonical_copy(dt, "Pupil", PUPIL_PREF)
  add_canonical_copy(dt, "Time",  TIME_PREF)
  
  # Basic types
  if ("Trial" %in% names(dt) && !is.integer(dt$Trial)) dt[, Trial := as.integer(Trial)]
  if ("Condition" %in% names(dt)) dt[, Condition := as.character(Condition)]
  if ("Audio" %in% names(dt))     dt[, Audio := as.character(Audio)]
  
  # Fix occasional typo in Condition
  if ("Condition" %in% names(dt)) {
    dt[Condition == "front_disharmonic_int", Condition := "front_disharmonic_inf"]
  }
  
  # Drop fillers if requested
  if (DROP_FILLERS) {
    if ("Audio" %in% names(dt))     dt <- dt[!startsWith(tolower(Audio), "filler_")]
    if ("Condition" %in% names(dt)) dt <- dt[!grepl("filler", Condition, ignore.case = TRUE)]
  }
  
  # Add Participant from filename (strip _clean + extension)
  base <- file_path_sans_ext(basename(fp))              # e.g., "VH_P02_sample_report_clean"
  participant <- sub("_clean$", "", base, perl = TRUE)  # -> "VH_P02_sample_report"
  dt[, Participant := participant]
  
  # Add HarmonyType & SuffixType (don’t overwrite if they already exist)
  if (!"HarmonyType" %in% names(dt) && "Condition" %in% names(dt)) {
    dt[, HarmonyType := str_match(Condition, "(harmonic|disharmonic)")[, 2]]
  }
  if (!"SuffixType" %in% names(dt) && "Condition" %in% names(dt)) {
    dt[, SuffixType  := str_match(Condition, "(der|inf)")[, 2]]
  }
  
  # Reorder (optional): move the three new columns near Condition/Audio, keep all others
  preferred_front <- c("Subject","Trial","Pupil","Time","Condition","Audio","Participant","HarmonyType","SuffixType")
  existing_front  <- intersect(preferred_front, names(dt))
  the_rest        <- setdiff(names(dt), existing_front)
  setcolorder(dt, c(existing_front, the_rest))
  
  dt
}

# ---- run ----
files <- list.files(IN_DIR, pattern = "_clean\\.csv$", full.names = TRUE)
if (!length(files)) stop("No *_clean.csv found in: ", IN_DIR)

merged <- rbindlist(lapply(files, clean_one), fill = TRUE, use.names = TRUE)

# Optional sort
ord_cols <- intersect(c("Subject","Trial","Time"), names(merged))
if (length(ord_cols)) setorder(merged, !!!rlang::syms(ord_cols)) else setorder(merged)

# Write
fwrite(merged, OUT_FILE)
message("Merged ", format(nrow(merged), big.mark=","), " rows from ", length(files),
        " files → ", OUT_FILE)

# Peek
print(head(merged, 6))
