# ============================================
# pupillometry_preprocess_eyelink.R
# Ali Çağan Kaya, 2025
# alicagankaya.com
# ============================================

suppressPackageStartupMessages({
  library(data.table)
  library(zoo)      # na.approx
  library(arrow)    # write_parquet (optional)
})

PARAMS <- list(
  in_dir          = "/Users/path/",          # folder with .txt
  out_dir         = "/Users/path/preprocessedR/",
  pattern         = "\\.txt$",       # your files are TSV .txt
  decimal_comma   = TRUE,            # your numbers use commas
  use_parquet     = TRUE,
  # Columns present in your sample (case insensitive, we normalise names)
  keep_cols       = c(
    "RECORDING_SESSION_LABEL",  # -> Subject
    "TRIAL_INDEX",              # -> Trial
    "AVERAGE_PUPIL_SIZE",       # -> Pupil (numeric, comma decimals)
    "AVERAGE_IN_BLINK",         # -> Blink flag (0/1)
    "TIMESTAMP",                # -> absolute sample time (ms)
    "TRIAL_START_TIME",         # -> trial start (ms)
    "condition", "audiofile"    # -> optional condition metadata
  ),
  rename_map      = c(
    RECORDING_SESSION_LABEL = "Subject",
    TRIAL_INDEX             = "Trial",
    AVERAGE_PUPIL_SIZE      = "Pupil",
    AVERAGE_IN_BLINK        = "InBlink",
    TIMESTAMP               = "Timestamp",
    TRIAL_START_TIME        = "TrialStart",
    condition               = "Condition",
    audiofile               = "Audio"
  ),
  trim_upper_q    = 0.999,
  min_pupil       = 0,
  baseline_ms     = c(-200, 0),   # relative to trial start; will auto-fallback
  analysis_ms     = c(700, 2000), # relative to trial start
  target_hz       = 100,
  max_interp_ms   = 75,
  write_qc_csv    = TRUE
)

to_lower_nospace <- function(x) gsub("[^A-Za-z0-9]+", "_", tolower(x))

read_one <- function(path, p) {
  f <- basename(path)
  message("Reading: ", f)
  # Tab-separated with decimal commas and '.' as NA
  dt <- fread(
    path, sep = "\t", quote = "",
    na.strings = c("NA", "", ".", "NaN"),
    showProgress = FALSE
  )
  # standardise names
  setnames(dt, names(dt), to_lower_nospace(names(dt)))
  
  # keep cols
  keep_std <- to_lower_nospace(p$keep_cols)
  present  <- intersect(keep_std, names(dt))
  if (!length(present)) stop("No expected columns in ", f)
  dt <- dt[, ..present]
  
  # Rename
  rn_from <- names(p$rename_map)
  rn_std  <- to_lower_nospace(rn_from)
  rn_to   <- unname(p$rename_map)
  for (i in seq_along(rn_std)) if (rn_std[i] %in% names(dt)) setnames(dt, rn_std[i], rn_to[i])
  
  # Decimal commas to numerics because ExpBuilder...
  # This took hours to figure out...
  if (p$decimal_comma) {
    numc <- intersect(c("Pupil","Timestamp","TrialStart","Trial"), names(dt))
    for (cc in numc) if (is.character(dt[[cc]])) dt[[cc]] <- as.numeric(sub(",", ".", dt[[cc]], fixed = TRUE))
    # Pupil might still be char due to commas; force numeric
    dt[, Pupil := as.numeric(sub(",", ".", as.character(Pupil), fixed = TRUE))]
    dt[, Timestamp := as.numeric(sub(",", ".", as.character(Timestamp), fixed = TRUE))]
    dt[, TrialStart := as.numeric(sub(",", ".", as.character(TrialStart), fixed = TRUE))]
    dt[, Trial := as.integer(Trial)]
  } else {
    dt[, `:=`(Pupil = as.numeric(Pupil), Timestamp = as.numeric(Timestamp),
              TrialStart = as.numeric(TrialStart), Trial = as.integer(Trial))]
  }
  
  # Factors
  if ("Subject" %in% names(dt))   dt[, Subject := as.factor(Subject)]
  if ("Condition" %in% names(dt)) dt[, Condition := as.factor(Condition)]
  if ("Audio" %in% names(dt))     dt[, Audio := as.factor(Audio)]
  
  # Blink flag (0/1) -> logical
  if ("InBlink" %in% names(dt)) dt[, InBlink := as.integer(InBlink)]
  
  # Basic hygiene
  dt <- dt[is.finite(Pupil) & is.finite(Timestamp) & is.finite(TrialStart) & !is.na(Trial)]
  dt <- dt[Pupil > p$min_pupil]
  
  # Light upper trim
  up <- as.numeric(quantile(dt$Pupil, p$trim_upper_q, na.rm = TRUE))
  if (is.finite(up)) dt <- dt[Pupil < up]
  
  # Trial-relative time from Timestamp and TrialStart
  # (not from IP_START_TIME!)
  setorderv(dt, c("Subject","Trial","Timestamp"))
  dt[, TimeRel := Timestamp - TrialStart]
  
  # Estimate Hz
  dt[, dT := Timestamp - shift(Timestamp), by = .(Subject, Trial)]
  med_dt <- median(dt$dT[is.finite(dt$dT) & dt$dT > 0], na.rm = TRUE)
  hz_est <- ifelse(is.finite(med_dt) && med_dt > 0, round(1000/med_dt), NA_real_)
  dt[, dT := NULL]
  attr(dt, "hz_est") <- hz_est
  attr(dt, "med_dt") <- med_dt
  
  dt
}

interp_baseline <- function(dt, p) {
  # Set blink samples to NA before interpolation (if available)
  if ("InBlink" %in% names(dt)) dt[InBlink == 1, Pupil := NA_real_]
  
  # Interpolate short gaps per trial
  med_dt <- attr(dt, "med_dt"); if (!is.finite(med_dt) || med_dt <= 0) med_dt <- 1 # fallback ~1000 Hz
  max_gap_n <- ceiling(p$max_interp_ms / med_dt)
  
  dt[, PupInterp :=
       { v <- Pupil; if (all(is.na(v))) v else na.approx(v, maxgap = max_gap_n, na.rm = FALSE) },
     by = .(Subject, Trial)]
  
  # Baseline window
  bl <- p$baseline_ms
  # If baseline window is negative but TimeRel starts at 0, fallback to [0, 200] ms
  if (bl[1] < 0 && min(dt$TimeRel, na.rm = TRUE) >= 0) bl <- c(0, 200)
  
  dt[, Base := median(PupInterp[TimeRel >= bl[1] & TimeRel <= bl[2]], na.rm = TRUE),
     by = .(Subject, Trial)]
  dt[, PupilBc := PupInterp - Base]
  
  keep <- intersect(c("Subject","Trial","Timestamp","TrialStart","TimeRel",
                      "Condition","Audio","Pupil","PupInterp","Base","PupilBc"),
                    names(dt))
  dt[, ..keep]
}

downsample_window <- function(dt, p) {
  bin_ms <- 1000 / p$target_hz
  dt[, bin := floor(TimeRel / bin_ms)]
  ds <- dt[, .(
    TimeRel   = mean(TimeRel, na.rm = TRUE),
    Timestamp = mean(Timestamp, na.rm = TRUE),
    Pupil     = mean(Pupil, na.rm = TRUE),
    PupInterp = mean(PupInterp, na.rm = TRUE),
    PupilBc   = mean(PupilBc, na.rm = TRUE),
    Base      = mean(Base, na.rm = TRUE)
  ), by = .(Subject, Trial, bin, Condition, Audio)]
  ds[, bin := NULL]
  
  aw <- p$analysis_ms
  ds <- ds[TimeRel >= aw[1] & TimeRel <= aw[2]]
  ds
}

qc_row <- function(raw_dt, proc_dt, p, file_id) {
  data.table(
    file      = file_id,
    rows_raw  = nrow(raw_dt),
    rows_proc = nrow(proc_dt),
    est_hz    = attr(raw_dt, "hz_est"),
    med_dt_ms = attr(raw_dt, "med_dt"),
    min_rel   = suppressWarnings(min(raw_dt$TimeRel, na.rm = TRUE)),
    max_rel   = suppressWarnings(max(raw_dt$TimeRel, na.rm = TRUE)),
    kept_from = p$analysis_ms[1],
    kept_to   = p$analysis_ms[2],
    target_hz = p$target_hz
  )
}

# ---------- Main ----------
dir.create(PARAMS$out_dir, showWarnings = FALSE, recursive = TRUE)
files <- list.files(PARAMS$in_dir, pattern = PARAMS$pattern, full.names = TRUE, ignore.case = TRUE)
if (!length(files)) stop("No input files found in: ", normalizePath(PARAMS$in_dir))

qc_all <- rbindlist(lapply(files, function(fp) {
  base <- tools::file_path_sans_ext(basename(fp))
  
  dt_raw <- read_one(fp, PARAMS)
  dt_ibc <- interp_baseline(dt_raw, PARAMS)
  ds     <- downsample_window(dt_ibc, PARAMS)
  
  # outputs
  csv_out <- file.path(PARAMS$out_dir, paste0(base, "_clean.csv"))
  fwrite(ds, csv_out)
  if (PARAMS$use_parquet) {
    pq_out <- file.path(PARAMS$out_dir, paste0(base, "_clean.parquet"))
    write_parquet(ds, pq_out)
  }
  message("✓ Wrote: ", csv_out)
  
  qc_row(dt_ibc, ds, PARAMS, basename(fp))
}), fill = TRUE)

if (PARAMS$write_qc_csv) {
  qc_file <- file.path(PARAMS$out_dir, "QC_summary.csv")
  fwrite(qc_all, qc_file)
  message("QC summary → ", qc_file)
}

message("Done. Files processed: ", length(files))
