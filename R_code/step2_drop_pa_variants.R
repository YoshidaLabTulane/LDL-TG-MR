##################################################################
###   Step 2: Drop palindromic high-MAF + low-quality variants  ###
###     Take the three csv files from step 1                  
###    Flags/drops: 
###                 high-MAF palindromic SNPs, (A/T, C/G, with MAF > 0.42)
###                 low INFO/low call-rate (if those columns exist)
###    write cleaned instrument set
###                                                             ###
###    Modified: 12-4-2025                                     ###
##################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
})

## ----------------------------------------------------------
## Paths
## ----------------------------------------------------------
base      <- "/Users/ylu6/Desktop/ADA_2026"

step1_dir <- file.path(base, "output", "step1")
step2_dir <- file.path(base, "output", "step2")
dir.create(step2_dir, showWarnings = FALSE, recursive = TRUE)

## ----------------------------------------------------------
## Flag palindromic high-MAF SNPs
## ----------------------------------------------------------
flag_palindromic <- function(dt, maf_thresh = 0.42) {
  if (!("eaf" %in% names(dt))) {
    stop("Column 'eaf' not found in input instruments; check Step 1 output.")
  }

  dt %>%
    mutate(
      maf = pmin(eaf, 1 - eaf),
      pal = (effect_allele %in% c("A","T") & other_allele %in% c("A","T")) |
            (effect_allele %in% c("C","G") & other_allele %in% c("C","G")),
      drop_pal = pal & maf > maf_thresh
    )
}

## ----------------------------------------------------------
## Apply INFO / call-rate filters if present
## ----------------------------------------------------------
apply_qc_filters <- function(dt,
                             info_min = 0.8,
                             call_min = 0.95) {
  has_info     <- "INFO"     %in% names(dt)
  has_callrate <- "CALLRATE" %in% names(dt)

  dt %>%
    mutate(
      drop_info = if (has_info)     INFO     < info_min  else FALSE,
      drop_call = if (has_callrate) CALLRATE < call_min  else FALSE
    )
}

## ----------------------------------------------------------
## Clean one instrument file and write to step2
## ----------------------------------------------------------
clean_instruments <- function(in_csv,
                              out_csv   = NULL,
                              maf_thresh = 0.42,
                              info_min   = 0.8,
                              call_min   = 0.95) {

  message("Cleaning instruments from: ", in_csv)
  dt <- fread(in_csv)
  n0 <- nrow(dt)

  dt2 <- dt %>%
    flag_palindromic(maf_thresh = maf_thresh) %>%
    apply_qc_filters(info_min = info_min, call_min = call_min)

  dt_keep <- dt2 %>%
    filter(!drop_pal, !drop_info, !drop_call)

  message("Input instruments:   ", n0)
  message("Drop palindromic:    ", sum(dt2$drop_pal))
  if ("drop_info" %in% names(dt2)) {
    message("Drop low INFO:      ", sum(dt2$drop_info))
  }
  if ("drop_call" %in% names(dt2)) {
    message("Drop low call-rate: ", sum(dt2$drop_call))
  }
  message("Final instruments:   ", nrow(dt_keep))

  if (is.null(out_csv)) {
    out_csv <- file.path(
      step2_dir,
      "bolt_stats_TG_LDL_alone_instruments_step2_cleaned.csv"
    )
  }

  fwrite(dt_keep, out_csv)
  message("Wrote cleaned instruments to: ", out_csv, "\n")

  invisible(dt_keep)
}

## ----------------------------------------------------------
## Run Step 2 on TG-LDL-alone instruments ONLY
## ----------------------------------------------------------
tg_step1_file <- file.path(step1_dir, "bolt_stats_TG_LDL_alone_instruments.csv")
if (!file.exists(tg_step1_file)) {
  stop("Step 1 TG_LDL_alone instrument file not found at: ", tg_step1_file)
}

message("Found TG-LDL Step 1 instrument file:")
print(tg_step1_file)

clean_instruments(tg_step1_file)
