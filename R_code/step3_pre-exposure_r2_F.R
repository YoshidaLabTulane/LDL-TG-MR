####################################################################
## Uses the Step 2 cleaned instrument files for each exposure:   ###
##   bolt_stats_TG_LDL_alone_instruments_step2_cleaned.csv       ###
## Reads the corresponding GWAS files:                           ###
##   bolt.stats.TG_LDL_alone.gz                                  ###
## Computes per-exposure strength metrics (K, R²_total,mean F,   ###
##   min F, overall F)                                           ### 
##   Tweak N if BOLT logs differ                                 ###
##    Update: 12-5-2025                                          ###
####################################################################


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readr)
})

# ------------------------------------------------------------------
# 0. Paths and output
# ------------------------------------------------------------------
base        <- "/Users/ylu6/Desktop/ADA_2026"

step2_dir   <- file.path(base, "output", "step2")
step3_dir  <- file.path(base, "output", "step3")
dir.create(step3_dir, showWarnings = FALSE, recursive = TRUE)

gwas_dir <- file.path(base, "resource_and_gwas_data", "UKB_inputs", "TG_LDL")

# ------------------------------------------------------------------
# 1. Exposure specification (TG-LDL only)
#    - inst_file: cleaned instrument set from Step 2
#    - gwas_file: GWAS summary file to pull beta / SE / EAF
#    - N        : sample size used in BOLT
# ------------------------------------------------------------------
exposure_specs <- tribble(
  ~exposure,      ~inst_file,                                                                 ~gwas_file,                           ~N,
  "TG_LDL_alone",
    file.path(step2_dir, "bolt_stats_TG_LDL_alone_instruments_step2_cleaned.csv"),
    file.path(gwas_dir,  "bolt.stats.TG_LDL_alone.gz"),
    272051L    # update if N_used in BOLT logs is different
)

print(exposure_specs)

# ------------------------------------------------------------------
# 2. Function to compute strength metrics for one exposure
# ------------------------------------------------------------------
compute_strength_metrics <- function(dt, N, exposure_name = "exposure") {
  dt2 <- dt %>%
    mutate(
      maf   = pmin(eaf, 1 - eaf),
      F_snp = (beta / se)^2,
      R2_snp = (2 * beta^2 * maf * (1 - maf)) /
               (2 * beta^2 * maf * (1 - maf) + se^2 * N)
    )

  K     <- nrow(dt2)
  R2tot <- sum(dt2$R2_snp, na.rm = TRUE)
  Fbar  <- mean(dt2$F_snp, na.rm = TRUE)
  Fmin  <- min(dt2$F_snp, na.rm = TRUE)

  F_overall <- ((N - K - 1) / K) * (R2tot / (1 - R2tot))

  tibble(
    exposure   = exposure_name,
    N          = N,
    K          = K,
    R2_total   = R2tot,
    F_mean     = Fbar,
    F_min      = Fmin,
    F_overall  = F_overall
  )
}

# ------------------------------------------------------------------
# 3. Helper: load exposure GWAS, slice to instrument SNPs, join back
#      *with explicit de-duplication per SNP*
# ------------------------------------------------------------------
load_exposure_for_instruments <- function(exposure_name,
                                          inst_file,
                                          gwas_file) {
  if (!file.exists(inst_file)) {
    stop("Instrument file for ", exposure_name, " not found at: ", inst_file)
  }
  if (!file.exists(gwas_file)) {
    stop("GWAS file for ", exposure_name, " not found at: ", gwas_file)
  }

  message("--------------------------------------------------")
  message("Exposure: ", exposure_name)
  message("  Instruments: ", inst_file)
  message("  GWAS file : ", gwas_file)

  inst <- fread(inst_file)
  if (!("SNP" %in% names(inst))) {
    stop("Instrument file for ", exposure_name, " must contain a 'SNP' column.")
  }

  # De-duplicate instruments by SNP (keep first occurrence)
  n_inst0 <- nrow(inst)
  inst <- inst %>% distinct(SNP, .keep_all = TRUE)
  if (nrow(inst) < n_inst0) {
    message("  [Note] Dropped ", n_inst0 - nrow(inst),
            " duplicate instrument rows (same SNP).")
  }

  snps <- inst$SNP
  message("  Unique instruments from Step 2: ", length(snps))

  gwas <- fread(gwas_file)

  needed_cols <- c("SNP", "BETA", "SE", "A1FREQ")
  missing_cols <- setdiff(needed_cols, names(gwas))
  if (length(missing_cols)) {
    stop("GWAS file for ", exposure_name,
         " is missing columns: ", paste(missing_cols, collapse = ", "))
  }

  gwas_sub <- gwas[SNP %in% snps,
                   .(SNP,
                     eaf  = A1FREQ,
                     beta = BETA,
                     se   = SE)]

  if (nrow(gwas_sub) == 0L) {
    stop("No overlapping SNPs between instruments and GWAS for ", exposure_name)
  }

  # De-duplicate GWAS slice by SNP as well
  n_g0 <- nrow(gwas_sub)
  gwas_sub <- gwas_sub %>% distinct(SNP, .keep_all = TRUE)
  if (nrow(gwas_sub) < n_g0) {
    message("  [Note] Dropped ", n_g0 - nrow(gwas_sub),
            " duplicate GWAS rows (same SNP).")
  }

  message("  Overlap with instruments (unique SNPs): ", nrow(gwas_sub), " SNPs")

  core_cols <- intersect(c("SNP", "CHR", "BP", "effect_allele", "other_allele"),
                         names(inst))

  inst_core <- inst[, ..core_cols]

  exp_inst <- inst_core %>%
    inner_join(gwas_sub, by = "SNP")

  # Final sanity check: one row per SNP
  if (any(duplicated(exp_inst$SNP))) {
    warning("After join, still found duplicated SNPs in exp_inst. ",
            "Consider inspecting inst_file and gwas_file directly.")
  } else {
    message("  Final instrument table rows: ", nrow(exp_inst),
            " (all unique SNPs)")
  }

  exp_inst
}

# ------------------------------------------------------------------
# 4. Build instrument table for TG-LDL and compute F/R²
# ------------------------------------------------------------------
strength_list <- list()

for (i in seq_len(nrow(exposure_specs))) {
  exposure_name <- exposure_specs$exposure[i]
  inst_file     <- exposure_specs$inst_file[i]
  gwas_file     <- exposure_specs$gwas_file[i]
  N_exposure    <- exposure_specs$N[i]

  exp_inst <- load_exposure_for_instruments(
    exposure_name = exposure_name,
    inst_file     = inst_file,
    gwas_file     = gwas_file
  )

  # Save exposure-specific instrument table with betas
  out_csv <- file.path(
    step3_dir,
    paste0("instruments_", exposure_name, "_with_betas_step3.csv")
  )
  fwrite(exp_inst, out_csv)
  message("  Wrote instruments+betas to: ", out_csv)

  # Compute strength metrics for this exposure
  metrics <- compute_strength_metrics(
    dt            = exp_inst,
    N             = N_exposure,
    exposure_name = exposure_name
  )

  print(metrics)
  strength_list[[exposure_name]] <- metrics
}

strength_summary <- bind_rows(strength_list)

message("--------------------------------------------------")
message("TG-LDL single-exposure instrument strength summary (Step 3):")
print(strength_summary)

# ------------------------------------------------------------------
# 5. Save summary to disk
# ------------------------------------------------------------------
out_csv <- file.path(step3_dir, "instrument_strength_step3_TG_LDL_alone.csv")
write_csv(strength_summary, out_csv)
message("Wrote Step 3 strength summary to: ", out_csv)
