####################################################################################
### Step 4 : harmonization for GCST90132314.tsv general cad outcome             ####
###  Hugh: 12-4-2025                                                            ####
####################################################################################


suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readr)
  library(TwoSampleMR)
})

## --------------------------------------------------------------
## 0. Paths
## --------------------------------------------------------------
base <- "/Users/ylu6/Desktop/ADA_2026"

step3_dir <- file.path(base, "output", "step3")
step4_dir <- file.path(base, "output", "step4")
dir.create(step4_dir, showWarnings = FALSE, recursive = TRUE)

# Step 3 exposure file (TG-LDL only)
exp_file <- file.path(step3_dir, "instruments_TG_LDL_alone_with_betas_step3.csv")
stopifnot(file.exists(exp_file))

# CAD outcome GWAS (new filename)
outcome_file <- file.path(base, "CAD", "GCST90132314.tsv")
stopifnot(file.exists(outcome_file))

# Sample sizes (update if you have exact N_used)
N_TG_LDL <- 272051L
N_CAD    <- 184305L   # placeholder from before

message("Exposure file (Step 3): ", exp_file)
message("CAD GWAS file        : ", outcome_file)

## --------------------------------------------------------------
## 1. CAD outcome -> TwoSampleMR format
##    Use rsid as SNP, fallback to markername if rsid is missing
## --------------------------------------------------------------
format_outcome_CAD <- function(file_path, trait_label, N) {
  o <- fread(file_path) %>% as.data.frame()

  if (!all(c("effect_allele", "other_allele", "beta",
             "standard_error", "p_value",
             "effect_allele_frequency") %in% names(o))) {
    stop("Outcome file is missing one or more required columns.")
  }
  if (!("rsid" %in% names(o)) && !("markername" %in% names(o))) {
    stop("Outcome file must have rsid or markername for SNP IDs.")
  }

  o <- o %>%
    mutate(
      SNP = dplyr::case_when(
        "rsid" %in% names(o) & !is.na(rsid) & rsid != "" ~ rsid,
        "markername" %in% names(o)                      ~ markername,
        TRUE                                            ~ NA_character_
      )
    )

  # Drop rows with missing SNP IDs
  o2 <- o %>% filter(!is.na(SNP) & SNP != "")
  message("Outcome rows with non-missing SNP IDs: ", nrow(o2))

  out <- TwoSampleMR::format_data(
    o2 %>% transmute(
      SNP,
      effect_allele  = effect_allele,
      other_allele   = other_allele,
      beta           = beta,
      se             = standard_error,
      pval           = p_value,
      eaf            = effect_allele_frequency
    ),
    type              = "outcome",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    phenotype_col     = trait_label,
    samplesize        = N
  )

  out$outcome <- trait_label
  out
}

outcome_dat <- format_outcome_CAD(outcome_file, "CAD", N_CAD)
message("Formatted CAD outcome SNPs: ", nrow(outcome_dat))

## --------------------------------------------------------------
## 2. Step 3 TG-LDL exposure -> TwoSampleMR exposure format
## --------------------------------------------------------------
format_exposure_from_step3 <- function(step3_file, label, N_exp) {
  dt <- fread(step3_file)

  needed <- c("SNP", "effect_allele", "other_allele", "beta", "se")
  miss   <- setdiff(needed, names(dt))
  if (length(miss)) {
    stop("Step 3 file missing columns: ", paste(miss, collapse = ", "),
         " in ", step3_file)
  }

  df <- dt %>%
    transmute(
      SNP,
      effect_allele,
      other_allele,
      beta,
      se,
      eaf  = if ("eaf" %in% names(dt)) eaf else NA_real_,
      pval = NA_real_
    )

  out <- TwoSampleMR::format_data(
    as.data.frame(df),
    type              = "exposure",
    snp_col           = "SNP",
    beta_col          = "beta",
    se_col            = "se",
    effect_allele_col = "effect_allele",
    other_allele_col  = "other_allele",
    eaf_col           = "eaf",
    pval_col          = "pval",
    phenotype_col     = label,
    samplesize        = N_exp
  )

  out$exposure <- label
  out
}

exp_TG <- format_exposure_from_step3(
  step3_file = exp_file,
  label      = "TG_LDL_alone",
  N_exp      = N_TG_LDL
)

message("Formatted TG-LDL exposure SNPs: ", nrow(exp_TG))

## --------------------------------------------------------------
## 3. Harmonise TG-LDL exposure with CAD outcome
## --------------------------------------------------------------
harm <- suppressMessages(
  TwoSampleMR::harmonise_data(
    exposure_dat = exp_TG,
    outcome_dat  = outcome_dat,
    action       = 2
  )
)

message("Harmonised SNPs (TG_LDL_alone vs CAD): ", nrow(harm))

## --------------------------------------------------------------
## 4. Save harmonised dataset
## --------------------------------------------------------------
harm_file <- file.path(step4_dir, "harm_TG_LDL_alone_CAD_step4.csv")
fwrite(harm, harm_file)
message("Wrote harmonised data to: ", harm_file)
