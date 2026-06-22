##################################################################
## Step 5 (CAD): Univariable MR  (TwoSampleMR)
## TG_LDL_alone → CAD
##
##   Uses harmonised SNPs from Step 4:
##     output/step4/harm_TG_LDL_alone_CAD_step4.csv
##
##   Outputs (in output/step5):
##     - univariable_MR_TG_LDL_alone_on_CAD_step5.csv
##     - Egger_intercept_TG_LDL_alone_on_CAD_step5.csv
##     - harm_TG_LDL_alone_CAD_step5.csv  (copy of harmonised data)
##
##   Updated 12.4.2025 – now includes exposure & outcome Ns
##################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readr)
  library(TwoSampleMR)
})

## --------------------------------------------------------------
## 0. Paths & sample sizes
## --------------------------------------------------------------
base <- "/Users/ylu6/Desktop/ADA_2026"

# Step 4 directory: harmonised SNPs
step4_dir <- file.path(base, "output", "step4")

# Step 5 MR output directory
step5_dir <- file.path(base, "output", "step5")
dir.create(step5_dir, showWarnings = FALSE, recursive = TRUE)

# Input: harmonised TG_LDL_alone → CAD from Step 4
harm_in_file <- file.path(step4_dir, "harm_TG_LDL_alone_CAD_step4.csv")
stopifnot(file.exists(harm_in_file))

# Exposure sample size (TG-LDL GWAS in UKB)
N_TG_LDL <- 272051L

# CAD meta-analysis sample sizes
# From the multi-consortium GWAS (GCST90132314):
# 181,522 CAD cases among 1,165,690 participants
N_CAD_cases    <- 181522L
N_CAD_total    <- 1165690L
N_CAD_controls <- N_CAD_total - N_CAD_cases

message("Reading harmonised CAD data from Step 4: ", harm_in_file)
message("Exposure N (TG_LDL): ", N_TG_LDL)
message("CAD GWAS N (cases / controls / total): ",
        N_CAD_cases, " / ", N_CAD_controls, " / ", N_CAD_total)

harm <- fread(harm_in_file) %>% as.data.frame()
message("Harmonised SNPs (Step 4): ", nrow(harm))

# Optional: quick sanity checks on expected columns
expected_cols <- c(
  "SNP", "beta.exposure", "se.exposure",
  "beta.outcome", "se.outcome",
  "effect_allele.exposure", "other_allele.exposure",
  "effect_allele.outcome", "other_allele.outcome",
  "exposure", "outcome"
)

missing_cols <- setdiff(expected_cols, names(harm))
if (length(missing_cols) > 0) {
  warning("The harmonised CAD file is missing some typical harmonise_data columns: ",
          paste(missing_cols, collapse = ", "),
          "\nAs long as TwoSampleMR::mr() works, this is not fatal.")
}

label_exp <- unique(harm$exposure)
label_out <- unique(harm$outcome)

message("--------------------------------------------------")
message("Running univariable MR for exposure: ", paste(label_exp, collapse = ", "))
message("Outcome: ", paste(label_out, collapse = ", "))

## --------------------------------------------------------------
## 1. Run MR (IVW, Egger, weighted median)
## --------------------------------------------------------------
mr_res <- TwoSampleMR::mr(
  harm,
  method_list = c(
    "mr_ivw",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)

# Add labels
mr_res$exposure_label <- if (length(label_exp) == 1) label_exp else NA_character_
mr_res$outcome_label  <- if (length(label_out) == 1) label_out else NA_character_

# Add OR + 95% CI + Ns
mr_res <- mr_res %>%
  as_tibble() %>%
  mutate(
    OR                 = exp(b),
    OR_LCI95           = exp(b - 1.96 * se),
    OR_UCI95           = exp(b + 1.96 * se),
    N_exposure         = N_TG_LDL,
    N_outcome_total    = N_CAD_total,
    N_outcome_cases    = N_CAD_cases,
    N_outcome_controls = N_CAD_controls
  )

print(mr_res)

## --------------------------------------------------------------
## 2. Egger intercept / directional pleiotropy test
## --------------------------------------------------------------
egger_res <- tryCatch(
  {
    er <- TwoSampleMR::mr_pleiotropy_test(harm) %>%
      as_tibble() %>%
      mutate(
        exposure_label     = if (length(label_exp) == 1) label_exp else NA_character_,
        outcome_label      = if (length(label_out) == 1) label_out else NA_character_,
        N_exposure         = N_TG_LDL,
        N_outcome_total    = N_CAD_total,
        N_outcome_cases    = N_CAD_cases,
        N_outcome_controls = N_CAD_controls
      )
    message("Egger intercept test completed for ", paste(label_exp, collapse = ", "))
    print(er)
    er
  },
  error = function(e) {
    message("mr_pleiotropy_test failed for ", paste(label_exp, collapse = ", "),
            ": ", e$message)
    tibble()
  }
)

## --------------------------------------------------------------
## 3. Save outputs
## --------------------------------------------------------------
# (a) MR results
out_mr_file <- file.path(step5_dir, "univariable_MR_TG_LDL_alone_on_CAD_step5.csv")
write_csv(mr_res, out_mr_file)
message("Wrote univariable MR summary (CAD) to: ", out_mr_file)

# (b) Egger intercept
if (nrow(egger_res) > 0) {
  pleio_file <- file.path(step5_dir, "Egger_intercept_TG_LDL_alone_on_CAD_step5.csv")
  write_csv(egger_res, pleio_file)
  message("Wrote Egger intercept (pleiotropy, CAD) results to: ", pleio_file)
} else {
  message("No Egger intercept results were written (empty table).")
}

# (c) Copy the harmonised file into step5 for convenience
harm_copy <- file.path(step5_dir, "harm_TG_LDL_alone_CAD_step5.csv")
fwrite(harm, harm_copy)
message("Copied harmonised CAD data from Step 4 to: ", harm_copy)

message("Step 5 (CAD two-sample MR with Egger intercept) complete.")


#########################################################################





































##################################################################
## Step 5 (CAD): Univariable MR  (TwoSampleMR)
## TG_LDL_alone → CAD
##
##   Uses harmonised SNPs from Step 4:
##     output/step4/harm_TG_LDL_alone_CAD_step4.csv
##
##   Outputs (in output/step5):
##     - univariable_MR_TG_LDL_alone_on_CAD_step5.csv
##     - Heterogeneity_Cochran_Q_TG_LDL_alone_on_CAD_step5.csv
##     - Egger_intercept_TG_LDL_alone_on_CAD_step5.csv
##     - harm_TG_LDL_alone_CAD_step5.csv
##
##   Updated: includes MR, Cochran's Q heterogeneity test,
##            and MR-Egger intercept pleiotropy test
##################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(readr)
  library(TwoSampleMR)
})

## --------------------------------------------------------------
## 0. Paths & sample sizes
## --------------------------------------------------------------
base <- "/Users/ylu6/Desktop/ADA_2026"

# Step 4 directory: harmonised SNPs
step4_dir <- file.path(base, "output", "step4")

# Step 5 MR output directory
step5_dir <- file.path(base, "output", "step5")
dir.create(step5_dir, showWarnings = FALSE, recursive = TRUE)

# Input: harmonised TG_LDL_alone → CAD from Step 4
harm_in_file <- file.path(step4_dir, "harm_TG_LDL_alone_CAD_step4.csv")
stopifnot(file.exists(harm_in_file))

# Exposure sample size: TG-LDL GWAS in UKB
N_TG_LDL <- 272051L

# CAD meta-analysis sample sizes
# From the multi-consortium GWAS GCST90132314:
# 181,522 CAD cases among 1,165,690 participants
N_CAD_cases    <- 181522L
N_CAD_total    <- 1165690L
N_CAD_controls <- N_CAD_total - N_CAD_cases

message("Reading harmonised CAD data from Step 4: ", harm_in_file)
message("Exposure N (TG_LDL): ", N_TG_LDL)
message("CAD GWAS N (cases / controls / total): ",
        N_CAD_cases, " / ", N_CAD_controls, " / ", N_CAD_total)

harm <- fread(harm_in_file) %>% as.data.frame()
message("Harmonised SNPs from Step 4: ", nrow(harm))

## --------------------------------------------------------------
## 0.1 Sanity checks
## --------------------------------------------------------------
expected_cols <- c(
  "SNP",
  "beta.exposure", "se.exposure",
  "beta.outcome", "se.outcome",
  "effect_allele.exposure", "other_allele.exposure",
  "effect_allele.outcome", "other_allele.outcome",
  "exposure", "outcome"
)

missing_cols <- setdiff(expected_cols, names(harm))

if (length(missing_cols) > 0) {
  warning(
    "The harmonised CAD file is missing some typical harmonise_data columns: ",
    paste(missing_cols, collapse = ", "),
    "\nAs long as TwoSampleMR::mr(), mr_heterogeneity(), and ",
    "mr_pleiotropy_test() work, this is not fatal."
  )
}

label_exp <- unique(harm$exposure)
label_out <- unique(harm$outcome)

message("--------------------------------------------------")
message("Running univariable MR for exposure: ", paste(label_exp, collapse = ", "))
message("Outcome: ", paste(label_out, collapse = ", "))

## Helper function to add labels and sample sizes
add_study_info <- function(df) {
  df %>%
    as_tibble() %>%
    mutate(
      exposure_label     = if (length(label_exp) == 1) label_exp else NA_character_,
      outcome_label      = if (length(label_out) == 1) label_out else NA_character_,
      N_exposure         = N_TG_LDL,
      N_outcome_total    = N_CAD_total,
      N_outcome_cases    = N_CAD_cases,
      N_outcome_controls = N_CAD_controls
    )
}

## --------------------------------------------------------------
## 1. Run MR: IVW, MR-Egger, weighted median
## --------------------------------------------------------------
mr_res <- TwoSampleMR::mr(
  harm,
  method_list = c(
    "mr_ivw",
    "mr_egger_regression",
    "mr_weighted_median"
  )
)

mr_res <- mr_res %>%
  as_tibble() %>%
  mutate(
    exposure_label     = if (length(label_exp) == 1) label_exp else NA_character_,
    outcome_label      = if (length(label_out) == 1) label_out else NA_character_,
    
    # For binary CAD outcome, exponentiate log-OR scale beta
    OR                 = exp(b),
    OR_LCI95           = exp(b - 1.96 * se),
    OR_UCI95           = exp(b + 1.96 * se),
    
    N_exposure         = N_TG_LDL,
    N_outcome_total    = N_CAD_total,
    N_outcome_cases    = N_CAD_cases,
    N_outcome_controls = N_CAD_controls
  )

message("MR results:")
print(mr_res)

## --------------------------------------------------------------
## 2. Cochran's Q heterogeneity test
##
## Important:
##   - Run this once.
##   - TwoSampleMR will return heterogeneity results for methods that support it.
##   - Usually you report:
##       * IVW Cochran's Q as the main heterogeneity test
##       * MR-Egger Q as a sensitivity heterogeneity test
##   - Weighted median does not usually have a Cochran's Q row.
## --------------------------------------------------------------
het_res <- tryCatch(
  {
    hr <- TwoSampleMR::mr_heterogeneity(
      harm,
      method_list = c(
        "mr_ivw",
        "mr_egger_regression"
      )
    ) %>%
      add_study_info()
    
    message("Cochran's Q heterogeneity test completed.")
    print(hr)
    hr
  },
  error = function(e) {
    message("mr_heterogeneity failed for ",
            paste(label_exp, collapse = ", "),
            " → ",
            paste(label_out, collapse = ", "),
            ": ",
            e$message)
    tibble()
  }
)

## Optional clean interpretation column for heterogeneity
if (nrow(het_res) > 0 && "Q_pval" %in% names(het_res)) {
  het_res <- het_res %>%
    mutate(
      heterogeneity_interpretation = case_when(
        is.na(Q_pval) ~ NA_character_,
        Q_pval < 0.05 ~ "Evidence of heterogeneity across SNP-specific MR estimates",
        TRUE ~ "No strong evidence of heterogeneity"
      )
    )
}

## --------------------------------------------------------------
## 3. MR-Egger intercept / directional horizontal pleiotropy test
##
## Interpretation:
##   - p < 0.05 suggests evidence of directional horizontal pleiotropy.
##   - p >= 0.05 means no strong evidence of directional horizontal pleiotropy.
## --------------------------------------------------------------
egger_res <- tryCatch(
  {
    er <- TwoSampleMR::mr_pleiotropy_test(harm) %>%
      add_study_info()
    
    message("MR-Egger intercept pleiotropy test completed.")
    print(er)
    er
  },
  error = function(e) {
    message("mr_pleiotropy_test failed for ",
            paste(label_exp, collapse = ", "),
            " → ",
            paste(label_out, collapse = ", "),
            ": ",
            e$message)
    tibble()
  }
)

## Optional clean interpretation column for pleiotropy
if (nrow(egger_res) > 0 && "pval" %in% names(egger_res)) {
  egger_res <- egger_res %>%
    mutate(
      pleiotropy_interpretation = case_when(
        is.na(pval) ~ NA_character_,
        pval < 0.05 ~ "Evidence of directional horizontal pleiotropy",
        TRUE ~ "No strong evidence of directional horizontal pleiotropy"
      )
    )
}

## --------------------------------------------------------------
## 4. Save outputs
## --------------------------------------------------------------

# 4a. MR results
out_mr_file <- file.path(
  step5_dir,
  "univariable_MR_TG_LDL_alone_on_CAD_step5.csv"
)

write_csv(mr_res, out_mr_file)
message("Wrote univariable MR summary to: ", out_mr_file)

# 4b. Cochran's Q heterogeneity results
if (nrow(het_res) > 0) {
  het_file <- file.path(
    step5_dir,
    "Heterogeneity_Cochran_Q_TG_LDL_alone_on_CAD_step5.csv"
  )
  
  write_csv(het_res, het_file)
  message("Wrote Cochran's Q heterogeneity results to: ", het_file)
} else {
  message("No heterogeneity results were written because het_res is empty.")
}

# 4c. MR-Egger intercept pleiotropy results
if (nrow(egger_res) > 0) {
  pleio_file <- file.path(
    step5_dir,
    "Egger_intercept_TG_LDL_alone_on_CAD_step5.csv"
  )
  
  write_csv(egger_res, pleio_file)
  message("Wrote MR-Egger intercept pleiotropy results to: ", pleio_file)
} else {
  message("No Egger intercept results were written because egger_res is empty.")
}

# 4d. Copy harmonised file into Step 5 for convenience
harm_copy <- file.path(
  step5_dir,
  "harm_TG_LDL_alone_CAD_step5.csv"
)

fwrite(harm, harm_copy)
message("Copied harmonised CAD data from Step 4 to: ", harm_copy)

message("--------------------------------------------------")
message("Step 5 complete: MR + Cochran's Q heterogeneity + Egger pleiotropy.")
message("Main files saved in: ", step5_dir)






