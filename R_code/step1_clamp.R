##################################################################
###   Step 1: Clamp the SNP region for TG-LDL alone              #
###
###   Hugh
###  Modified 12-4-2025
#####################################################################

suppressPackageStartupMessages({
  library(data.table)
  library(dplyr)
  library(tibble)
  library(stringr)
})

# ---- paths (edit if yours differ) ----
base <- "/Users/ylu6/Desktop/ADA_2026"

ref_prefix <- file.path(
  base,
  "resource_and_gwas_data", "1000G_inputs",
  "1000G.EUR.QC"
)

stats_dir  <- file.path(
  base,
  "resource_and_gwas_data", "UKB_inputs",
  "TG_LDL"
)

out_dir   <- file.path(base, "output", "step1")
plink_bin <- file.path(base, "plink_mac", "plink")

# ---- clumping knobs (for TG-LDL alone) ----
p_threshold <- 1e-6      # loosen (e.g. 5e-6 / 1e-5 / 5e-5) if you want more IVs
clump_r2    <- 0.01
clump_kb    <- 100       # can try 250/1000 as sensitivity
clump_p2    <- 1e-4
# --------------------------------------

# Basic checks
stopifnot(
  file.exists(paste0(ref_prefix, ".bed")),
  file.exists(paste0(ref_prefix, ".bim")),
  file.exists(paste0(ref_prefix, ".fam"))
)
stopifnot(file.exists(plink_bin))

dir.create(out_dir, showWarnings = FALSE, recursive = TRUE)

# --------------------------------------------------
# Helper: detect p-value column from BOLT-LMM output
# --------------------------------------------------
detect_pcol <- function(nm) {
  cand <- c("P_BOLT_LMM_INF", "P_BOLT_LMM", "P_BOLT", "P", "p_value", "PVAL")
  first <- intersect(cand, nm)
  if (!length(first)) stop("No recognized p-value column in file header")
  first[1]
}

# --------------------------------------------------
# Helper: read a GWAS file into a standardized format
# --------------------------------------------------
read_gwas_standard <- function(fpath) {
  nm <- basename(fpath)
  message("Reading GWAS file: ", nm)

  gwas <- fread(fpath)
  pcol <- detect_pcol(names(gwas))

  gwas %>%
    transmute(
      SNP           = SNP,
      CHR           = CHR,
      BP            = BP,
      effect_allele = ALLELE1,
      other_allele  = ALLELE0,
      eaf           = A1FREQ,
      beta          = BETA,
      se            = SE,
      pval          = .data[[pcol]]
    ) %>%
    distinct(SNP, .keep_all = TRUE)
}

# --------------------------------------------------
# Helper: PLINK clumping for TG-LDL (single trait)
# --------------------------------------------------
clump_single_trait <- function(raw, stem) {
  message("Pre-clump hits (p < ", p_threshold, "): ", nrow(raw))
  if (nrow(raw) == 0) {
    warning("No variants pass p < ", p_threshold, " for ", stem)
    return(invisible(NULL))
  }

  tmpdir <- tempfile("clump_")
  dir.create(tmpdir)

  in_txt <- file.path(tmpdir, "clump_input.txt")
  fwrite(raw %>% select(SNP, P = pval), in_txt, sep = "\t")

  out_prefix <- file.path(tmpdir, "clump_out")
  cmd <- sprintf(
    '"%s" --bfile "%s" --clump "%s" --clump-p1 %.3g --clump-p2 %.3g --clump-r2 %.3f --clump-kb %d --out "%s"',
    plink_bin, ref_prefix, in_txt,
    p_threshold, clump_p2, clump_r2, clump_kb, out_prefix
  )
  message("Running PLINK clumping for ", stem, ":\n", cmd)
  status <- system(cmd)
  if (status != 0) stop("PLINK clumping failed for ", stem)

  clumped <- fread(paste0(out_prefix, ".clumped"), fill = TRUE)
  lead_snps <- unique(na.omit(clumped$SNP))

  keep <- raw %>%
    semi_join(tibble(SNP = lead_snps), by = "SNP") %>%
    # Drop ambiguous A/T or C/G with missing EAF
    filter(
      !(effect_allele %in% c("A", "T") & other_allele %in% c("A", "T") & is.na(eaf)) &
      !(effect_allele %in% c("C", "G") & other_allele %in% c("C", "G") & is.na(eaf))
    )

  out_csv <- file.path(out_dir, paste0(stem, "_instruments.csv"))
  fwrite(keep, out_csv)
  message("Wrote: ", out_csv, " | kept instruments: ", nrow(keep))
}

# --------------------------------------------------
# Main: load TG-LDL GWAS file and run clumping
# --------------------------------------------------

tg_file <- file.path(stats_dir, "bolt.stats.TG_LDL_alone.gz")
if (!file.exists(tg_file)) stop("TG-LDL GWAS file not found: ", tg_file)

message("Processing TG-LDL_alone GWAS file: ", basename(tg_file))
tg_raw  <- read_gwas_standard(tg_file)
tg_hits <- tg_raw %>% filter(pval < p_threshold)

if (nrow(tg_hits) > 0) {
  clump_single_trait(
    raw  = tg_hits,
    stem = "bolt_stats_TG_LDL_alone"
  )
} else {
  warning("No TG-LDL_alone variants pass p < ", p_threshold)
}



