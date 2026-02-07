# 04_prepare_prsice_phenotype.R
# Prepare PRSice-2 input files by merging phenotype data with PCs

library(tidyverse)

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  cat("Usage: Rscript 04_prepare_prsice_phenotype.R <outcome_variable> <phenotype_file.csv>\n")
  cat("Example: Rscript 04_prepare_prsice_phenotype.R \"GDS5\" temp/gds5_y1.csv\n")
  quit(status = 1)
}

# Set paths
outcome_var <- args[1]
pheno_file <- args[2]
pc_file <- "temp/adni1_pcs.eigenvec"

# Generate output filename: replace .csv with .pheno
output_file <- sub("\\.csv$", ".pheno", pheno_file)

cat("===========================================\n")
cat("Step 4: Prepare PRSice-2 Phenotype File\n")
cat("===========================================\n\n")
cat("Input file:", pheno_file, "\n")
cat("Output file:", output_file, "\n")
cat("Outcome variable:", outcome_var, "\n\n")

# Read phenotype file
cat("Reading phenotype file...\n")
pheno <- read.csv(pheno_file, stringsAsFactors = FALSE)
cat("Samples:", nrow(pheno), "\n")
cat("Columns:", paste(colnames(pheno), collapse = ", "), "\n\n")

# Check if outcome variable exists
if (!outcome_var %in% colnames(pheno)) {
  cat("ERROR: Outcome variable", outcome_var, "not found in phenotype file\n")
  quit(status = 1)
}

# Read PC file
cat("Reading PC file...\n")
# PLINK2 eigenvec format: FID IID PC1 PC2 ... PC10
pcs <- read.table(pc_file, header = FALSE, stringsAsFactors = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:10))
cat("PC file:", nrow(pcs), "samples\n\n")

# Merge with PCs
cat("Merging phenotype with PCs...\n")
merged <- pheno %>%
  left_join(pcs %>% select(-FID), by = "IID")

cat("After merging:", nrow(merged), "samples\n")

# Check for missing PCs
missing_pcs <- sum(is.na(merged$PC1))
if (missing_pcs > 0) {
  cat("WARNING:", missing_pcs, "samples missing PC data\n")
  merged <- merged %>% filter(!is.na(PC1))
  cat("After removing samples without PCs:", nrow(merged), "samples\n")
}

# Convert PTGENDER to binary (Male=1, Female=0)
cat("\nConverting categorical variables...\n")
merged <- merged %>%
  mutate(PTGENDER = ifelse(PTGENDER == "Male", 1, 0))

cat("Gender distribution (0=Female, 1=Male):\n")
print(table(merged$PTGENDER))

# Create dummy variables for DX.bl (CN as reference)
cat("\nBaseline diagnosis distribution:\n")
print(table(merged$DX.bl))

merged <- merged %>%
  mutate(
    DX_MCI = ifelse(DX.bl == "MCI", 1, 0),
    DX_AD = ifelse(DX.bl == "AD", 1, 0)
  )

# Note: SITE variable dropped due to too many categories
cat("\nNote: SITE variable not included (too many categories)\n")

# Create FID column (set to 0 for all, as IID is unique)
merged <- merged %>%
  mutate(FID = 0)

# Rename outcome variable to standardized name for PRSice-2
merged <- merged %>%
  rename(PHENO = !!sym(outcome_var))

# Select and reorder columns: FID, IID first, then phenotype and covariates
# Include timepoint_censored if it exists (for EVER outcomes)
base_cols <- c("FID", "IID", "PHENO", "AGE", "PTGENDER", "APOE4", "PTEDUCAT", 
               "DX_MCI", "DX_AD")
pc_cols <- paste0("PC", 1:10)

if ("timepoint_censored" %in% colnames(merged)) {
  final_cols <- c(base_cols, "timepoint_censored", pc_cols)
} else {
  final_cols <- c(base_cols, pc_cols)
}

final_df <- merged %>%
  select(all_of(final_cols))

cat("\nFinal dimensions:", nrow(final_df), "samples x", ncol(final_df), "variables\n")
cat("Columns:", paste(colnames(final_df), collapse = ", "), "\n")

# Check for missing data
cat("\nMissing data summary:\n")
missing <- colSums(is.na(final_df))
if (any(missing > 0)) {
  print(missing[missing > 0])
} else {
  cat("No missing data\n")
}

# Outcome distribution
cat("\nOutcome (PHENO) distribution:\n")
print(table(final_df$PHENO, useNA = "ifany"))

# Save output
cat("\nSaving output file...\n")
write.table(final_df, output_file, 
            row.names = FALSE, quote = FALSE, sep = "\t")

cat("\n===========================================\n")
cat("COMPLETE!\n")
cat("===========================================\n")
cat("Output file:", output_file, "\n")
cat("Samples:", nrow(final_df), "\n")
cat("\nFor PRSice-2, use:\n")
cat("  --pheno", output_file, "\n")
cat("  --pheno-col PHENO\n")
cat("  --cov-col AGE,PTGENDER,APOE4,PTEDUCAT,DX_MCI,DX_AD,PC1-PC4\n")
if ("timepoint_censored" %in% colnames(final_df)) {
  cat("  Note: timepoint_censored available for survival analysis\n")
}
cat("===========================================\n")