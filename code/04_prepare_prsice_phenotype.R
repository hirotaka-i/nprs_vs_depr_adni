# Prepare PRSice-2 input files by merging phenotype data with PCs

library(tidyverse)

# Set paths
pheno_y1_file <- "temp/depr1_y1.csv"
pheno_y2_file <- "temp/depr1_y2.csv"
pc_file <- "temp/adni1_pcs.eigenvec"
output_y1 <- "temp/depr1_y1_prsice.txt"
output_y2 <- "temp/depr1_y2_prsice.txt"

cat("=============================================\n")
cat("Step 4: Prepare PRSice-2 Phenotype Files\n")
cat("=============================================\n\n")

# Read phenotype files
cat("Reading phenotype files...\n")
pheno_y1 <- read.csv(pheno_y1_file, stringsAsFactors = FALSE)
pheno_y2 <- read.csv(pheno_y2_file, stringsAsFactors = FALSE)

cat("Year 1 phenotype: ", nrow(pheno_y1), "samples\n")
cat("Year 2 phenotype: ", nrow(pheno_y2), "samples\n\n")

# Read PC file
cat("Reading PC file...\n")
# PLINK2 eigenvec format: FID IID PC1 PC2 ... PC10
pcs <- read.table(pc_file, header = FALSE, stringsAsFactors = FALSE)
colnames(pcs) <- c("FID", "IID", paste0("PC", 1:10))

cat("PC file: ", nrow(pcs), "samples\n\n")

# Function to prepare phenotype file
prepare_prsice_pheno <- function(pheno_df, pcs_df, cohort_name) {
  
  cat("Processing", cohort_name, "...\n")
  
  # Merge with PCs
  merged <- pheno_df %>%
    left_join(pcs_df, by = "IID")
  
  cat("  After merging with PCs:", nrow(merged), "samples\n")
  
  # Check for missing PCs
  missing_pcs <- sum(is.na(merged$PC1))
  if (missing_pcs > 0) {
    cat("  WARNING:", missing_pcs, "samples missing PC data\n")
    merged <- merged %>% filter(!is.na(PC1))
    cat("  After removing samples without PCs:", nrow(merged), "samples\n")
  }
  
  # Convert PTGENDER to binary (Male=1, Female=0)
  merged <- merged %>%
    mutate(PTGENDER = ifelse(PTGENDER == "Male", 1, 0))
  
  cat("  Gender distribution (0=Female, 1=Male):\n")
  print(table(merged$PTGENDER))
  
  # Create dummy variables for DX.bl (CN as reference)
  cat("  Baseline diagnosis distribution:\n")
  print(table(merged$DX.bl))
  
  merged <- merged %>%
    mutate(
      DX_MCI = ifelse(DX.bl == "MCI", 1, 0),
      DX_AD = ifelse(DX.bl == "AD", 1, 0)
    )
  
  # Create dummy variables for SITE (first site alphabetically as reference)
  cat("  Site distribution:\n")
  site_table <- table(merged$SITE)
  print(site_table)
  
  # Get unique sites and create dummies (exclude first as reference)
#   sites <- sort(unique(merged$SITE))
#   cat("  Using", sites[1], "as reference site\n")
  
#   for (i in 2:length(sites)) {
#     site_name <- paste0("SITE_", sites[i])
#     merged[[site_name]] <- ifelse(merged$SITE == sites[i], 1, 0)
#   }
  # Too many sites! Drop the sites variable for this analysis. 
  
  # Create FID column (same as IID for unrelated individuals)
  merged <- merged %>%
    mutate(FID = 0)
  
  # Select and reorder columns: FID, IID first, then phenotype and covariates
  site_cols <- grep("^SITE_", colnames(merged), value = TRUE)
  
  final_df <- merged %>%
    select(FID, IID, DEPR1, AGE, PTGENDER, APOE4, PTEDUCAT, 
           DX_MCI, DX_AD, all_of(site_cols),
           PC1, PC2, PC3, PC4, PC5, PC6, PC7, PC8, PC9, PC10)
  
  cat("  Final dimensions:", nrow(final_df), "samples x", ncol(final_df), "variables\n")
  cat("  Columns:", paste(colnames(final_df), collapse=", "), "\n")
  
  return(final_df)
}

# Process Year 1
cat("\n")
cat("==========================================\n")
prsice_y1 <- prepare_prsice_pheno(pheno_y1, pcs, "Year 1")

# Process Year 2
cat("\n")
cat("==========================================\n")
prsice_y2 <- prepare_prsice_pheno(pheno_y2, pcs, "Year 2")

# Check for missing data
cat("\n")
cat("==========================================\n")
cat("Missing Data Summary\n")
cat("==========================================\n")

cat("\nYear 1:\n")
missing_y1 <- colSums(is.na(prsice_y1))
if (any(missing_y1 > 0)) {
  print(missing_y1[missing_y1 > 0])
} else {
  cat("  No missing data\n")
}

cat("\nYear 2:\n")
missing_y2 <- colSums(is.na(prsice_y2))
if (any(missing_y2 > 0)) {
  print(missing_y2[missing_y2 > 0])
} else {
  cat("  No missing data\n")
}

# Save outputs
cat("\n")
cat("==========================================\n")
cat("Saving output files...\n")
cat("==========================================\n")

write.table(prsice_y1, output_y1, 
            row.names = FALSE, quote = FALSE, sep = "\t")
cat("Year 1 phenotype saved to:", output_y1, "\n")
cat("  Samples:", nrow(prsice_y1), "\n")

write.table(prsice_y2, output_y2, 
            row.names = FALSE, quote = FALSE, sep = "\t")
cat("Year 2 phenotype saved to:", output_y2, "\n")
cat("  Samples:", nrow(prsice_y2), "\n")

cat("\n")
cat("===========================================\n")
cat("SUMMARY\n")
cat("===========================================\n")
cat("Phenotype files ready for PRSice-2\n")
