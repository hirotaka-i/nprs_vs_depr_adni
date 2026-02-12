#!/usr/bin/env Rscript
# cv_create_folds.R
# Create 5 stratified folds for cross-validation
# Usage: Rscript cv_create_folds.R <phenotype.pheno> <output_dir>

library(tidyverse)
library(caret)  # For createFolds

# Get command line arguments
args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 2) {
  cat("Usage: Rscript cv_create_folds.R <phenotype.pheno> <output_dir>\n")
  quit(status = 1)
}

pheno_file <- args[1]
output_dir <- args[2]

cat("========================================\n")
cat("Creating 5 Stratified Folds\n")
cat("========================================\n\n")

# Create output directory
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Read phenotype file
cat("Reading phenotype file:", pheno_file, "\n")
pheno <- read.table(pheno_file, header = TRUE, stringsAsFactors = FALSE)

cat("Samples:", nrow(pheno), "\n")
cat("Columns:", paste(colnames(pheno), collapse = ", "), "\n\n")

# Check required columns
required_cols <- c("IID", "PHENO", "PTGENDER", "AGE")
missing_cols <- setdiff(required_cols, colnames(pheno))
if (length(missing_cols) > 0) {
  cat("ERROR: Missing required columns:", paste(missing_cols, collapse = ", "), "\n")
  quit(status = 1)
}

# Create stratification variable
# Combine phenotype, sex, and diagnosis
pheno <- pheno %>%
  mutate(
    age_quartile = cut(AGE, 
                       breaks = quantile(AGE, probs = seq(0, 1, 0.25), na.rm = TRUE),
                       labels = c("Q1", "Q2", "Q3", "Q4"),
                       include.lowest = TRUE),
    strata = paste(PHENO, PTGENDER, DX_MCI, DX_AD, age_quartile, sep = "_")
  )

cat("Stratification summary:\n")
print(table(pheno$strata))
cat("\n")

# Create 5 folds stratified by combined variable
set.seed(42)  # For reproducibility

# Use createFolds from caret package
fold_indices <- createFolds(pheno$strata, k = 5, list = TRUE, returnTrain = FALSE)

# Assign fold numbers to each individual
pheno$fold <- NA
for (i in 1:5) {
  pheno$fold[fold_indices[[i]]] <- i
}

# Check fold distribution
cat("Fold distribution:\n")
print(table(pheno$fold))
cat("\n")

cat("Phenotype distribution by fold:\n")
print(table(pheno$fold, pheno$PHENO))
cat("\n")

cat("Sex distribution by fold:\n")
print(table(pheno$fold, pheno$PTGENDER))
cat("\n")

cat("Diagnosis distribution by fold (MCI):\n")
print(table(pheno$fold, pheno$DX_MCI))
cat("\n")

cat("Diagnosis distribution by fold (AD):\n")
print(table(pheno$fold, pheno$DX_AD))
cat("\n")

cat("Age quartile distribution by fold:\n")
print(table(pheno$fold, pheno$age_quartile))
cat("\n")

cat("Age quartile distribution by fold:\n")
# Save fold assignments
fold_file <- file.path(output_dir, "fold_assignments.csv")
pheno %>%
  select(FID, IID, fold, PHENO, PTGENDER, DX_MCI, DX_AD, AGE) %>%
  write.csv(fold_file, row.names = FALSE, quote = FALSE)

cat("Fold assignments saved to:", fold_file, "\n\n")


# Create keep files for each fold (FID and IID format)
for (k in 1:5) {
  # Training set: all folds except k
  train_data <- pheno %>% 
    filter(fold != k) %>% 
    select(FID, IID)
  train_file <- file.path(output_dir, paste0("train_fold_", k, ".keep"))
  write.table(train_data, train_file, 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  # Test set: only fold k
  test_data <- pheno %>% 
    filter(fold == k) %>% 
    select(FID, IID)
  test_file <- file.path(output_dir, paste0("test_fold_", k, ".keep"))
  write.table(test_data, test_file, 
              row.names = FALSE, col.names = FALSE, quote = FALSE, sep = "\t")
  
  cat(sprintf("Fold %d: %d train, %d test\n", k, nrow(train_data), nrow(test_data)))
}


cat("\n========================================\n")
cat("Fold Creation Complete!\n")
cat("========================================\n")
cat("Output directory:", output_dir, "\n")
cat("Files created:\n")
cat("  - fold_assignments.csv\n")
cat("  - train_fold_1.keep ... train_fold_5.keep\n")
cat("  - test_fold_1.keep ... test_fold_5.keep\n")
