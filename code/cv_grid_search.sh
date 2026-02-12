#!/bin/bash
# cv_grid_search.sh
# Master script for 5-fold cross-validation grid search for PRS optimization
# Usage: bash code/cv_grid_search.sh <phenotype.pheno>

set -e  # Exit on error

# Check arguments
if [ $# -ne 1 ]; then
    echo "Usage: bash code/cv_grid_search.sh <phenotype.pheno>"
    echo "Example: bash code/cv_grid_search.sh temp/npi_dpr_ever.pheno"
    exit 1
fi

PHENO_FILE=$1

# Check if phenotype file exists
if [ ! -f "${PHENO_FILE}" ]; then
    echo "ERROR: Phenotype file not found: ${PHENO_FILE}"
    exit 1
fi

# Extract outcome name from filename
OUTCOME=$(basename ${PHENO_FILE} .pheno)
CV_DIR="temp/cv5/${OUTCOME}"

echo "========================================================================"
echo "PRSice 5-Fold Cross-Validation Grid Search"
echo "========================================================================"
echo "Phenotype file: ${PHENO_FILE}"
echo "Outcome: ${OUTCOME}"
echo "CV directory: ${CV_DIR}"
echo ""

# Create output directory
mkdir -p ${CV_DIR}/{folds,training,testing,results,logs}

# Fixed paths
GWAS_FILE="temp/neuroticism_gwas_qc.txt.gz"
GENO_PREFIX="temp/adni1_qc_final"

# Check required files
if [ ! -f "${GWAS_FILE}" ]; then
    echo "ERROR: GWAS file not found: ${GWAS_FILE}"
    exit 1
fi

if [ ! -f "${GENO_PREFIX}.bed" ]; then
    echo "ERROR: Genotype files not found: ${GENO_PREFIX}.*"
    exit 1
fi

# ============================================================================
# Step 1: Create stratified folds
# ============================================================================
echo "========================================================================"
echo "Step 1: Creating 5 stratified folds"
echo "========================================================================"

Rscript code/cv_create_folds.R ${PHENO_FILE} ${CV_DIR}/folds

if [ $? -ne 0 ]; then
    echo "ERROR: Fold creation failed"
    exit 1
fi

echo "Folds created successfully"
echo ""

# ============================================================================
# Step 2-3: Training phase (grid search)
# ============================================================================
echo "========================================================================"
echo "Step 2-3: Training phase (hyperparameter tuning)"
echo "========================================================================"

bash code/cv_train_grid.sh ${PHENO_FILE} ${CV_DIR} ${GWAS_FILE} ${GENO_PREFIX}

if [ $? -ne 0 ]; then
    echo "ERROR: Training phase failed"
    exit 1
fi

echo "Training phase completed successfully"
echo ""

# ============================================================================
# Step 4: Test phase (held-out evaluation)
# ============================================================================
echo "========================================================================"
echo "Step 4: Test phase (held-out evaluation)"
echo "========================================================================"

bash code/cv_test_best.sh ${PHENO_FILE} ${CV_DIR} ${GWAS_FILE} ${GENO_PREFIX}

if [ $? -ne 0 ]; then
    echo "ERROR: Test phase failed"
    exit 1
fi

echo "Test phase completed successfully"
echo ""

# ============================================================================
# Step 5: Summarize results
# ============================================================================
echo "========================================================================"
echo "Step 5: Summarizing CV results"
echo "========================================================================"

Rscript code/cv_summarize.R ${CV_DIR}

if [ $? -ne 0 ]; then
    echo "ERROR: Summarization failed"
    exit 1
fi

echo ""
echo "========================================================================"
echo "CROSS-VALIDATION COMPLETE!"
echo "========================================================================"
echo ""
echo "Results saved to: ${CV_DIR}/results/"
echo ""
echo "Review the following files:"
echo "  - ${CV_DIR}/results/cv_summary.txt"
echo "  - ${CV_DIR}/results/recommended_params.txt"
echo "  - ${CV_DIR}/results/cv_results_aggregated.csv"
echo ""
echo "Next step: Train final model with chosen parameters"
echo "Command: bash code/cv_final_model.sh ${PHENO_FILE} <clump_r2> <p_threshold>"
echo ""
echo "========================================================================"
