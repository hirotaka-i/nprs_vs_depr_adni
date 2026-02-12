#!/bin/bash
# cv_final_model.sh
# Train final PRS model on full cohort with user-specified parameters
# Usage: bash cv_final_model.sh <phenotype.pheno> <clump_r2> <p_threshold>

set -e

# Check arguments
if [ $# -ne 3 ]; then
    echo "Usage: bash code/cv_final_model.sh <phenotype.pheno> <clump_r2> <p_threshold>"
    echo "Example: bash code/cv_final_model.sh temp/npi_dpr_ever.pheno 0.1 0.005"
    echo ""
    echo "Run cv_grid_search.sh first to determine optimal parameters"
    exit 1
fi

PHENO_FILE=$1
CLUMP_R2=$2
P_THRESHOLD=$3

# Extract outcome name
OUTCOME=$(basename ${PHENO_FILE} .pheno)
CV_DIR="temp/cv5/${OUTCOME}"
FINAL_DIR="${CV_DIR}/final"

# Create output directory
mkdir -p ${FINAL_DIR}

# Fixed paths
GWAS_FILE="temp/neuroticism_gwas_qc.txt.gz"
GENO_PREFIX="temp/adni1_qc_final"

echo "========================================================================"
echo "Training Final PRS Model"
echo "========================================================================"
echo "Outcome: ${OUTCOME}"
echo "Phenotype file: ${PHENO_FILE}"
echo "Parameters:"
echo "  clump-r2: ${CLUMP_R2}"
echo "  p-threshold: ${P_THRESHOLD}"
echo "Output directory: ${FINAL_DIR}"
echo ""

# Check if files exist
if [ ! -f "${PHENO_FILE}" ]; then
    echo "ERROR: Phenotype file not found: ${PHENO_FILE}"
    exit 1
fi

if [ ! -f "${GWAS_FILE}" ]; then
    echo "ERROR: GWAS file not found: ${GWAS_FILE}"
    exit 1
fi

if [ ! -f "${GENO_PREFIX}.bed" ]; then
    echo "ERROR: Genotype files not found: ${GENO_PREFIX}.*"
    exit 1
fi

# Check if PRSice is available
if ! command -v PRSice &> /dev/null; then
    echo "ERROR: PRSice not found in PATH"
    exit 1
fi

# Validate parameters
if ! [[ ${CLUMP_R2} =~ ^[0-9]*\.?[0-9]+$ ]]; then
    echo "ERROR: Invalid clump-r2 value: ${CLUMP_R2}"
    exit 1
fi

if ! [[ ${P_THRESHOLD} =~ ^[0-9]*\.?[0-9]+([eE][-+]?[0-9]+)?$ ]]; then
    echo "ERROR: Invalid p-threshold value: ${P_THRESHOLD}"
    exit 1
fi

OUT_PREFIX="${FINAL_DIR}/full_cohort"

# Record parameters
PARAM_FILE="${FINAL_DIR}/parameters_used.txt"
cat > ${PARAM_FILE} << EOF
Final PRS Model Parameters
==========================

Date: $(date)
Outcome: ${OUTCOME}
Phenotype file: ${PHENO_FILE}

PRSice Parameters:
  GWAS file: ${GWAS_FILE}
  Genotype file: ${GENO_PREFIX}
  Clumping r²: ${CLUMP_R2}
  Clumping kb: 250
  P-value threshold: ${P_THRESHOLD}
  
Covariates:
  AGE, PTGENDER, APOE4, PTEDUCAT, DX_MCI, DX_AD, PC1, PC2, PC3, PC4

Notes:
  - Trained on full cohort (all individuals)
  - Parameters selected from 5-fold cross-validation
  - See ${CV_DIR}/results/ for CV results
  
WARNING: Performance metrics from this run are optimistic.
Use cross-validation test results for performance claims.
EOF

echo "Parameters recorded in: ${PARAM_FILE}"
echo ""

# Run PRSice on full cohort
echo "========================================================================"
echo "Running PRSice on Full Cohort..."
echo "========================================================================"

PRSice.R \
    --prsice $(which PRSice) \
    --base ${GWAS_FILE} \
    --target ${GENO_PREFIX} \
    --snp SNP \
    --chr CHR \
    --bp BP \
    --A1 A1 \
    --A2 A2 \
    --stat BETA \
    --beta \
    --binary-target T \
    --pvalue P \
    --pheno ${PHENO_FILE} \
    --pheno-col PHENO \
    --cov ${PHENO_FILE} \
    --cov-col AGE,PTGENDER,APOE4,PTEDUCAT,DX_MCI,DX_AD,PC1,PC2,PC3,PC4 \
    --clump-kb 250 \
    --clump-r2 ${CLUMP_R2} \
    --clump-p 1 \
    --bar-levels ${P_THRESHOLD} \
    --fastscore \
    --all-score \
    --print-snp \
    --thread 4 \
    --out ${OUT_PREFIX}

# Check if PRSice completed successfully
if [ $? -ne 0 ]; then
    echo ""
    echo "ERROR: PRSice failed!"
    echo "Check log file: ${OUT_PREFIX}.log"
    exit 1
fi

echo ""
echo "========================================================================"
echo "Final Model Training Complete!"
echo "========================================================================"
echo ""
echo "Output files:"
echo "  ${OUT_PREFIX}.best - PRS scores for all individuals"
echo "  ${OUT_PREFIX}.summary - Performance summary (OPTIMISTIC - for reference only)"
echo "  ${OUT_PREFIX}.prsice - Detailed results"
echo "  ${OUT_PREFIX}.snp - SNPs included in PRS"
echo "  ${OUT_PREFIX}.log - PRSice log file"
echo "  ${PARAM_FILE} - Parameter record"
echo ""

# Extract and display key information
if [ -f "${OUT_PREFIX}.best" ]; then
    N_INDIVIDUALS=$(tail -n +2 ${OUT_PREFIX}.best | wc -l)
    echo "PRS calculated for ${N_INDIVIDUALS} individuals"
fi

if [ -f "${OUT_PREFIX}.snp" ]; then
    # Count SNPs at the chosen threshold
    N_SNPS=$(tail -n +2 ${OUT_PREFIX}.snp | wc -l)
    echo "Number of SNPs in final PRS: ${N_SNPS}"
fi

echo ""
echo "========================================================================"
echo "IMPORTANT NOTES"
echo "========================================================================"
echo ""
echo "1. Performance metrics from this run are OPTIMISTIC (trained on full data)"
echo "   Use CV test results from ${CV_DIR}/results/ for performance claims"
echo ""
echo "2. The .best file contains final PRS scores ready for:"
echo "   - External validation studies"
echo "   - Association analyses"
echo "   - Publication/sharing"
echo ""
echo "3. Record the parameters used:"
echo "   clump-r2 = ${CLUMP_R2}"
echo "   p-threshold = ${P_THRESHOLD}"
echo ""
echo "========================================================================"
