#!/bin/bash
# cv_test_best.sh
# Test phase: Evaluate best parameters on held-out test sets
# Usage: bash cv_test_best.sh <pheno_file> <cv_dir> <gwas_file> <geno_prefix>

set -e

PHENO_FILE=$1
CV_DIR=$2
GWAS_FILE=$3
GENO_PREFIX=$4

echo "========================================================================"
echo "Test Phase: Held-Out Evaluation"
echo "========================================================================"

# Read best parameters file
BEST_PARAMS_FILE="${CV_DIR}/training/best_params_per_fold.csv"

if [ ! -f "${BEST_PARAMS_FILE}" ]; then
  echo "ERROR: Best parameters file not found: ${BEST_PARAMS_FILE}"
  exit 1
fi

# Check if PRSice is available
if ! command -v PRSice &> /dev/null; then
    echo "ERROR: PRSice not found in PATH"
    exit 1
fi

# Create output file for test results
TEST_RESULTS_FILE="${CV_DIR}/testing/test_results.csv"
echo "fold,test_r2,test_p,test_R2,test_P,test_AUC,num_snp" > ${TEST_RESULTS_FILE}

# Process each fold
fold=0
while IFS=, read -r fold_id best_r2 best_p best_R2_train best_P_train; do
  # Skip header
  if [ "${fold_id}" = "fold" ]; then
    continue
  fi
  
  fold=$((fold + 1))
  
  echo ""
  echo "========================================================================"
  echo "Testing Fold ${fold_id}"
  echo "========================================================================"
  echo "Best parameters from training: r²=${best_r2}, p=${best_p}"
  
  KEEP_FILE="${CV_DIR}/folds/test_fold_${fold_id}.keep"
  
  if [ ! -f "${KEEP_FILE}" ]; then
    echo "ERROR: Keep file not found: ${KEEP_FILE}"
    exit 1
  fi
  
  OUT_PREFIX="${CV_DIR}/testing/fold_${fold_id}_test"
  
  # Run PRSice with best parameters on test set
  PRSice.R \
    --prsice $(which PRSice) \
    --base ${GWAS_FILE} \
    --target ${GENO_PREFIX} \
    --keep ${KEEP_FILE} \
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
    --clump-r2 ${best_r2} \
    --clump-p 1 \
    --bar-levels ${best_p} \
    --fastscore \
    --thread 2 \
    --out ${OUT_PREFIX} \
    > ${CV_DIR}/logs/test_fold_${fold_id}.log 2>&1
  
  if [ $? -ne 0 ]; then
    echo "ERROR: PRSice failed for test fold ${fold_id}"
    echo "Check log: ${CV_DIR}/logs/test_fold_${fold_id}.log"
    exit 1
  fi
  
  # Extract test results from .prsice file
  prsice_file="${OUT_PREFIX}.prsice"
  
  if [ ! -f "${prsice_file}" ]; then
    echo "ERROR: Output file not found: ${prsice_file}"
    exit 1
  fi
  
  # Read the second line (first data line after header)
  test_results=$(tail -n 1 ${prsice_file})
  
  # Parse results (format: Pheno Set Threshold R2 P Coefficient Standard.Error Num_SNP)
  test_R2=$(echo ${test_results} | awk '{print $4}')
  test_P=$(echo ${test_results} | awk '{print $5}')
  num_snp=$(echo ${test_results} | awk '{print $8}')
  
  # Try to extract AUC if available (PRSice doesn't always report this in .prsice)
  # We'll set it to NA for now and let the summarize script handle it
  test_AUC="NA"
  
  echo "${fold_id},${best_r2},${best_p},${test_R2},${test_P},${test_AUC},${num_snp}" >> ${TEST_RESULTS_FILE}
  
  echo "Test results: R²=${test_R2}, P=${test_P}, SNPs=${num_snp}"
  
done < ${BEST_PARAMS_FILE}

echo ""
echo "========================================================================"
echo "Test Phase Complete!"
echo "========================================================================"
echo "Results saved to: ${TEST_RESULTS_FILE}"
echo ""
cat ${TEST_RESULTS_FILE}
echo ""
