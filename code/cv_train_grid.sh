#!/bin/bash
# cv_train_grid.sh
# Training phase: Grid search over clump-r2 and p-thresholds
# Usage: bash cv_train_grid.sh <pheno_file> <cv_dir> <gwas_file> <geno_prefix>

set -e

PHENO_FILE=$1
CV_DIR=$2
GWAS_FILE=$3
GENO_PREFIX=$4

echo "========================================================================"
echo "Training Phase: Grid Search"
echo "========================================================================"
echo "Phenotype: ${PHENO_FILE}"
echo "CV directory: ${CV_DIR}"
echo ""

# Parameter grid
CLUMP_R2_VALUES=(0.01 0.05 0.1 0.2 0.5)
P_THRESHOLDS="5e-8,1e-6,1e-5,1e-4,1e-3,5e-3,1e-2,5e-2,1e-1,0.5,1"

# Check if PRSice is available
if ! command -v PRSice &> /dev/null; then
    echo "ERROR: PRSice not found in PATH"
    exit 1
fi

# Training loop
total_runs=$((5 * ${#CLUMP_R2_VALUES[@]}))
current_run=0

for fold in {1..5}; do
  KEEP_FILE="${CV_DIR}/folds/train_fold_${fold}.keep"
  
  if [ ! -f "${KEEP_FILE}" ]; then
    echo "ERROR: Keep file not found: ${KEEP_FILE}"
    exit 1
  fi
  
  for r2 in "${CLUMP_R2_VALUES[@]}"; do
    current_run=$((current_run + 1))
    
    # Format r2 for filename (replace . with p)
    r2_label=$(echo ${r2} | sed 's/\./_/g')
    
    OUT_PREFIX="${CV_DIR}/training/fold_${fold}_r2_${r2_label}"
    
    echo ""
    echo "[$current_run/$total_runs] Running fold ${fold}, r²=${r2}"
    echo "Output: ${OUT_PREFIX}"
    
    # Run PRSice
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
      --clump-r2 ${r2} \
      --clump-p 1 \
      --bar-levels ${P_THRESHOLDS} \
      --fastscore \
      --thread 2 \
      --out ${OUT_PREFIX} \
      > ${CV_DIR}/logs/train_fold_${fold}_r2_${r2_label}.log 2>&1
    
    if [ $? -ne 0 ]; then
      echo "ERROR: PRSice failed for fold ${fold}, r²=${r2}"
      echo "Check log: ${CV_DIR}/logs/train_fold_${fold}_r2_${r2_label}.log"
      exit 1
    fi
    
    echo "Completed fold ${fold}, r²=${r2}"
  done
done

echo ""
echo "========================================================================"
echo "Training Phase Complete!"
echo "========================================================================"
echo "Total runs: ${total_runs}"
echo "Results saved to: ${CV_DIR}/training/"
echo ""

# Now select best parameters for each fold
echo "========================================================================"
echo "Selecting Best Parameters Per Fold"
echo "========================================================================"

# Create a summary file for best parameters per fold
BEST_PARAMS_FILE="${CV_DIR}/training/best_params_per_fold.csv"
echo "fold,best_r2,best_p,best_R2,best_P" > ${BEST_PARAMS_FILE}

for fold in {1..5}; do
  echo "Processing fold ${fold}..."
  
  # Find all .prsice files for this fold
  best_R2=0
  best_r2=""
  best_p=""
  best_P=""
  
  for r2 in "${CLUMP_R2_VALUES[@]}"; do
    r2_label=$(echo ${r2} | sed 's/\./_/g')
    prsice_file="${CV_DIR}/training/fold_${fold}_r2_${r2_label}.prsice"
    
    if [ ! -f "${prsice_file}" ]; then
      echo "WARNING: File not found: ${prsice_file}"
      continue
    fi
    
    # Parse .prsice file (skip header, find max R2)
    while IFS=$'\t' read -r pheno set threshold R2 P coef se num_snp; do
      # Skip header and non-numeric lines
      if [[ ${R2} =~ ^[0-9] ]]; then
        # Compare R2 values (bash doesn't do floating point, use awk)
        is_better=$(awk -v r2="${R2}" -v best="${best_R2}" 'BEGIN {print (r2 > best)}')
        if [ "${is_better}" = "1" ]; then
          best_R2=${R2}
          best_r2=${r2}
          best_p=${threshold}
          best_P=${P}
        fi
      fi
    done < <(tail -n +2 ${prsice_file})
  done
  
  if [ -z "${best_r2}" ]; then
    echo "ERROR: Could not find best parameters for fold ${fold}"
    exit 1
  fi
  
  echo "${fold},${best_r2},${best_p},${best_R2},${best_P}" >> ${BEST_PARAMS_FILE}
  echo "Fold ${fold}: best r²=${best_r2}, p=${best_p}, R²=${best_R2}"
done

echo ""
echo "Best parameters saved to: ${BEST_PARAMS_FILE}"
cat ${BEST_PARAMS_FILE}
echo ""
