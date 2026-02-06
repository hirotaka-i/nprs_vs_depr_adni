#!/bin/bash
# 05_run_prsice2.sh
# Calculate Neuroticism PRS in ADNI1 and test association with depression

# Read paths from environment variables
source .env

# Set paths
GWAS_FILE="temp/neuroticism_gwas_qc.txt.gz"
GENO_PREFIX="temp/adni1_qc_final"
PHENO_FILE="temp/depr1_y2_prsice.txt"
OUT_DIR="temp/prsice"
OUT_PREFIX="${OUT_DIR}/nprs_depr1_y2"

# Create output directory
mkdir -p ${OUT_DIR}


# Run PRSice-2
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
    --binary-target T\
    --pvalue P \
    --pheno ${PHENO_FILE} \
    --pheno-col DEPR1 \
    --cov ${PHENO_FILE} \
    --cov-col AGE,PTGENDER,APOE4,PTEDUCAT,DX_MCI,DX_AD,PC1,PC2 \
    --clump-kb 250 \
    --clump-r2 0.1 \
    --clump-p 1 \
    --bar-levels 5e-8,5e-7,5e-6,5e-5,5e-4,5e-3,5e-2,5e-1,1 \
    --fastscore \
    --all-score \
    --print-snp \
    --thread 4 \
    --out ${OUT_PREFIX}

cp ${OUT_PREFIX}*.png report/
cp ${OUT_PREFIX}.log report/
cp ${OUT_PREFIX}.summary report/
cp ${OUT_PREFIX}.prsice report/
