#!/bin/bash
# 05_run_prsice2.sh
# Calculate Neuroticism PRS in ADNI1 and test association with outcomes

# Check arguments
if [ $# -eq 0 ]; then
    echo "Usage: bash 05_run_prsice2.sh <phenotype_file.pheno>"
    echo "Example: bash 05_run_prsice2.sh temp/gds5_y1.pheno"
    exit 1
fi

# Read paths from environment variables if .env exists
if [ -f .env ]; then
    source .env
fi

# Input phenotype file
PHENO_FILE=$1

# Extract base name without extension for output prefix
PHENO_BASE=$(basename ${PHENO_FILE} .pheno)
OUT_PREFIX="${PHENO_FILE%.pheno}.prs"

# Set other paths
GWAS_FILE="temp/neuroticism_gwas_qc.txt.gz"
GENO_PREFIX="temp/adni1_qc_final"

echo "========================================================================"
echo "Step 5: Calculate Neuroticism PRS and Test Association"
echo "========================================================================"
echo ""
echo "GWAS summary statistics: ${GWAS_FILE}"
echo "Target genotype data: ${GENO_PREFIX}"
echo "Phenotype file: ${PHENO_FILE}"
echo "Output prefix: ${OUT_PREFIX}"
echo ""

# Check if files exist
if [ ! -f "${GWAS_FILE}" ]; then
    echo "ERROR: GWAS file not found: ${GWAS_FILE}"
    exit 1
fi

if [ ! -f "${GENO_PREFIX}.bed" ]; then
    echo "ERROR: Genotype files not found: ${GENO_PREFIX}.*"
    exit 1
fi

if [ ! -f "${PHENO_FILE}" ]; then
    echo "ERROR: Phenotype file not found: ${PHENO_FILE}"
    exit 1
fi

# Check if PRSice is available
if ! command -v PRSice &> /dev/null; then
    echo "ERROR: PRSice not found in PATH"
    echo "Please install PRSice or add it to your PATH"
    exit 1
fi

# Run PRSice-2
echo "========================================================================"
echo "Running PRSice-2..."
echo "========================================================================"
echo ""

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
    --clump-r2 0.1 \
    --clump-p 1 \
    --bar-levels 5e-8,5e-7,5e-6,5e-5,5e-4,5e-3,5e-2,5e-1,1 \
    --fastscore \
    --all-score \
    --print-snp \
    --thread 2 \
    --out ${OUT_PREFIX}

# Check if PRSice completed successfully
if [ $? -eq 0 ]; then
    echo ""
    echo "========================================================================"
    echo "PRSice-2 COMPLETE!"
    echo "========================================================================"
    echo ""
    echo "Output files:"
    echo "  ${OUT_PREFIX}.summary - Summary statistics"
    echo "  ${OUT_PREFIX}.prsice - Best PRS results per threshold"
    echo "  ${OUT_PREFIX}.best - Best threshold PRS per individual"
    echo "  ${OUT_PREFIX}.all_score - All threshold scores"
    echo "  ${OUT_PREFIX}.snp - SNPs included at each threshold"
    echo "  ${OUT_PREFIX}.log - PRSice log file"
    echo "  ${OUT_PREFIX}_BARPLOT*.png - Bar plot of results"
    echo "  ${OUT_PREFIX}_HIGH-RES_PLOT*.png - High-resolution plot"
    echo ""
    
    # Display quick results if summary exists
    if [ -f "${OUT_PREFIX}.summary" ]; then
        echo "========================================================================"
        echo "Quick Results Preview"
        echo "========================================================================"
        echo ""
        head -2 ${OUT_PREFIX}.summary | column -t
        echo ""
    fi
else
    echo ""
    echo "========================================================================"
    echo "ERROR: PRSice-2 failed!"
    echo "========================================================================"
    echo "Check ${OUT_PREFIX}.log for details"
    exit 1
fi

echo "========================================================================"
echo "Analysis complete for: ${PHENO_BASE}"
echo "========================================================================"