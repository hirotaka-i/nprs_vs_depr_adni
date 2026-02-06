#!/bin/bash
# QC for ADNI1 imputed genotype data and calculate PCs for population stratification

# Read paths from environment variables
source .env

# Set paths
GENO_PREFIX="${INPUT_PLINK_SUFFIX}"
WORK_DIR="temp/step3_working"
OUT_DIR="temp"
OUT_PREFIX="${OUT_DIR}/adni1_qc_final"
PC_FILE="${OUT_DIR}/adni1_pcs.eigenvec"

# Create directories
mkdir -p ${WORK_DIR}
mkdir -p ${OUT_DIR}

echo "========================================================================"
echo "Step 3: Genotype QC and PCA for ADNI1"
echo "========================================================================"
echo ""
echo "Input: ${GENO_PREFIX}"
echo "Working directory: ${WORK_DIR}"
echo "Output prefix: ${OUT_PREFIX}"
echo "PC output: ${PC_FILE}"
echo ""

# Check if input files exist
if [ ! -f "${GENO_PREFIX}.bed" ] || [ ! -f "${GENO_PREFIX}.bim" ] || [ ! -f "${GENO_PREFIX}.fam" ]; then
    echo "ERROR: Input genotype files not found!"
    exit 1
fi

# Count initial samples and SNPs
echo "Initial dataset:"
plink2 --bfile ${GENO_PREFIX} --freq --out ${WORK_DIR}/initial_stats
N_SAMPLES=$(wc -l < ${GENO_PREFIX}.fam)
N_SNPS=$(wc -l < ${GENO_PREFIX}.bim)
echo "  Samples: ${N_SAMPLES}"
echo "  SNPs: ${N_SNPS}"
echo ""

# Step 1: Keep only autosomes (chr 1-22) and apply HWE filter
echo "========================================================================"
echo "Step 1: Filter to autosomes (chr 1-22) and apply HWE filter (p > 1e-6)"
echo "========================================================================"
plink2 \
    --bfile ${GENO_PREFIX} \
    --chr 1-22 \
    --hwe 1e-6 \
    --make-bed \
    --out ${WORK_DIR}/01_autosome_hwe

N_SNPS_AUTO=$(wc -l < ${WORK_DIR}/01_autosome_hwe.bim)
echo "SNPs after filtering: ${N_SNPS_AUTO}"
echo ""

# Step 2: LD pruning for PCA
echo "========================================================================"
echo "Step 2: LD pruning (--indep-pairwise 50 5 0.2)"
echo "========================================================================"
plink2 \
    --bfile ${WORK_DIR}/01_autosome_hwe \
    --indep-pairwise 50 5 0.2 \
    --out ${WORK_DIR}/02_ld_prune

N_PRUNED=$(wc -l < ${WORK_DIR}/02_ld_prune.prune.in)
echo "LD-pruned SNPs: ${N_PRUNED}"
echo ""

# Step 3: Calculate PCA using LD-pruned SNPs
echo "========================================================================"
echo "Step 3: Calculate principal components (PC1-10)"
echo "========================================================================"
plink2 \
    --bfile ${WORK_DIR}/01_autosome_hwe \
    --extract ${WORK_DIR}/02_ld_prune.prune.in \
    --pca 10 \
    --out ${WORK_DIR}/03_pca

echo "PCA complete. Output files:"
echo "  Eigenvalues: ${WORK_DIR}/03_pca.eigenval"
echo "  Eigenvectors: ${WORK_DIR}/03_pca.eigenvec"
echo ""

# Step 4: Copy final outputs
echo "========================================================================"
echo "Step 4: Create final QC'd genotype files"
echo "========================================================================"
# Use the autosome + HWE filtered data as final output
cp ${WORK_DIR}/01_autosome_hwe.bed ${OUT_PREFIX}.bed
cp ${WORK_DIR}/01_autosome_hwe.bim ${OUT_PREFIX}.bim
cp ${WORK_DIR}/01_autosome_hwe.fam ${OUT_PREFIX}.fam
cp ${WORK_DIR}/03_pca.eigenvec ${PC_FILE}

echo "Final QC'd genotype files:"
echo "  ${OUT_PREFIX}.bed"
echo "  ${OUT_PREFIX}.bim"
echo "  ${OUT_PREFIX}.fam"
echo ""
echo "PC file:"
echo "  ${PC_FILE}"
echo ""

# Summary statistics
echo "========================================================================"
echo "SUMMARY"
echo "========================================================================"
echo "Final dataset:"
N_FINAL_SAMPLES=$(wc -l < ${OUT_PREFIX}.fam)
N_FINAL_SNPS=$(wc -l < ${OUT_PREFIX}.bim)
echo "  Samples: ${N_FINAL_SAMPLES}"
echo "  SNPs: ${N_FINAL_SNPS}"
echo "  PCs calculated: 10"
echo ""
echo "SNPs removed:"
echo "  Non-autosomal + HWE: $((N_SNPS - N_FINAL_SNPS))"
echo ""
echo "PC variance explained:"
head -10 ${WORK_DIR}/03_pca.eigenval | nl
echo ""
echo "========================================================================"
echo "QC and PCA COMPLETE!"
echo "========================================================================"
echo ""
echo "Next step: Merge PCs with phenotype file (step 4)"