# Neuroticism PRS and Depression in ADNI

## Overview
Calculate neuroticism polygenic risk scores (N-PRS) in ADNI1 and test association with depression phenotypes at 1-year and 2-year follow-up.

## Input Data
Paths are configured in `.env`:
- **Neuroticism GWAS** (hg38): `$SUMSTATS`
- **ADNI1 Genotypes** (hg38): `$INPUT_PLINK_SUFFIX`
- **ADNI Phenotypes**: `$ADNIMERGE_DATA`

## Analysis Pipeline

### Step 1: Create Phenotype Files
**Script:** `01_create_adni1_phenotype_file.R`

- Extract depression phenotypes from ADNI (adnimerge + gdscale packages)
- Depression defined by Geriatric Depression Scale (GDS):
  - DEPR1: GDS ≥5 (mild threshold)
  - DEPR2: GDS ≥10 (moderate threshold)
- Match with genetic data (fam file) by IID
- Extract cross-sectional data at 12-month and 24-month timepoints
- **Outputs:** `temp/depr1_y1.csv`, `temp/depr1_y2.csv`
  - Columns: IID, DEPR1, AGE, PTGENDER, APOE4, PTEDUCAT, SITE, DX.bl

### Step 2: QC GWAS Summary Statistics
**Script:** `02_qc_gwas_sumstats.py`

- Standard QC filters:
  - Autosomes + X/Y only (remove alt contigs)
  - INFO > 0.8, MAF > 0.01
  - Remove ambiguous SNPs (A/T, G/C)
  - Remove invalid/extreme values
  - Remove duplicates
- Calculate BETA and SE from Z-scores
- **Output:** `temp/neuroticism_gwas_qc.txt.gz` (PRSice-2 format)

### Step 3: QC Genotypes and Calculate PCs
**Script:** `03_qc_genotype_pca.sh`

- Filter to autosomes (chr 1-22)
- Apply HWE filter (p > 1e-6)
- LD pruning for PCA (--indep-pairwise 50 5 0.2)
- Calculate PC1-10 for population stratification
- **Outputs:** 
  - `temp/adni1_qc_final.*` (bed/bim/fam)
  - `temp/adni1_pcs.eigenvec`

### Step 4: Prepare PRSice-2 Input
**Script:** `04_prepare_prsice_phenotype.R`

- Merge phenotype data with PCs (PC1-10)
- Convert categorical variables:
  - PTGENDER: Male=1, Female=0
  - DX.bl: Dummy code (CN=reference, MCI, AD)
  - SITE: Dummy code (first site=reference)
- Format: FID, IID as first columns (PRSice-2 requirement)
- **Outputs:** `temp/prsice_pheno_y1.txt`, `temp/prsice_pheno_y2.txt`

### Step 5: Calculate PRS and Test Association
**Script:** `05_run_prsice_y*.sh`

- Run PRSice-2 with:
  - Clumping: r² < 0.1, 250kb window
  - P-value thresholds: 5e-8 to 1.0 (9 levels)
  - Covariates: AGE, PTGENDER, APOE4, PTEDUCAT, DX_MCI, DX_AD, SITE dummies, PC1-2
  - Binary outcome: DEPR1
- **Outputs:** `temp/prsice/nprs_depr1_y*.*` (summary, scores, plots)

## Running the Pipeline
```bash
# Load environment variables
source .env

module load R plink/6-alpha prsice

# Step 1
Rscript code/01_create_adni1_phenotype_file.R

# Step 2
python code/02_qc_gwas_sumstats.py

# Step 3
bash code/03_qc_genotype_pca.sh

# Step 4
Rscript code/04_prepare_prsice_phenotype.R

# Step 5
bash code/05_prsice_y1.sh
bash code/05_prsice_y2.sh
```

## Key Results
Key results are also copied at `report/` for easy access.

## Notes
- Analysis uses European ancestry samples only (pre-QC'd)
- Year 1 (12-month) data is primary analysis; Year 2 available for replication
- Depression prevalence may be low in ADNI (cognitively normal/MCI enriched cohort)