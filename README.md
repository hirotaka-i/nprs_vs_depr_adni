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

- Extract depression and neuropsychiatric phenotypes from ADNI
- Phenotype definitions:
  - **GDS5/GDS10**: Self-reported depression (Geriatric Depression Scale ≥5 or ≥10)
  - **NPI_DPR**: Clinician-rated depression/dysphoria (Neuropsychiatric Inventory)
  - **NPI_ANX**: Clinician-rated elation/euphoria (exploratory)
- Match with genetic data (fam file) by IID
- Create cross-sectional phenotype files:
  - Y1 (12-month), Y2 (24-month) for GDS
  - Y6 (72-month) for NPI cross-sectional
- Create time-to-event phenotypes (NPI_DPR_EVER, NPI_ANX_EVER):
  - Captures first occurrence if present, or last follow-up if never occurred
  - Includes `timepoint_censored` for survival analysis
- **Outputs:** 
  - `temp/gds5_y1.csv`, `temp/gds5_y2.csv` - Cross-sectional GDS at 12m/24m
  - `temp/npi_dpr_y6.csv`, `temp/npi_anx_y6.csv` - Cross-sectional NPI at 72m
  - `temp/npi_dpr_ever.csv`, `temp/npi_anx_ever.csv` - Time-to-event outcomes
  - Columns: IID, [outcome], AGE, PTGENDER, APOE4, PTEDUCAT, SITE, DX.bl (+timepoint_censored for EVER outcomes)

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
**Script:** `04_prepare_prsice_phenotype.R <phenotype_csv>`

- Merge phenotype data with PCs (PC1-10)
- Convert categorical variables:
  - PTGENDER: Male=1, Female=0
  - DX.bl: Dummy code (CN=reference, MCI, AD)
  - SITE: Dropped (too many categories)
- Rename outcome to standardized "PHENO" column
- Format: FID=0, IID, PHENO, covariates, PC1-10 (PRSice-2 requirement)
- **Outputs:** `<input_basename>.pheno` (e.g., `temp/gds5_y1.pheno`, `temp/npi_dpr_ever.pheno`)

### Step 5: Calculate PRS and Test Association
**Script:** `05_run_prsice2.sh <phenotype_file.pheno>`

- Run PRSice-2 with:
  - Clumping: r² < 0.1, 250kb window
  - P-value thresholds: 5e-8 to 1.0 (9 levels)
  - Covariates: AGE, PTGENDER, APOE4, PTEDUCAT, DX_MCI, DX_AD, PC1-PC4
  - Outcome: PHENO column (binary)
- **Outputs:** `<input_basename>.prs.*` (e.g., `temp/gds5_y1.prs.summary`, `temp/gds5_y1.prs.best`, etc.)
  - `.summary` - Association statistics and optimal P-value threshold
  - `.best` - Individual PRS at best threshold
  - `.all_score` - PRS at all thresholds
  - `.prsice` - Full results table
  - `_BARPLOT*.png`, `_HIGH-RES_PLOT*.png` - Visualizations

## Running the Pipeline
```bash
# Load environment variables
source .env

module load R plink/6-alpha prsice python/3.11

# Step 1
Rscript code/01_create_adni1_phenotype_file.R

# Step 2
python code/02_qc_gwas_sumstats.py

# Step 3
bash code/03_qc_genotype_pca.sh

# Step 4
Rscript code/04_prepare_prsice_phenotype.R "GDS5" temp/gds5_y1.csv
Rscript code/04_prepare_prsice_phenotype.R "GDS5" temp/gds5_y2.csv
Rscript code/04_prepare_prsice_phenotype.R "NPI_DPR" temp/npi_dpr_y6.csv
Rscript code/04_prepare_prsice_phenotype.R "NPI_ANX" temp/npi_anx_y6.csv
Rscript code/04_prepare_prsice_phenotype.R "NPI_DPR_EVER" temp/npi_dpr_ever.csv
Rscript code/04_prepare_prsice_phenotype.R "NPI_ANX_EVER" temp/npi_anx_ever.csv


# Step 5
bash code/05_run_prsice2.sh temp/gds5_y1.pheno
bash code/05_run_prsice2.sh temp/gds5_y2.pheno
bash code/05_run_prsice2.sh temp/npi_dpr_y6.pheno
bash code/05_run_prsice2.sh temp/npi_anx_y6.pheno
bash code/05_run_prsice2.sh temp/npi_dpr_ever.pheno
bash code/05_run_prsice2.sh temp/npi_anx_ever.pheno

# Step 6: KM plots for time-to-event outcomes
python code/06_km_curves_ps_anx.py
python code/06_km_curves_ps_dpr.py

# Copy key results to the report folder
cp temp/*.prs.log report/
cp temp/*.prs.prsice report/
cp temp/*.prs.summary report/
cp temp/*.prs*.png report/

```

## Key Results
Key results are also copied at `report/` for easy access.

## Notes
- Analysis uses European ancestry samples only (pre-QC'd)
- Year 1 (12-month) data is primary analysis; Year 2 available for replication
- Depression prevalence may be low in ADNI (cognitively normal/MCI enriched cohort)