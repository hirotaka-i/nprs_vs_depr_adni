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

----
----

# 5-Fold Cross-Validation for PRS Optimization

### Overview
Perform hyperparameter tuning for PRSice-2 using 5-fold cross-validation to identify optimal clumping r² and p-value threshold parameters while avoiding overfitting (winner's curse).

### Input Requirements
- **Phenotype file**: `temp/<outcome>.pheno` (output from Step 4)
- **GWAS summary statistics**: `temp/neuroticism_gwas_qc.txt.gz` (from Step 2)
- **Genotype data**: `temp/adni1_qc_final.*` (from Step 3)

### Scripts

#### **cv_grid_search.sh** (Master script - Steps 1-5)
Runs the complete CV pipeline automatically.

**Usage:**
```bash
bash code/cv_grid_search.sh temp/npi_dpr_ever.pheno
```

**What it does:**
1. Creates stratified folds
2. Runs grid search on training sets
3. Evaluates on held-out test sets
4. Summarizes results with visualizations
5. Recommends optimal parameters

---

#### **Step 1: cv_create_folds.R**
Creates 5 stratified folds balanced by phenotype, sex, baseline diagnosis, and age.

**Outputs:**
- `fold_assignments.csv` - Subject-to-fold mapping
- `train_fold_k.keep` - Training set IDs (k=1..5)
- `test_fold_k.keep` - Test set IDs (k=1..5)

**Stratification variables:**
- PHENO (case/control)
- PTGENDER (sex)
- DX.bl (baseline diagnosis: CN/MCI/AD)
- AGE (quartiles)

---

#### **Step 2-3: cv_train_grid.sh**
Grid search over hyperparameters on training sets.

**Parameter grid:**
- **Clumping r²**: 0.01, 0.05, 0.1, 0.2, 0.5
- **P-value thresholds**: 5e-8, 1e-6, 1e-5, 1e-4, 1e-3, 5e-3, 1e-2, 5e-2, 1e-1, 0.5, 1

**Total runs:** 5 folds × 5 r² values = 25 PRSice runs

**Outputs:**
- `training/fold_k_r2_*.prsice` - Training results for each parameter combination
- `training/best_params_per_fold.csv` - Best parameters selected per fold

---

#### **Step 4: cv_test_best.sh**
Evaluates best parameters on held-out test sets.

**For each fold:**
- Uses parameters selected from training
- Runs PRSice once on test set
- Records R², P-value, number of SNPs

**Outputs:**
- `testing/fold_k_test.*` - Test set results
- `testing/test_results.csv` - Summary of test performance

---

#### **Step 5: cv_summarize.R**
Aggregates results and generates comprehensive summary with figures.

**If all the step1-4 are done, we can rerun this with the following command:**
```bash
Rscript code/cv_summarize.R temp/npi_dpr_ever.pheno
```

**Outputs:**

**Tables:**
- `results/cv_summary.txt` - Human-readable summary
- `results/recommended_params.txt` - Parameter recommendations
- `results/cv_results_aggregated.csv` - Performance by parameter combination
- `results/cv_results_test.csv` - Test fold results
- `results/cv_results_full.csv` - All training results
- `results/forest_plot_data.csv` - Meta-analysis data
- `results/auc_by_fold.csv` - AUC per fold
- `results/quantile_summary.csv` - Prevalence by PRS quartile
- `results/quantile_regression_coefs.csv` - Adjusted associations

**Figures:**

1. **figure1_forest_plot.png**
   - OR (95% CI) per SD increase in PRS for each test fold
   - Meta-analyzed effect across folds (red)
   - Shows heterogeneity across folds

2. **figure2_roc_curve.png**
   - Pooled ROC curve from all test folds
   - Overall AUC with 95% CI
   - Reference diagonal (chance line)

3. **figure3_roc_by_fold.png**
   - Individual ROC curves for each test fold (gray)
   - Mean ROC curve (dark blue)
   - Pointwise 95% confidence band (shaded)
   - Demonstrates consistency across folds

4. **figure4_quantile_analysis.png** (3 panels)
   - **Panel A**: PRS distribution by quartile (boxplot)
   - **Panel B**: Unadjusted depression prevalence by quartile (bar chart)
   - **Panel C**: Adjusted log odds ratios (beta coefficients)
     - Adjusted for: AGE, PTGENDER, APOE4, PTEDUCAT, DX_MCI, DX_AD, PC1-4
     - Shows dose-response relationship

---

#### **Step 6: cv_final_model.sh** (User decision)
Train final model on full cohort with user-selected parameters.

**Usage:**
```bash
bash code/cv_final_model.sh temp/npi_dpr_ever.pheno <clump_r2> <p_threshold>

# Example (using recommended parameters):
bash code/cv_final_model.sh temp/npi_dpr_ever.pheno 0.1 0.005
```

**Outputs:**
- `final/full_cohort.prs.best` - Final PRS scores for all individuals
- `final/full_cohort.prs.summary` - Performance on full cohort
- `final/full_cohort.prs.snp` - SNPs included in final PRS
- `final/parameters_used.txt` - Record of parameters and rationale

**⚠️ Important:** Performance metrics from Step 6 are **optimistically biased** (trained on full data). Always report **CV test results** from Step 5 for performance claims.

---

### Output Directory Structure
```
temp/cv5/<outcome>/
├── folds/
│   ├── fold_assignments.csv
│   ├── train_fold_1.keep ... train_fold_5.keep
│   └── test_fold_1.keep ... test_fold_5.keep
├── training/
│   ├── fold_1_r2_0_01.prsice ... (all parameter combinations)
│   └── best_params_per_fold.csv
├── testing/
│   ├── fold_1_test.best ... fold_5_test.best
│   └── test_results.csv
├── results/
│   ├── cv_summary.txt
│   ├── recommended_params.txt
│   ├── cv_results_*.csv (3 files)
│   ├── forest_plot_data.csv
│   ├── auc_by_fold.csv
│   ├── quantile_summary.csv
│   ├── quantile_regression_coefs.csv
│   └── figure*.png (4 figures)
├── final/
│   ├── full_cohort.prs.best
│   ├── full_cohort.prs.summary
│   ├── full_cohort.prs.snp
│   └── parameters_used.txt
└── logs/
    └── *.log (PRSice log files)
```



### Citation

If using this CV framework, cite:
- PRSice-2: Choi & O'Reilly, GigaScience 2019

---

### Notes
- **Final PRS** from Step 6 is for external validation only