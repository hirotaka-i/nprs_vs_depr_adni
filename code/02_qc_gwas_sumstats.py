# 02_qc_gwas_sumstats.py
# QC for Neuroticism GWAS summary statistics for PRSice-2

import pandas as pd
import numpy as np
from scipy import stats
import gzip
import os

# Load environment variables from .env file
if os.path.exists(".env"):
    with open(".env") as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                key, value = line.split("=", 1)
                os.environ[key.strip()] = value.strip().strip('"')

# Read paths from environment variables
input_file = os.getenv("SUMSTATS")
output_file = "temp/neuroticism_gwas_qc.txt.gz"

print("=" * 80)
print("QC for Neuroticism GWAS Summary Statistics")
print("=" * 80)

# Read the data
print("\nReading summary statistics...")
df = pd.read_csv(input_file, sep='\t', compression='gzip')
print(f"Initial number of SNPs: {len(df):,}")

# Select relevant columns
print("\nSelecting relevant columns...")
df = df[['RSID', 'CHR', 'POS', 'A1', 'A2', 'EAF_UKB', 'MAF_UKB', 'Z', 'P', 'N', 'INFO_UKB']]
print(f"Columns selected: {list(df.columns)}")

# Remove rows with missing values
print("\nRemoving SNPs with missing values...")
initial_count = len(df)
df = df.dropna()
print(f"Removed {initial_count - len(df):,} SNPs with missing values")
print(f"Remaining SNPs: {len(df):,}")

# Keep only standard chromosomes (1-22, X, Y)
print("\nKeeping only standard chromosomes (chr1-chr22, chrX, chrY)...")
initial_count = len(df)
# Remove 'chr' prefix for comparison and filtering
df['CHR_clean'] = df['CHR'].astype(str).str.replace('chr', '')
valid_chrs = [str(i) for i in range(1, 23)] + ['X', 'Y']
df = df[df['CHR_clean'].isin(valid_chrs)]
# Convert chromosome to integer for 1-22, keep X and Y as strings
def convert_chr(x):
    if x in ['X', 'Y']:
        return x
    else:
        return int(x)
df['CHR'] = df['CHR_clean'].apply(convert_chr)
df = df.drop('CHR_clean', axis=1)
print(f"Removed {initial_count - len(df):,} SNPs on non-standard chromosomes")
print(f"Remaining SNPs: {len(df):,}")

# Convert POS to integer (important for PRSice-2)
df['POS'] = df['POS'].astype(int)
df['N'] = df['N'].astype(int)

# Filter by INFO score
print("\nFiltering by INFO > 0.8...")
initial_count = len(df)
df = df[df['INFO_UKB'] > 0.8]
print(f"Removed {initial_count - len(df):,} SNPs with INFO <= 0.8")
print(f"Remaining SNPs: {len(df):,}")

# Filter by MAF
print("\nFiltering by MAF > 0.01...")
initial_count = len(df)
df = df[df['MAF_UKB'] > 0.01]
print(f"Removed {initial_count - len(df):,} SNPs with MAF <= 0.01")
print(f"Remaining SNPs: {len(df):,}")

# Remove ambiguous SNPs (A/T and G/C)
print("\nRemoving strand-ambiguous SNPs (A/T and G/C)...")
initial_count = len(df)
ambiguous = ((df['A1'] == 'A') & (df['A2'] == 'T')) | \
            ((df['A1'] == 'T') & (df['A2'] == 'A')) | \
            ((df['A1'] == 'G') & (df['A2'] == 'C')) | \
            ((df['A1'] == 'C') & (df['A2'] == 'G'))
df = df[~ambiguous]
print(f"Removed {initial_count - len(df):,} ambiguous SNPs")
print(f"Remaining SNPs: {len(df):,}")

# Check for valid P-values
print("\nChecking P-value range...")
initial_count = len(df)
df = df[(df['P'] > 0) & (df['P'] <= 1)]
print(f"Removed {initial_count - len(df):,} SNPs with invalid P-values")
print(f"Remaining SNPs: {len(df):,}")

# Check for extreme Z-scores
print("\nChecking for extreme Z-scores (|Z| < 50)...")
initial_count = len(df)
df = df[np.abs(df['Z']) < 50]
print(f"Removed {initial_count - len(df):,} SNPs with extreme Z-scores")
print(f"Remaining SNPs: {len(df):,}")

# Calculate BETA from Z-score
print("\nCalculating BETA from Z-score...")
# BETA = Z / sqrt(2 * EAF * (1 - EAF) * (N + Z^2))
df['BETA'] = df['Z'] / np.sqrt(2 * df['EAF_UKB'] * (1 - df['EAF_UKB']) * (df['N'] + df['Z']**2))

# Calculate SE (standard error)
print("Calculating SE (standard error)...")
df['SE'] = 1 / np.sqrt(2 * df['EAF_UKB'] * (1 - df['EAF_UKB']) * (df['N'] + df['Z']**2))

# Check for extreme BETA values
print("\nChecking for extreme BETA values (|BETA| < 10)...")
initial_count = len(df)
df = df[np.abs(df['BETA']) < 10]
print(f"Removed {initial_count - len(df):,} SNPs with extreme BETA values")
print(f"Remaining SNPs: {len(df):,}")

# Remove duplicates (keep first occurrence)
print("\nRemoving duplicate RSIDs...")
initial_count = len(df)
df = df.drop_duplicates(subset='RSID', keep='first')
print(f"Removed {initial_count - len(df):,} duplicate SNPs")
print(f"Remaining SNPs: {len(df):,}")

# Also check for duplicate positions (same CHR:POS but different RSID)
print("\nRemoving duplicate positions (CHR:POS)...")
initial_count = len(df)
df['CHR_POS'] = df['CHR'].astype(str) + ':' + df['POS'].astype(str)
df = df.drop_duplicates(subset='CHR_POS', keep='first')
df = df.drop('CHR_POS', axis=1)
print(f"Removed {initial_count - len(df):,} SNPs with duplicate positions")
print(f"Remaining SNPs: {len(df):,}")

# Prepare final output with PRSice-2 compatible columns
print("\nPreparing final output...")
output_df = df[['RSID', 'CHR', 'POS', 'A1', 'A2', 'BETA', 'SE', 'P', 'MAF_UKB', 'INFO_UKB', 'EAF_UKB']].copy()

# Rename columns for PRSice-2 compatibility
output_df.columns = ['SNP', 'CHR', 'BP', 'A1', 'A2', 'BETA', 'SE', 'P', 'MAF', 'INFO', 'A1_FREQ']

# Ensure correct data types for PRSice-2
output_df['BP'] = output_df['BP'].astype(int)
output_df['CHR'] = output_df['CHR'].astype(str)  # Convert all to string for consistency

# Summary statistics
print("\n" + "=" * 80)
print("SUMMARY STATISTICS")
print("=" * 80)
print(f"Final number of SNPs: {len(output_df):,}")
print(f"\nBETA statistics:")
print(f"  Mean: {output_df['BETA'].mean():.6f}")
print(f"  Median: {output_df['BETA'].median():.6f}")
print(f"  SD: {output_df['BETA'].std():.6f}")
print(f"  Min: {output_df['BETA'].min():.6f}")
print(f"  Max: {output_df['BETA'].max():.6f}")
print(f"\nP-value statistics:")
print(f"  Min: {output_df['P'].min():.2e}")
print(f"  Max: {output_df['P'].max():.2e}")
print(f"  SNPs with P < 5e-8: {(output_df['P'] < 5e-8).sum():,}")
print(f"\nMAF statistics:")
print(f"  Mean: {output_df['MAF'].mean():.4f}")
print(f"  Min: {output_df['MAF'].min():.4f}")
print(f"  Max: {output_df['MAF'].max():.4f}")
print(f"\nChromosome distribution:")
chr_dist = output_df['CHR'].value_counts().sort_index()
for chr_num, count in chr_dist.items():
    print(f"  Chr {chr_num}: {count:,}")

# Save output - specify format to ensure integers are saved correctly
print(f"\nSaving QC'd summary statistics to: {output_file}")
output_df.to_csv(output_file, sep='\t', index=False, compression='gzip', 
                 float_format='%.10g')  # Use general format to avoid .0 for integers

print("\n" + "=" * 80)
print("QC COMPLETE!")
print("=" * 80)
print(f"\nOutput file: {output_file}")
print(f"Number of SNPs: {len(output_df):,}")
print("\nColumn names in output:")
print(list(output_df.columns))