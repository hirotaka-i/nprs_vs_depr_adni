#!/usr/bin/env python3
"""
Kaplan-Meier Survival Analysis for Depression based on PRS, stratified by baseline diagnosis
This script creates KM curves comparing high vs low PRS groups, separately for CN, MCI, and AD.
"""

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from lifelines import KaplanMeierFitter
from lifelines.statistics import logrank_test, multivariate_logrank_test
import seaborn as sns

# Set style
sns.set_style("whitegrid")
plt.rcParams['figure.figsize'] = (14, 10)

# Load PRS data
print("Loading PRS data...")
prs_data = pd.read_csv('temp/npi_dpr_ever.prs.best', sep=r'\s+')
print(f"PRS data shape: {prs_data.shape}")
print(f"PRS columns: {prs_data.columns.tolist()}")

# Load phenotype data
print("\nLoading phenotype data...")
pheno_data = pd.read_csv('temp/npi_dpr_ever.csv')
print(f"Phenotype data shape: {pheno_data.shape}")
print(f"Phenotype columns: {pheno_data.columns.tolist()}")

# Merge datasets on IID
print("\nMerging datasets...")
merged_data = pheno_data.merge(prs_data[['IID', 'PRS']], on='IID', how='inner')
print(f"Merged data shape: {merged_data.shape}")

# Check for missing values
print("\nMissing values:")
print(merged_data[['timepoint_censored', 'NPI_DPR_EVER', 'PRS', 'DX.bl']].isnull().sum())

# Remove rows with missing values in key variables
merged_data = merged_data.dropna(subset=['timepoint_censored', 'NPI_DPR_EVER', 'PRS', 'DX.bl'])
print(f"Data shape after removing missing values: {merged_data.shape}")

# Print baseline diagnosis distribution
print("\nBaseline diagnosis distribution:")
print(merged_data['DX.bl'].value_counts())

# Create PRS groups (high vs low based on median split)
prs_median = merged_data['PRS'].median()
merged_data['PRS_group'] = merged_data['PRS'].apply(
    lambda x: 'High PRS' if x >= prs_median else 'Low PRS'
)

print(f"\nPRS median: {prs_median:.6f}")
print(f"PRS range: [{merged_data['PRS'].min():.6f}, {merged_data['PRS'].max():.6f}]")

# Create combined group variable for plotting
merged_data['Strata'] = merged_data['DX.bl'] + ' - ' + merged_data['PRS_group']

print("\nStrata distribution:")
print(merged_data['Strata'].value_counts().sort_index())

# Print summary statistics
print("\nEvent status by DX.bl and PRS group:")
summary_table = pd.crosstab(
    [merged_data['DX.bl'], merged_data['PRS_group']], 
    merged_data['NPI_DPR_EVER'], 
    margins=True, 
    margins_name='Total'
)
print(summary_table)

# =============================================================================
# Plot 1: All strata on one plot
# =============================================================================
print("\nCreating combined KM plot...")
fig, ax = plt.subplots(figsize=(14, 10))

kmf = KaplanMeierFitter()

# Define colors for each diagnosis
dx_colors = {'CN': ['#2E86AB', '#A23B72'], 'MCI': ['#F18F01', '#C73E1D'], 'AD': ['#6A994E', '#BC4B51']}
dx_order = ['CN', 'MCI', 'AD']

for dx in dx_order:
    if dx not in merged_data['DX.bl'].values:
        continue
    
    for i, prs_group in enumerate(['Low PRS', 'High PRS']):
        mask = (merged_data['DX.bl'] == dx) & (merged_data['PRS_group'] == prs_group)
        
        if mask.sum() > 0:
            kmf.fit(
                durations=merged_data.loc[mask, 'timepoint_censored'],
                event_observed=merged_data.loc[mask, 'NPI_DPR_EVER'],  # CHANGED
                label=f'{dx} - {prs_group}'
            )
            
            linestyle = '--' if prs_group == 'High PRS' else '-'
            kmf.plot_survival_function(
                ax=ax, 
                ci_show=True, 
                linewidth=2.5,
                linestyle=linestyle,
                color=dx_colors[dx][i]
            )

ax.set_xlabel('Time (months)', fontsize=13)
ax.set_ylabel('Depression-Free Survival Probability', fontsize=13)  # CHANGED
ax.set_title('Kaplan-Meier Curves: Depression Risk by PRS and Baseline Diagnosis',  # CHANGED
             fontsize=15, fontweight='bold', pad=20)
ax.legend(loc='best', fontsize=10, framealpha=0.9)
ax.grid(True, alpha=0.3)
plt.tight_layout()

output_file_combined = 'report/km_curve_npi_dpr_ever_prs_stratified.png'  # CHANGED
plt.savefig(output_file_combined, dpi=300, bbox_inches='tight')
print(f"Combined KM curve saved to: {output_file_combined}")

# =============================================================================
# Plot 2: Separate subplots for each diagnosis
# =============================================================================
print("\nCreating separate KM plots by diagnosis...")
fig, axes = plt.subplots(1, 3, figsize=(18, 6), sharey=True)

for idx, dx in enumerate(dx_order):
    ax = axes[idx]
    
    if dx not in merged_data['DX.bl'].values:
        ax.text(0.5, 0.5, f'No {dx} data', ha='center', va='center', fontsize=14)
        ax.set_title(dx, fontsize=14, fontweight='bold')
        continue
    
    kmf_low = KaplanMeierFitter()
    kmf_high = KaplanMeierFitter()
    
    # Low PRS
    mask_low = (merged_data['DX.bl'] == dx) & (merged_data['PRS_group'] == 'Low PRS')
    if mask_low.sum() > 0:
        kmf_low.fit(
            durations=merged_data.loc[mask_low, 'timepoint_censored'],
            event_observed=merged_data.loc[mask_low, 'NPI_DPR_EVER'],  # CHANGED
            label='Low PRS'
        )
        kmf_low.plot_survival_function(ax=ax, ci_show=True, linewidth=2.5, color='#2E86AB')
    
    # High PRS
    mask_high = (merged_data['DX.bl'] == dx) & (merged_data['PRS_group'] == 'High PRS')
    if mask_high.sum() > 0:
        kmf_high.fit(
            durations=merged_data.loc[mask_high, 'timepoint_censored'],
            event_observed=merged_data.loc[mask_high, 'NPI_DPR_EVER'],  # CHANGED
            label='High PRS'
        )
        kmf_high.plot_survival_function(ax=ax, ci_show=True, linewidth=2.5, color='#C73E1D', linestyle='--')
    
    # Log-rank test for this diagnosis
    if mask_low.sum() > 0 and mask_high.sum() > 0:
        lr_result = logrank_test(
            durations_A=merged_data.loc[mask_low, 'timepoint_censored'],
            durations_B=merged_data.loc[mask_high, 'timepoint_censored'],
            event_observed_A=merged_data.loc[mask_low, 'NPI_DPR_EVER'],  # CHANGED
            event_observed_B=merged_data.loc[mask_high, 'NPI_DPR_EVER']   # CHANGED
        )
        
        # Add p-value to plot
        ax.text(0.05, 0.05, f'Log-rank p = {lr_result.p_value:.4f}', 
                transform=ax.transAxes, fontsize=11, 
                bbox=dict(boxstyle='round', facecolor='wheat', alpha=0.5))
    
    ax.set_xlabel('Time (months)', fontsize=12)
    if idx == 0:
        ax.set_ylabel('Depression-Free Survival Probability', fontsize=12)  # CHANGED
    ax.set_title(f'{dx} (n={mask_low.sum() + mask_high.sum()})', fontsize=14, fontweight='bold')
    ax.legend(loc='best', fontsize=10)
    ax.grid(True, alpha=0.3)

plt.tight_layout()

output_file_separate = 'report/km_curve_npi_dpr_ever_prs_by_diagnosis.png'  # CHANGED
plt.savefig(output_file_separate, dpi=300, bbox_inches='tight')
print(f"Separate KM curves saved to: {output_file_separate}")

# =============================================================================
# Statistical Tests
# =============================================================================
print("\n" + "="*70)
print("STATISTICAL TESTS")
print("="*70)

# Overall log-rank test (PRS groups, ignoring diagnosis)
print("\n1. Overall Log-Rank Test (PRS groups, all diagnoses combined):")
print("-"*70)
low_prs_mask = merged_data['PRS_group'] == 'Low PRS'
high_prs_mask = merged_data['PRS_group'] == 'High PRS'

overall_lr = logrank_test(
    durations_A=merged_data.loc[low_prs_mask, 'timepoint_censored'],
    durations_B=merged_data.loc[high_prs_mask, 'timepoint_censored'],
    event_observed_A=merged_data.loc[low_prs_mask, 'NPI_DPR_EVER'],  # CHANGED
    event_observed_B=merged_data.loc[high_prs_mask, 'NPI_DPR_EVER']   # CHANGED
)
print(overall_lr)

# Stratified log-rank tests by diagnosis
print("\n2. Log-Rank Tests Stratified by Baseline Diagnosis:")
print("-"*70)

for dx in dx_order:
    if dx not in merged_data['DX.bl'].values:
        continue
        
    mask_low = (merged_data['DX.bl'] == dx) & (merged_data['PRS_group'] == 'Low PRS')
    mask_high = (merged_data['DX.bl'] == dx) & (merged_data['PRS_group'] == 'High PRS')
    
    if mask_low.sum() > 0 and mask_high.sum() > 0:
        lr_result = logrank_test(
            durations_A=merged_data.loc[mask_low, 'timepoint_censored'],
            durations_B=merged_data.loc[mask_high, 'timepoint_censored'],
            event_observed_A=merged_data.loc[mask_low, 'NPI_DPR_EVER'],  # CHANGED
            event_observed_B=merged_data.loc[mask_high, 'NPI_DPR_EVER']   # CHANGED
        )
        
        print(f"\n{dx}:")
        print(f"  Low PRS:  n={mask_low.sum()}, events={merged_data.loc[mask_low, 'NPI_DPR_EVER'].sum()}")  # CHANGED
        print(f"  High PRS: n={mask_high.sum()}, events={merged_data.loc[mask_high, 'NPI_DPR_EVER'].sum()}")  # CHANGED
        print(f"  Test statistic: {lr_result.test_statistic:.4f}")
        print(f"  P-value: {lr_result.p_value:.4f}")
        print(f"  Significant (α=0.05): {'Yes' if lr_result.p_value < 0.05 else 'No'}")

# Multivariate log-rank test (test for interaction)
print("\n3. Multivariate Log-Rank Test (all strata):")
print("-"*70)
mv_lr = multivariate_logrank_test(
    durations=merged_data['timepoint_censored'],
    groups=merged_data['Strata'],
    event_observed=merged_data['NPI_DPR_EVER']  # CHANGED
)
print(mv_lr)

# =============================================================================
# Save Summary
# =============================================================================
summary_file = 'report/km_summary_npi_dpr_ever_prs_stratified.txt'  # CHANGED
with open(summary_file, 'w') as f:
    f.write("Kaplan-Meier Analysis Summary: Depression Risk by PRS (Stratified by DX.bl)\n")  # CHANGED
    f.write("="*70 + "\n\n")
    
    f.write(f"Total sample size: {len(merged_data)}\n")
    f.write(f"PRS median split: {prs_median:.6f}\n\n")
    
    f.write("Sample sizes by baseline diagnosis and PRS group:\n")
    f.write(summary_table.to_string())
    f.write("\n\n")
    
    f.write("Overall Log-rank test (PRS groups, all diagnoses):\n")
    f.write(f"  Test statistic: {overall_lr.test_statistic:.4f}\n")
    f.write(f"  P-value: {overall_lr.p_value:.4f}\n\n")
    
    f.write("Stratified Log-rank tests:\n")
    for dx in dx_order:
        if dx not in merged_data['DX.bl'].values:
            continue
            
        mask_low = (merged_data['DX.bl'] == dx) & (merged_data['PRS_group'] == 'Low PRS')
        mask_high = (merged_data['DX.bl'] == dx) & (merged_data['PRS_group'] == 'High PRS')
        
        if mask_low.sum() > 0 and mask_high.sum() > 0:
            lr_result = logrank_test(
                durations_A=merged_data.loc[mask_low, 'timepoint_censored'],
                durations_B=merged_data.loc[mask_high, 'timepoint_censored'],
                event_observed_A=merged_data.loc[mask_low, 'NPI_DPR_EVER'],  # CHANGED
                event_observed_B=merged_data.loc[mask_high, 'NPI_DPR_EVER']   # CHANGED
            )
            f.write(f"\n  {dx}:\n")
            f.write(f"    n(Low/High): {mask_low.sum()}/{mask_high.sum()}\n")
            f.write(f"    P-value: {lr_result.p_value:.4f}\n")
    
    f.write(f"\nMultivariate log-rank test:\n")
    f.write(f"  Test statistic: {mv_lr.test_statistic:.4f}\n")
    f.write(f"  P-value: {mv_lr.p_value:.4f}\n")

print(f"\nSummary statistics saved to: {summary_file}")
print("\nAnalysis complete!")