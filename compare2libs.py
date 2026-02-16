#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spectral Library Comparison Tool
=================================

This script compares two spectral libraries by calculating Pearson correlation
coefficients (PCC) between matching peptide spectra. This is useful for:
- Comparing predicted vs experimental libraries
- Comparing different prediction models
- Assessing library quality

The script:
1. Matches peptides between two libraries based on sequence and charge
2. Matches fragment ions within a specified m/z tolerance
3. Calculates PCC for matched intensities
4. Generates distribution plots and summary statistics

Author: lafields2
Created: Tue Apr 15 14:55:48 2025
"""

import os
import csv
import pandas as pd
import numpy as np
from scipy.stats import pearsonr
import matplotlib.pyplot as plt
import seaborn as sns

# ============================================================================
# Configuration
# ============================================================================

# Library 1 (e.g., predicted library)
LIB1_NAME = 'UniSpec QE'
LIB1_PATH = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_building\Round1_NoModifications\Koina\UniSpec_QE\predicted_library.csv"

# Library 2 (e.g., experimental library)
LIB2_NAME = 'Exp'
LIB2_PATH = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\real_library\formatted_normalized.csv"

# Output directory
OUTPUT_FOLDER = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_building\Round1_NoModifications\models_comparison\exp_compare"

# Peak matching tolerance (m/z units)
TOLERANCE = 0.1

# Plot styling
PLOT_DPI = 300


# ============================================================================
# Core Functions
# ============================================================================

def calculate_pearson_correlation(lib1, lib2, lib1_mz_col, lib1_int_col,
                                  lib2_mz_col, lib2_int_col, tolerance):
    """
    Calculate Pearson correlation coefficient between matched peaks.
    
    This function matches fragment ions between two spectra based on m/z values
    within a specified tolerance, then computes the correlation of their intensities.
    
    Parameters
    ----------
    lib1 : pd.DataFrame
        First spectrum with m/z and intensity columns
    lib2 : pd.DataFrame
        Second spectrum with m/z and intensity columns
    lib1_mz_col : str
        Column name for m/z in first spectrum
    lib1_int_col : str
        Column name for intensity in first spectrum
    lib2_mz_col : str
        Column name for m/z in second spectrum
    lib2_int_col : str
        Column name for intensity in second spectrum
    tolerance : float
        m/z tolerance for matching peaks (e.g., 0.1 m/z units)
        
    Returns
    -------
    float
        Pearson correlation coefficient of matched intensities.
        Returns np.nan if no peaks matched or correlation cannot be computed.
        
    Notes
    -----
    The function uses a simple nearest-neighbor matching approach:
    For each peak in lib1, it finds peaks in lib2 within ±tolerance
    and uses the first match if available.
    """
    matched_lib1_intensities = []
    matched_lib2_intensities = []
    
    # Iterate through first library's peaks
    for _, lib1_peak in lib1.iterrows():
        lib1_mz = lib1_peak[lib1_mz_col]
        lib1_intensity = lib1_peak[lib1_int_col]
        
        # Find matching peaks in second library within tolerance
        matches = lib2[
            (lib2[lib2_mz_col] >= lib1_mz - tolerance) &
            (lib2[lib2_mz_col] <= lib1_mz + tolerance)
        ]
        
        # If a match is found, use the first match
        if not matches.empty:
            matched_lib1_intensities.append(lib1_intensity)
            matched_lib2_intensities.append(matches[lib2_int_col].iloc[0])
    
    # Convert to pandas Series for correlation calculation
    matched_lib1_intensities = pd.Series(matched_lib1_intensities, dtype=float)
    matched_lib2_intensities = pd.Series(matched_lib2_intensities, dtype=float)
    
    # Return NaN if no peaks matched
    if len(matched_lib1_intensities) == 0:
        return np.nan
    
    # Calculate Pearson correlation coefficient
    try:
        corr_coeff, _ = pearsonr(matched_lib1_intensities, matched_lib2_intensities)
    except ValueError:
        # This can occur if all intensities are identical (zero variance)
        corr_coeff = np.nan
    
    return corr_coeff


def create_peptide_key(df, peptide_col, charge_col):
    """
    Create a unique key for each peptide-charge combination.
    
    Parameters
    ----------
    df : pd.DataFrame
        Library dataframe
    peptide_col : str
        Column name for peptide sequence
    charge_col : str
        Column name for precursor charge
        
    Returns
    -------
    pd.DataFrame
        Dataframe with added 'key' column
    """
    df['key'] = df[peptide_col] + '_' + df[charge_col].astype(str)
    return df


def compare_libraries(lib1_df, lib2_df, mz_col='mz', int_col='intensities',
                     tolerance=0.1):
    """
    Compare two complete libraries and calculate PCC for all matching peptides.
    
    Parameters
    ----------
    lib1_df : pd.DataFrame
        First library with peptide sequences, charges, m/z, and intensities
    lib2_df : pd.DataFrame
        Second library with peptide sequences, charges, m/z, and intensities
    mz_col : str
        Column name for fragment m/z values
    int_col : str
        Column name for fragment intensities
    tolerance : float
        m/z tolerance for peak matching
        
    Returns
    -------
    pd.DataFrame
        Summary with peptide keys and their PCC values
    """
    # Get unique peptide-charge combinations from both libraries
    all_keys = set(lib1_df['key'].values.tolist() + lib2_df['key'].values.tolist())
    
    key_storage = []
    pcc_storage = []
    
    print(f"Comparing {len(all_keys)} unique peptide-charge combinations...")
    
    # Calculate PCC for each peptide
    for i, key in enumerate(all_keys):
        if (i + 1) % 100 == 0:
            print(f"  Processed {i + 1}/{len(all_keys)} peptides...")
        
        # Filter both libraries for this peptide
        lib1_filtered = lib1_df[lib1_df['key'] == key]
        lib2_filtered = lib2_df[lib2_df['key'] == key]
        
        # Skip if peptide not in both libraries
        if lib1_filtered.empty or lib2_filtered.empty:
            continue
        
        # Calculate correlation
        pcc = calculate_pearson_correlation(
            lib1=lib1_filtered,
            lib2=lib2_filtered,
            lib1_mz_col=mz_col,
            lib1_int_col=int_col,
            lib2_mz_col=mz_col,
            lib2_int_col=int_col,
            tolerance=tolerance
        )
        
        # Round to 2 decimal places
        pcc_round = round(pcc, 2) if not np.isnan(pcc) else pcc
        
        pcc_storage.append(pcc_round)
        key_storage.append(key)
    
    # Create results dataframe
    results_df = pd.DataFrame({
        'key': key_storage,
        'pcc': pcc_storage
    })
    
    return results_df


def plot_pcc_distribution(results_df, lib1_name, lib2_name, tolerance,
                          output_folder, xlim=(0, 1.2)):
    """
    Generate histogram of PCC distribution.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Results with PCC values
    lib1_name : str
        Name of first library
    lib2_name : str
        Name of second library
    tolerance : float
        m/z tolerance used
    output_folder : str
        Directory for saving plot
    xlim : tuple
        x-axis limits (min, max)
        
    Returns
    -------
    str
        Path to saved plot
    """
    # Set plot style
    sns.set(style='whitegrid', context='talk')
    
    # Create figure
    plt.figure(figsize=(10, 10))
    ax = sns.histplot(
        results_df['pcc'],
        bins=50,
        kde=True,
        color='skyblue',
        edgecolor='black'
    )
    
    # Add title and labels
    plt.title(
        f'Distribution of Pearson Correlation Coefficients\n'
        f'{lib1_name} vs. {lib2_name}\n'
        f'Tolerance = {tolerance} m/z',
        fontsize=18,
        pad=15
    )
    plt.xlabel('Pearson Correlation Coefficient', fontsize=14)
    plt.ylabel('Count', fontsize=14)
    
    # Set x-axis limits
    plt.xlim(xlim)
    
    # Generate filename
    xlim_suffix = '' if xlim == (0, 1.2) else '_extended'
    plot_filename = f'{lib1_name}_v_{lib2_name}_tol{str(tolerance)}{xlim_suffix}.png'
    plot_path = os.path.join(output_folder, plot_filename)
    
    # Save figure
    plt.tight_layout()
    plt.savefig(plot_path, dpi=PLOT_DPI)
    plt.close()
    
    return plot_path


def print_summary_statistics(results_df):
    """
    Print summary statistics for PCC distribution.
    
    Parameters
    ----------
    results_df : pd.DataFrame
        Results with PCC values
    """
    # Remove NaN values for statistics
    pcc_values = results_df['pcc'].dropna()
    
    print(f"\n{'='*60}")
    print("Summary Statistics")
    print(f"{'='*60}")
    print(f"Total comparisons: {len(results_df)}")
    print(f"Valid correlations: {len(pcc_values)}")
    print(f"Failed correlations: {len(results_df) - len(pcc_values)}")
    print(f"\nPCC Distribution:")
    print(f"  Mean: {pcc_values.mean():.3f}")
    print(f"  Median: {pcc_values.median():.3f}")
    print(f"  Std Dev: {pcc_values.std():.3f}")
    print(f"  Min: {pcc_values.min():.3f}")
    print(f"  Max: {pcc_values.max():.3f}")
    print(f"\nPercentiles:")
    print(f"  25th: {pcc_values.quantile(0.25):.3f}")
    print(f"  50th: {pcc_values.quantile(0.50):.3f}")
    print(f"  75th: {pcc_values.quantile(0.75):.3f}")
    print(f"  90th: {pcc_values.quantile(0.90):.3f}")
    print(f"{'='*60}\n")


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    """
    Main execution block.
    
    Compares two spectral libraries and generates correlation analysis.
    """
    
    print("="*60)
    print("Spectral Library Comparison Tool")
    print("="*60)
    print(f"Library 1: {LIB1_NAME}")
    print(f"  Path: {LIB1_PATH}")
    print(f"Library 2: {LIB2_NAME}")
    print(f"  Path: {LIB2_PATH}")
    print(f"m/z tolerance: {TOLERANCE}")
    print(f"Output folder: {OUTPUT_FOLDER}")
    print("="*60)
    
    # Create output folder if needed
    if not os.path.exists(OUTPUT_FOLDER):
        os.makedirs(OUTPUT_FOLDER, exist_ok=True)
        print(f"Created output directory: {OUTPUT_FOLDER}")
    
    # Load libraries
    print("\nLoading libraries...")
    lib1 = pd.read_csv(LIB1_PATH)
    lib2 = pd.read_csv(LIB2_PATH)
    
    print(f"Library 1: {len(lib1)} entries")
    print(f"Library 2: {len(lib2)} entries")
    
    # Create peptide keys
    lib1 = create_peptide_key(lib1, 'peptide', 'precursor_charges')
    lib2 = create_peptide_key(lib2, 'peptide', 'precursor_charges')
    
    print(f"Library 1: {lib1['key'].nunique()} unique peptides")
    print(f"Library 2: {lib2['key'].nunique()} unique peptides")
    
    # Compare libraries
    results_summarized = compare_libraries(
        lib1_df=lib1,
        lib2_df=lib2,
        mz_col='mz',
        int_col='intensities',
        tolerance=TOLERANCE
    )
    
    # Print summary statistics
    print_summary_statistics(results_summarized)
    
    # Save results
    results_path = os.path.join(
        OUTPUT_FOLDER,
        f'{LIB1_NAME}_v_{LIB2_NAME}_tol{str(TOLERANCE)}.csv'
    )
    results_summarized.to_csv(results_path, index=False)
    print(f"Results saved to: {results_path}")
    
    # Generate plots
    print("\nGenerating plots...")
    
    # Standard plot (0 to 1.2)
    plot_path_1 = plot_pcc_distribution(
        results_summarized,
        LIB1_NAME,
        LIB2_NAME,
        TOLERANCE,
        OUTPUT_FOLDER,
        xlim=(0, 1.2)
    )
    print(f"Saved plot: {plot_path_1}")
    
    # Extended plot (-1.2 to 1.2)
    plot_path_2 = plot_pcc_distribution(
        results_summarized,
        LIB1_NAME,
        LIB2_NAME,
        TOLERANCE,
        OUTPUT_FOLDER,
        xlim=(-1.2, 1.2)
    )
    print(f"Saved plot: {plot_path_2}")
    
    print("\n" + "="*60)
    print("Comparison complete!")
    print("="*60)
