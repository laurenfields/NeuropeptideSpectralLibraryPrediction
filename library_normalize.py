#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Spectral Library Intensity Normalization
=========================================

This script normalizes fragment ion intensities in spectral libraries.
Normalization is performed on a per-peptide basis, scaling intensities
relative to the maximum intensity for each unique peptide-charge combination.

This is important for:
- Comparing spectra across different experiments
- Ensuring consistent intensity scales
- Preparing libraries for downstream analysis tools

Author: lafields2
Created: Wed Apr 16 12:04:00 2025
"""

import os
import csv
import pandas as pd

# ============================================================================
# Configuration
# ============================================================================

# Input library path
NON_NORM_PATH = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\real_library\formatted_non_normalized.csv"

# Output library path
OUTPUT_PATH = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\real_library\formatted_normalized.csv"


# ============================================================================
# Core Functions
# ============================================================================

def normalize_intensities(df, peptide_col, charge_col, intensity_col):
    """
    Normalize fragment intensities per peptide-charge combination.
    
    For each unique peptide with a given charge state, this function:
    1. Finds the maximum fragment intensity
    2. Divides all fragment intensities by this maximum
    3. Results in normalized intensities ranging from 0 to 1
    
    Parameters
    ----------
    df : pd.DataFrame
        Input library dataframe
    peptide_col : str
        Column name containing peptide sequences (with modifications)
    charge_col : str
        Column name containing precursor charge states
    intensity_col : str
        Column name containing fragment intensities to normalize
        
    Returns
    -------
    pd.DataFrame
        Dataframe with added 'NormalizedIntensity' column
        
    Examples
    --------
    >>> library = pd.DataFrame({
    ...     'peptide': ['PEPTIDE', 'PEPTIDE', 'PEPTIDE'],
    ...     'charge': [2, 2, 2],
    ...     'intensity': [1000, 500, 750]
    ... })
    >>> normalized = normalize_intensities(library, 'peptide', 'charge', 'intensity')
    >>> normalized['NormalizedIntensity'].tolist()
    [1.0, 0.5, 0.75]
    """
    # Group by peptide sequence and charge
    group_keys = [peptide_col, charge_col]
    
    # Apply normalization: divide each intensity by the max in its group
    df['NormalizedIntensity'] = df.groupby(group_keys)[intensity_col].transform(
        lambda x: x / x.max()
    )
    
    return df


def validate_normalization(df):
    """
    Validate that normalization was successful.
    
    Checks:
    1. All normalized intensities are between 0 and 1
    2. Each peptide-charge group has at least one intensity of 1.0
    3. No NaN values were introduced
    
    Parameters
    ----------
    df : pd.DataFrame
        Normalized library dataframe
        
    Returns
    -------
    bool
        True if validation passes, False otherwise
    """
    # Check for NaN values
    if df['NormalizedIntensity'].isna().any():
        print("Warning: NaN values found in normalized intensities")
        return False
    
    # Check value range
    min_val = df['NormalizedIntensity'].min()
    max_val = df['NormalizedIntensity'].max()
    
    if min_val < 0 or max_val > 1:
        print(f"Warning: Normalized intensities out of range [0,1]: min={min_val}, max={max_val}")
        return False
    
    # Check that max is 1.0 (accounting for floating point precision)
    if abs(max_val - 1.0) > 1e-6:
        print(f"Warning: Maximum normalized intensity is {max_val}, expected 1.0")
        return False
    
    print("Normalization validation passed")
    return True


def print_normalization_summary(original_df, normalized_df):
    """
    Print summary statistics comparing original and normalized intensities.
    
    Parameters
    ----------
    original_df : pd.DataFrame
        Original library before normalization
    normalized_df : pd.DataFrame
        Library after normalization
    """
    print(f"\n{'='*60}")
    print("Normalization Summary")
    print(f"{'='*60}")
    print(f"Total entries: {len(normalized_df)}")
    print(f"Unique peptides: {normalized_df['peptide_sequences'].nunique()}")
    print(f"Unique peptide-charge combinations: {normalized_df.groupby(['peptide_sequences', 'precursor_charges']).ngroups}")
    
    print(f"\nOriginal Intensity Statistics:")
    print(f"  Min: {original_df['intensities'].min():.2e}")
    print(f"  Max: {original_df['intensities'].max():.2e}")
    print(f"  Mean: {original_df['intensities'].mean():.2e}")
    print(f"  Median: {original_df['intensities'].median():.2e}")
    
    print(f"\nNormalized Intensity Statistics:")
    print(f"  Min: {normalized_df['NormalizedIntensity'].min():.4f}")
    print(f"  Max: {normalized_df['NormalizedIntensity'].max():.4f}")
    print(f"  Mean: {normalized_df['NormalizedIntensity'].mean():.4f}")
    print(f"  Median: {normalized_df['NormalizedIntensity'].median():.4f}")
    print(f"{'='*60}\n")


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    """
    Main execution block.
    
    Loads a spectral library, normalizes intensities, validates the result,
    and saves the normalized library.
    """
    
    print("="*60)
    print("Spectral Library Intensity Normalization")
    print("="*60)
    print(f"Input: {NON_NORM_PATH}")
    print(f"Output: {OUTPUT_PATH}")
    print("="*60)
    
    # Load non-normalized library
    print("\nLoading library...")
    non_norm = pd.read_csv(NON_NORM_PATH)
    print(f"Loaded {len(non_norm)} entries")
    
    # Check required columns
    required_cols = ['peptide_sequences', 'precursor_charges', 'intensities']
    missing_cols = [col for col in required_cols if col not in non_norm.columns]
    
    if missing_cols:
        print(f"Error: Missing required columns: {missing_cols}")
        print(f"Available columns: {non_norm.columns.tolist()}")
        exit(1)
    
    # Perform normalization
    print("Normalizing intensities...")
    norm_lib = normalize_intensities(
        df=non_norm,
        peptide_col='peptide_sequences',
        charge_col='precursor_charges',
        intensity_col='intensities'
    )
    
    # Validate normalization
    validation_passed = validate_normalization(norm_lib)
    
    if not validation_passed:
        print("Warning: Normalization validation failed. Please review the data.")
    
    # Print summary
    print_normalization_summary(non_norm, norm_lib)
    
    # Create output directory if needed
    output_dir = os.path.dirname(OUTPUT_PATH)
    if output_dir and not os.path.exists(output_dir):
        os.makedirs(output_dir, exist_ok=True)
        print(f"Created output directory: {output_dir}")
    
    # Save normalized library
    print(f"Saving normalized library to: {OUTPUT_PATH}")
    norm_lib.to_csv(OUTPUT_PATH, index=False)
    
    print("\n" + "="*60)
    print("Normalization complete!")
    print("="*60)
