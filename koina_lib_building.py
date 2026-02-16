#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Koina Spectral Library Generation
==================================

This script generates predicted spectral libraries using various prediction models
available through the Koina API. It processes peptide sequences and generates
intensity predictions for different fragmentation instruments and models.

The script supports multiple prediction models:
- AlphaPeptDeep (QE and LUMOS instruments)
- Prosit variants (2019, 2020, 2024, 2025)
- UniSpec (multiple instrument types)
- ms2pip variants

Author: lafields2
Created: Tue Apr 15 11:58:04 2025
"""

import os
import re
import csv
import pandas as pd
import numpy as np
from koinapy import Koina
from more_itertools import sliced

# ============================================================================
# Configuration
# ============================================================================

# Directory paths
OUTPUT_DIR = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_building\Round2_FullDB"
LIBRARY_PATH = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_building\Round2_FullDB\crustaceanDB_wMods.csv"

# Processing parameters
CHUNK_SIZE = 300  # Number of peptides to process in each batch (for large datasets)

# Koina server configuration
KOINA_SERVER = "koina.wilhelmlab.org:443"


# ============================================================================
# Helper Functions
# ============================================================================

def filter_modifications(library_df, modifications_to_remove):
    """
    Remove peptides containing specific UNIMOD modifications.
    
    Parameters
    ----------
    library_df : pd.DataFrame
        Input library with peptide sequences
    modifications_to_remove : list
        List of UNIMOD modification strings to filter out (e.g., ['[UNIMOD:2]'])
        
    Returns
    -------
    pd.DataFrame
        Filtered library without specified modifications
    """
    pattern = '|'.join(map(re.escape, modifications_to_remove))
    filtered_df = library_df[~library_df['peptide_sequences'].str.contains(
        pattern, regex=True, na=False)].reset_index(drop=True)
    return filtered_df


def reformat_modifications(library_df, replacements):
    """
    Reformat modification notation to match model requirements.
    
    Some models require N-terminal or C-terminal modifications to be placed
    after the amino acid rather than before.
    
    Parameters
    ----------
    library_df : pd.DataFrame
        Input library with peptide sequences
    replacements : dict
        Dictionary mapping old patterns to new patterns (e.g., {r'\[UNIMOD:28\]Q': r'Q[UNIMOD:28]'})
        
    Returns
    -------
    pd.DataFrame
        Library with reformatted modification notation
    """
    for old, new in replacements.items():
        library_df['peptide_sequences'] = library_df['peptide_sequences'].str.replace(
            old, new, regex=True)
    return library_df


def save_predictions(predictions_df, output_dir, model_name, instrument_type=None):
    """
    Save prediction results to CSV file.
    
    Parameters
    ----------
    predictions_df : pd.DataFrame
        Predictions from Koina model
    output_dir : str
        Base output directory
    model_name : str
        Name of the prediction model
    instrument_type : str, optional
        Instrument type (if applicable)
        
    Returns
    -------
    str
        Path to saved file
    """
    # Create subdirectory name
    if instrument_type:
        subdir_name = f'{model_name}_{instrument_type}'
    else:
        subdir_name = model_name
        
    output_subdir = os.path.join(output_dir, subdir_name)
    
    # Create output directory if it doesn't exist
    if not os.path.exists(output_subdir):
        os.makedirs(output_subdir, exist_ok=True)
    
    # Save predictions
    output_path = os.path.join(output_subdir, 'predicted_library.csv')
    predictions_df.to_csv(output_path, index=False)
    
    print(f"Saved predictions to: {output_path}")
    return output_path


def predict_in_chunks(library_df, model_name, chunk_size=CHUNK_SIZE):
    """
    Process library in chunks to handle large datasets.
    
    Parameters
    ----------
    library_df : pd.DataFrame
        Input library for prediction
    model_name : str
        Name of the Koina model to use
    chunk_size : int
        Number of peptides to process per chunk
        
    Returns
    -------
    pd.DataFrame
        Combined predictions from all chunks
    """
    index_slices = sliced(range(len(library_df)), chunk_size)
    predictions = pd.DataFrame()
    
    for i, index_slice in enumerate(index_slices):
        print(f"Processing chunk {i+1}...")
        chunk = library_df.iloc[index_slice]
        model = Koina(model_name, KOINA_SERVER)
        predictions_chunk = model.predict(chunk, debug=True, mode='sync')
        
        # Concatenate results
        if len(predictions) > 0:
            predictions = pd.concat([predictions, predictions_chunk], ignore_index=True)
        else:
            predictions = predictions_chunk
    
    return predictions


# ============================================================================
# Model-Specific Prediction Functions
# ============================================================================

def predict_alphapeptdeep_generic(library_path, output_dir, instrument_type='QE'):
    """
    Generate predictions using AlphaPeptDeep_ms2_generic model.
    
    Limitations:
    - Maximum peptide length: 32 amino acids
    - Does not support [UNIMOD:2] modifications
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    instrument_type : str
        Instrument type ('QE' or 'LUMOS')
    """
    model_name = "AlphaPeptDeep_ms2_generic"
    print(f"\n{'='*60}")
    print(f"Predicting with {model_name} ({instrument_type})")
    print(f"{'='*60}")
    
    # Load and prepare library
    library = pd.read_csv(library_path)
    library['instrument_types'] = instrument_type
    
    # Filter out unsupported modifications
    library = filter_modifications(library, ['[UNIMOD:2]'])
    
    # Reformat modification notation
    replacements = {
        r'\[UNIMOD:28\]Q': r'Q[UNIMOD:28]',
        r'\[UNIMOD:27\]E': r'E[UNIMOD:27]'
    }
    library = reformat_modifications(library, replacements)
    
    # Apply length filter
    library = library[library['peptide'].str.len() <= 32].reset_index(drop=True)
    
    print(f"Processing {len(library)} peptides...")
    
    # Generate predictions in chunks
    predictions = predict_in_chunks(library, model_name)
    
    # Save results
    save_predictions(predictions, output_dir, model_name, instrument_type)


def predict_prosit_2019(library_path, output_dir):
    """
    Generate predictions using Prosit_2019_intensity model.
    
    Limitations:
    - Maximum peptide length: 30 amino acids
    - Maximum precursor charge: 6
    - Does not support PTMs ([UNIMOD:40], [UNIMOD:28], [UNIMOD:27], [UNIMOD:2])
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    """
    model_name = "Prosit_2019_intensity"
    print(f"\n{'='*60}")
    print(f"Predicting with {model_name}")
    print(f"{'='*60}")
    
    # Load and prepare library
    library = pd.read_csv(library_path)
    
    # Reformat modification notation
    replacements = {
        r'\[UNIMOD:28\]Q': r'Q[UNIMOD:28]',
        r'\[UNIMOD:27\]E': r'E[UNIMOD:27]'
    }
    library = reformat_modifications(library, replacements)
    
    # Filter out unsupported modifications
    modifications_to_remove = ['[UNIMOD:40]', '[UNIMOD:28]', '[UNIMOD:27]', '[UNIMOD:2]']
    library = filter_modifications(library, modifications_to_remove)
    
    # Apply length and charge filters
    library = library[library['peptide'].str.len() <= 30].reset_index(drop=True)
    library = library[library['precursor_charges'] <= 6].reset_index(drop=True)
    
    print(f"Processing {len(library)} peptides...")
    
    # Generate predictions
    model = Koina(model_name, KOINA_SERVER)
    predictions = model.predict(library)
    
    # Save results
    save_predictions(predictions, output_dir, model_name)


def predict_prosit_2020_hcd(library_path, output_dir):
    """
    Generate predictions using Prosit_2020_intensity_HCD model.
    
    Limitations:
    - Maximum peptide length: 30 amino acids
    - Maximum precursor charge: 6
    - Does not support PTMs ([UNIMOD:40], [UNIMOD:28], [UNIMOD:27], [UNIMOD:2])
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    """
    model_name = "Prosit_2020_intensity_HCD"
    print(f"\n{'='*60}")
    print(f"Predicting with {model_name}")
    print(f"{'='*60}")
    
    # Load and prepare library
    library = pd.read_csv(library_path)
    
    # Reformat modification notation
    replacements = {
        r'\[UNIMOD:28\]Q': r'Q[UNIMOD:28]',
        r'\[UNIMOD:27\]E': r'E[UNIMOD:27]'
    }
    library = reformat_modifications(library, replacements)
    
    # Filter out unsupported modifications
    modifications_to_remove = ['[UNIMOD:40]', '[UNIMOD:28]', '[UNIMOD:27]', '[UNIMOD:2]']
    library = filter_modifications(library, modifications_to_remove)
    
    # Apply length and charge filters
    library = library[library['peptide'].str.len() <= 30].reset_index(drop=True)
    library = library[library['precursor_charges'] <= 6].reset_index(drop=True)
    
    print(f"Processing {len(library)} peptides...")
    
    # Generate predictions
    model = Koina(model_name, KOINA_SERVER)
    predictions = model.predict(library)
    
    # Save results
    save_predictions(predictions, output_dir, model_name)


def predict_prosit_2024_ptms(library_path, output_dir):
    """
    Generate predictions using Prosit_2024_intensity_PTMs_gl model.
    
    This model supports post-translational modifications.
    
    Limitations:
    - Maximum peptide length: 30 amino acids
    - Maximum precursor charge: 6
    - Does not support some PTMs ([UNIMOD:40], [UNIMOD:28], [UNIMOD:27], [UNIMOD:2])
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    """
    model_name = "Prosit_2024_intensity_PTMs_gl"
    print(f"\n{'='*60}")
    print(f"Predicting with {model_name}")
    print(f"{'='*60}")
    
    # Load and prepare library
    library = pd.read_csv(library_path)
    library['fragmentation_types'] = 'HCD'
    
    # Reformat modification notation
    replacements = {
        r'\[UNIMOD:28\]Q': r'Q[UNIMOD:28]',
        r'\[UNIMOD:27\]E': r'E[UNIMOD:27]'
    }
    library = reformat_modifications(library, replacements)
    
    # Filter out unsupported modifications
    modifications_to_remove = ['[UNIMOD:40]', '[UNIMOD:28]', '[UNIMOD:27]', '[UNIMOD:2]']
    library = filter_modifications(library, modifications_to_remove)
    
    # Apply length and charge filters
    library = library[library['peptide'].str.len() <= 30].reset_index(drop=True)
    library = library[library['precursor_charges'] <= 6].reset_index(drop=True)
    
    print(f"Processing {len(library)} peptides...")
    
    # Generate predictions
    model = Koina(model_name, KOINA_SERVER)
    predictions = model.predict(library)
    
    # Save results
    save_predictions(predictions, output_dir, model_name)


def predict_prosit_2025_multifrag(library_path, output_dir):
    """
    Generate predictions using Prosit_2025_intensity_MultiFrag model.
    
    Limitations:
    - Maximum peptide length: 30 amino acids
    - Does not support some PTMs ([UNIMOD:40], [UNIMOD:28], [UNIMOD:27], [UNIMOD:2])
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    """
    model_name = "Prosit_2025_intensity_MultiFrag"
    print(f"\n{'='*60}")
    print(f"Predicting with {model_name}")
    print(f"{'='*60}")
    
    # Load and prepare library
    library = pd.read_csv(library_path)
    library['fragmentation_types'] = 'HCD'
    
    # Reformat modification notation
    replacements = {
        r'\[UNIMOD:28\]Q': r'Q[UNIMOD:28]',
        r'\[UNIMOD:27\]E': r'E[UNIMOD:27]'
    }
    library = reformat_modifications(library, replacements)
    
    # Filter out unsupported modifications
    modifications_to_remove = ['[UNIMOD:40]', '[UNIMOD:28]', '[UNIMOD:27]', '[UNIMOD:2]']
    library = filter_modifications(library, modifications_to_remove)
    
    # Apply length filter
    library = library[library['peptide'].str.len() <= 30].reset_index(drop=True)
    
    print(f"Processing {len(library)} peptides...")
    
    # Generate predictions
    model = Koina(model_name, KOINA_SERVER)
    predictions = model.predict(library)
    
    # Save results
    save_predictions(predictions, output_dir, model_name)


def predict_unispec(library_path, output_dir, instrument_type='QE'):
    """
    Generate predictions using UniSpec model.
    
    Limitations vary by instrument type:
    - QE/QEHFX: Maximum precursor charge 5, Maximum peptide length 40 AA
    - LUMOS/NONE: Maximum peptide length 40 AA
    - Does not support [UNIMOD:40] and [UNIMOD:2] modifications
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    instrument_type : str
        Instrument type ('QE', 'QEHFX', 'LUMOS', or 'NONE')
    """
    model_name = "UniSpec"
    print(f"\n{'='*60}")
    print(f"Predicting with {model_name} ({instrument_type})")
    print(f"{'='*60}")
    
    # Load and prepare library
    library = pd.read_csv(library_path)
    library['instrument_types'] = instrument_type
    
    # Reformat modification notation
    replacements = {
        r'\[UNIMOD:28\]Q': r'Q[UNIMOD:28]',
        r'\[UNIMOD:27\]E': r'E[UNIMOD:27]'
    }
    library = reformat_modifications(library, replacements)
    
    # Apply instrument-specific filters
    if instrument_type in ['QE', 'QEHFX']:
        library = library[library['precursor_charges'] <= 5]
    
    # Apply length filter
    library = library[library['peptide'].str.len() <= 40].reset_index(drop=True)
    
    # Filter out unsupported modifications
    modifications_to_remove = ['[UNIMOD:40]', '[UNIMOD:2]']
    library = filter_modifications(library, modifications_to_remove)
    
    print(f"Processing {len(library)} peptides...")
    
    # Generate predictions
    model = Koina(model_name, KOINA_SERVER)
    predictions = model.predict(library, debug=True)
    
    # Save results
    save_predictions(predictions, output_dir, model_name, instrument_type)


def predict_ms2pip_hcd2021(library_path, output_dir):
    """
    Generate predictions using ms2pip_HCD2021 model.
    
    Limitations:
    - Maximum peptide length: 30 amino acids
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    """
    model_name = "ms2pip_HCD2021"
    print(f"\n{'='*60}")
    print(f"Predicting with {model_name}")
    print(f"{'='*60}")
    
    # Load and prepare library
    library = pd.read_csv(library_path)
    
    # Reformat modification notation
    replacements = {
        r'\[UNIMOD:28\]Q': r'Q[UNIMOD:28]',
        r'\[UNIMOD:27\]E': r'E[UNIMOD:27]'
    }
    library = reformat_modifications(library, replacements)
    
    # Apply length filter
    library = library[library['peptide'].str.len() <= 30].reset_index(drop=True)
    
    print(f"Processing {len(library)} peptides...")
    
    # Generate predictions
    model = Koina(model_name, KOINA_SERVER)
    predictions = model.predict(library)
    
    # Save results
    save_predictions(predictions, output_dir, model_name)


def predict_ms2pip_immuno_hcd(library_path, output_dir):
    """
    Generate predictions using ms2pip_Immuno_HCD model.
    
    Optimized for immunopeptide analysis.
    
    Limitations:
    - Maximum peptide length: 30 amino acids
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    """
    model_name = "ms2pip_Immuno_HCD"
    print(f"\n{'='*60}")
    print(f"Predicting with {model_name}")
    print(f"{'='*60}")
    
    # Load and prepare library
    library = pd.read_csv(library_path)
    
    # Reformat modification notation
    replacements = {
        r'\[UNIMOD:28\]Q': r'Q[UNIMOD:28]',
        r'\[UNIMOD:27\]E': r'E[UNIMOD:27]'
    }
    library = reformat_modifications(library, replacements)
    
    # Apply length filter
    library = library[library['peptide'].str.len() <= 30].reset_index(drop=True)
    
    print(f"Processing {len(library)} peptides...")
    
    # Generate predictions
    model = Koina(model_name, KOINA_SERVER)
    predictions = model.predict(library)
    
    # Save results
    save_predictions(predictions, output_dir, model_name)


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    """
    Main execution block.
    
    Uncomment the desired prediction function(s) to run.
    Each function is independent and can be run separately or in sequence.
    """
    
    print("="*60)
    print("Koina Spectral Library Generation")
    print("="*60)
    print(f"Input library: {LIBRARY_PATH}")
    print(f"Output directory: {OUTPUT_DIR}")
    print("="*60)
    
    # AlphaPeptDeep predictions
    # predict_alphapeptdeep_generic(LIBRARY_PATH, OUTPUT_DIR, instrument_type='QE')
    # predict_alphapeptdeep_generic(LIBRARY_PATH, OUTPUT_DIR, instrument_type='LUMOS')
    
    # Prosit predictions
    # predict_prosit_2019(LIBRARY_PATH, OUTPUT_DIR)
    # predict_prosit_2020_hcd(LIBRARY_PATH, OUTPUT_DIR)
    # predict_prosit_2024_ptms(LIBRARY_PATH, OUTPUT_DIR)
    # predict_prosit_2025_multifrag(LIBRARY_PATH, OUTPUT_DIR)
    
    # UniSpec predictions (multiple instrument types)
    # predict_unispec(LIBRARY_PATH, OUTPUT_DIR, instrument_type='QE')
    # predict_unispec(LIBRARY_PATH, OUTPUT_DIR, instrument_type='QEHFX')
    # predict_unispec(LIBRARY_PATH, OUTPUT_DIR, instrument_type='LUMOS')
    # predict_unispec(LIBRARY_PATH, OUTPUT_DIR, instrument_type='NONE')
    
    # ms2pip predictions
    # predict_ms2pip_hcd2021(LIBRARY_PATH, OUTPUT_DIR)
    # predict_ms2pip_immuno_hcd(LIBRARY_PATH, OUTPUT_DIR)
    
    print("\nProcessing complete!")
