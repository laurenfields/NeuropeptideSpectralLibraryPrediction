#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Koina Retention Time Prediction
================================

This script generates retention time (RT) predictions using various RT prediction
models available through the Koina API. These predictions can be used to enhance
spectral library quality for DIA proteomics analysis.

Supported RT prediction models:
- AlphaPeptDeep_rt_generic
- Chronologer_RT
- Deeplc_hela_hf
- Prosit_2019_irt
- Prosit_2024_irt_PTMs_gl

Author: lafields2
Created: Tue Apr 22 15:42:02 2025
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
OUTPUT_DIR = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_building\RT_eval"
LIBRARY_PATH = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_building\Round2_FullDB\crustaceanDB_wMods.csv"

# Koina server configuration
KOINA_SERVER = "koina.wilhelmlab.org:443"


# ============================================================================
# Helper Functions
# ============================================================================

def prepare_library_for_model(library_df, model_name):
    """
    Prepare library dataframe based on model-specific requirements.
    
    Different RT models have different requirements regarding:
    - Supported modifications
    - Maximum peptide length
    - Minimum peptide length
    - Modification notation format
    
    Parameters
    ----------
    library_df : pd.DataFrame
        Input library dataframe
    model_name : str
        Name of the RT prediction model
        
    Returns
    -------
    pd.DataFrame
        Prepared library meeting model requirements
    """
    library = library_df.copy()
    
    # Standard modification reformatting (N-terminal modifications)
    # Some models require modifications after the amino acid
    replacements = {
        r'\[UNIMOD:28\]Q': r'Q[UNIMOD:28]',  # Pyroglutamate from Q
        r'\[UNIMOD:27\]E': r'E[UNIMOD:27]'   # Pyroglutamate from E
    }
    for old, new in replacements.items():
        library['peptide_sequences'] = library['peptide_sequences'].str.replace(
            old, new, regex=True)
    
    # Model-specific filtering
    if model_name == "AlphaPeptDeep_rt_generic":
        # Does not support acetylation at protein N-terminus
        modifications_to_remove = ['[UNIMOD:2]']
        pattern = '|'.join(map(re.escape, modifications_to_remove))
        library = library[~library['peptide_sequences'].str.contains(
            pattern, regex=True, na=False)].reset_index(drop=True)
    
    elif model_name == "Deeplc_hela_hf":
        # Length constraints: 4-60 amino acids
        library = library[library['peptide'].str.len() <= 60].reset_index(drop=True)
        library = library[library['peptide'].str.len() >= 4].reset_index(drop=True)
    
    elif model_name in ["Prosit_2019_irt", "Prosit_2024_irt_PTMs_gl"]:
        # Maximum length: 30 amino acids
        library = library[library['peptide'].str.len() <= 30].reset_index(drop=True)
        
        # Prosit 2019 does not support PTMs
        if model_name == "Prosit_2019_irt":
            modifications_to_remove = [
                '[UNIMOD:40]',  # Phosphorylation
                '[UNIMOD:28]',  # Pyroglutamate from Q
                '[UNIMOD:27]',  # Pyroglutamate from E
                '[UNIMOD:2]'    # Acetylation
            ]
            pattern = '|'.join(map(re.escape, modifications_to_remove))
            library = library[~library['peptide_sequences'].str.contains(
                pattern, regex=True, na=False)].reset_index(drop=True)
        
        # Prosit 2024 currently does not support these specific PTMs
        elif model_name == "Prosit_2024_irt_PTMs_gl":
            modifications_to_remove = [
                '[UNIMOD:40]',  # Phosphorylation
                '[UNIMOD:28]',  # Pyroglutamate from Q
                '[UNIMOD:27]',  # Pyroglutamate from E
                '[UNIMOD:2]'    # Acetylation
            ]
            pattern = '|'.join(map(re.escape, modifications_to_remove))
            library = library[~library['peptide_sequences'].str.contains(
                pattern, regex=True, na=False)].reset_index(drop=True)
    
    return library


def predict_and_save(library_df, model_name, output_dir, debug=False):
    """
    Generate RT predictions and save results.
    
    Parameters
    ----------
    library_df : pd.DataFrame
        Prepared library for prediction
    model_name : str
        Name of the Koina RT model
    output_dir : str
        Base output directory
    debug : bool, optional
        Enable debug mode for predictions (default: False)
        
    Returns
    -------
    str
        Path to saved predictions
    """
    print(f"\n{'='*60}")
    print(f"Predicting with {model_name}")
    print(f"Processing {len(library_df)} peptides...")
    print(f"{'='*60}")
    
    # Initialize model and generate predictions
    model = Koina(model_name, KOINA_SERVER)
    
    # Different models may require different prediction modes
    if model_name == "Deeplc_hela_hf":
        # DeepLC requires synchronous mode
        predictions = model.predict(library_df, debug=True, mode='sync')
    else:
        predictions = model.predict(library_df)
    
    # Create output subdirectory
    output_subdir = os.path.join(output_dir, model_name)
    if not os.path.exists(output_subdir):
        os.makedirs(output_subdir, exist_ok=True)
    
    # Save predictions
    output_path = os.path.join(output_subdir, 'predicted_library.csv')
    predictions.to_csv(output_path, index=False)
    
    print(f"Saved predictions to: {output_path}")
    print(f"Prediction columns: {predictions.columns.tolist()}")
    
    return output_path


# ============================================================================
# Model-Specific Prediction Functions
# ============================================================================

def predict_alphapeptdeep_rt(library_path, output_dir):
    """
    Generate RT predictions using AlphaPeptDeep_rt_generic model.
    
    This is a generic retention time model suitable for various LC conditions.
    
    Limitations:
    - Does not support [UNIMOD:2] (acetylation at protein N-terminus)
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    """
    model_name = "AlphaPeptDeep_rt_generic"
    library = pd.read_csv(library_path)
    library = prepare_library_for_model(library, model_name)
    predict_and_save(library, model_name, output_dir)


def predict_chronologer_rt(library_path, output_dir):
    """
    Generate RT predictions using Chronologer_RT model.
    
    Chronologer is a deep learning model for retention time prediction
    with broad applicability across different chromatographic conditions.
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    """
    model_name = "Chronologer_RT"
    library = pd.read_csv(library_path)
    library = prepare_library_for_model(library, model_name)
    predict_and_save(library, model_name, output_dir)


def predict_deeplc_hela(library_path, output_dir):
    """
    Generate RT predictions using Deeplc_hela_hf model.
    
    This model is trained on HeLa cell lysate data and is optimized for
    similar biological matrices. Uses high-fidelity mode.
    
    Limitations:
    - Peptide length must be between 4-60 amino acids
    - Requires synchronous prediction mode
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    """
    model_name = "Deeplc_hela_hf"
    library = pd.read_csv(library_path)
    library = prepare_library_for_model(library, model_name)
    predict_and_save(library, model_name, output_dir, debug=True)


def predict_prosit_2019_irt(library_path, output_dir):
    """
    Generate iRT predictions using Prosit_2019_irt model.
    
    Predicts indexed retention time (iRT) values, which are normalized
    retention times relative to standard peptides.
    
    Limitations:
    - Maximum peptide length: 30 amino acids
    - Does not support post-translational modifications
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    """
    model_name = "Prosit_2019_irt"
    library = pd.read_csv(library_path)
    library = prepare_library_for_model(library, model_name)
    predict_and_save(library, model_name, output_dir)


def predict_prosit_2024_irt_ptms(library_path, output_dir):
    """
    Generate iRT predictions using Prosit_2024_irt_PTMs_gl model.
    
    Enhanced version with support for post-translational modifications.
    Predicts indexed retention time (iRT) values.
    
    Limitations:
    - Maximum peptide length: 30 amino acids
    - Limited PTM support (excludes some modifications)
    
    Parameters
    ----------
    library_path : str
        Path to input library CSV
    output_dir : str
        Output directory for predictions
    """
    model_name = "Prosit_2024_irt_PTMs_gl"
    library = pd.read_csv(library_path)
    library = prepare_library_for_model(library, model_name)
    predict_and_save(library, model_name, output_dir)


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    """
    Main execution block.
    
    Uncomment the desired prediction function(s) to run.
    Each function generates RT predictions using a different model.
    """
    
    print("="*60)
    print("Koina Retention Time Prediction")
    print("="*60)
    print(f"Input library: {LIBRARY_PATH}")
    print(f"Output directory: {OUTPUT_DIR}")
    print("="*60)
    
    # Run RT predictions for each model
    # Uncomment the models you want to use
    
    # predict_alphapeptdeep_rt(LIBRARY_PATH, OUTPUT_DIR)
    # predict_chronologer_rt(LIBRARY_PATH, OUTPUT_DIR)
    # predict_deeplc_hela(LIBRARY_PATH, OUTPUT_DIR)
    # predict_prosit_2019_irt(LIBRARY_PATH, OUTPUT_DIR)
    # predict_prosit_2024_irt_ptms(LIBRARY_PATH, OUTPUT_DIR)
    
    print("\nRT prediction complete!")
