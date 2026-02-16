#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DIA-NN Automated Library Evaluation
====================================

This script automates the evaluation of spectral libraries using DIA-NN.
It processes multiple libraries against DIA mass spectrometry data and
generates summary statistics for comparison.

The workflow:
1. Identifies all libraries in the specified directory
2. Runs DIA-NN for each library with consistent parameters
3. Collects identification results
4. Generates a summary report comparing all libraries

Author: lafields2
Created: Sat Apr 26 14:35:22 2025
"""

import os
import subprocess
import pandas as pd

# ============================================================================
# Configuration
# ============================================================================

# Directory containing libraries to evaluate
LIB_DIRECTORY = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_building\libraries_for_evaluation"

# DIA-NN executable path
DIANN_PATH = r"C:\DIA-NN\1.8.1\DiaNN.exe"

# Output directory for results
OUT_FOLDER = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_evaluation\figure_analyses\SG3"

# DIA mass spectrometry data files
SPECTRA_FILES = [
    r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_evaluation\RawData\Refinement_data\480\SG_DIA_newmeth_TR3.mzML"
]

# Alternative dataset (uncomment to use)
# SPECTRA_FILES = [
#     r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_evaluation\RawData\Refinement_data\480\PO_DIA_newmeth_TR3.mzML"
# ]

# Modification format for DIA-NN
UNIMOD_TONE = 'UniMod'

# DIA-NN search parameters
THREADS = 4
QVALUE_THRESHOLD = 0.05


# ============================================================================
# Helper Functions
# ============================================================================

def build_spectra_command(spectra_file_list):
    """
    Build the spectra file command string for DIA-NN.
    
    Parameters
    ----------
    spectra_file_list : list
        List of paths to mzML or raw files
        
    Returns
    -------
    str
        Formatted command string with all spectra files
    """
    spectra_command = ''
    for filepath in spectra_file_list:
        spectra_command += f' --f {filepath}'
    return spectra_command


def get_library_list(lib_directory):
    """
    Get list of all CSV libraries in the directory.
    
    Parameters
    ----------
    lib_directory : str
        Directory containing library CSV files
        
    Returns
    -------
    list
        List of full paths to library files
    """
    lib_list = []
    
    for filename in os.listdir(lib_directory):
        if filename.endswith(".csv"):
            lib_path = os.path.join(lib_directory, filename)
            lib_list.append(lib_path)
    
    return lib_list


def build_diann_command(diann_path, spectra_command, library_path,
                        output_report, output_lib, threads, qvalue,
                        unimod_tone):
    """
    Build complete DIA-NN command with all parameters.
    
    Parameters
    ----------
    diann_path : str
        Path to DIA-NN executable
    spectra_command : str
        Pre-formatted spectra file command string
    library_path : str
        Path to spectral library
    output_report : str
        Path for output report.tsv
    output_lib : str
        Path for generated library output
    threads : int
        Number of threads to use
    qvalue : float
        Q-value threshold for filtering
    unimod_tone : str
        Modification notation format
        
    Returns
    -------
    str
        Complete DIA-NN command string
    """
    command = (
        f'{diann_path} {spectra_command} '
        f'--lib {library_path} '
        f'--threads {threads} '
        f'--verbose 1 '
        f'--out {output_report} '
        f'--qvalue {qvalue} '
        f'--matrices '
        f'--out-lib {output_lib} '
        f'--gen-spec-lib '
        f'--reanalyse '
        f'--relaxed-prot-inf '
        f'--id-profiling '
        f'--peak-center '
        f'--no-ifs-removal '
        f'--global-norm '
        # Enzyme specificity: exclude cleavage after these residues
        f'--cut !*A,!*G,!*I,!*L,!*P,!*V,!*F,!*W,!*Y,!*D,!*E,!*R,!*H,!*K,!*S,!*T,!*C,!*M,!*N,!*Q '
        # Variable modifications
        f'--var-mod {unimod_tone}:2,-0.984016,*c '  # Acetylation with deamidation
        f'--var-mod {unimod_tone}:40,79.956815,STY '  # Phosphorylation
        f'--var-mod {unimod_tone}:35,15.994915,M '  # Oxidation
        f'--var-mod {unimod_tone}:28,-17.026549,Q*n '  # Pyroglutamate from Q
        f'--var-mod {unimod_tone}:27,-18.010565,E*n '  # Pyroglutamate from E
        # Additional parameters
        f'--global-mass-cal '
        f'--min-pep-len 3 '
        f'--int-removal 0 '
        f'--il-eq '  # Treat I and L as equivalent
        f'--strip-unknown-mods '
        f'--no-swissprot '
        f'--peak-translation '
        f'--relaxed-prot-inf'
    )
    
    return command


def run_diann_for_library(library_path, lib_directory, diann_path,
                           spectra_command, out_folder, threads,
                           qvalue, unimod_tone):
    """
    Run DIA-NN analysis for a single library.
    
    Parameters
    ----------
    library_path : str
        Path to spectral library CSV
    lib_directory : str
        Base library directory (for extracting library name)
    diann_path : str
        Path to DIA-NN executable
    spectra_command : str
        Pre-formatted spectra command string
    out_folder : str
        Base output folder
    threads : int
        Number of threads
    qvalue : float
        Q-value threshold
    unimod_tone : str
        Modification format
        
    Returns
    -------
    str
        Path to output subdirectory
    """
    # Extract library name from path
    file_name = library_path.replace(lib_directory, '')
    file_name = file_name.replace('.csv', '')
    file_name = file_name.lstrip('\\').lstrip('/')
    
    # Create output subdirectory
    output_subdir = os.path.join(out_folder, file_name)
    if not os.path.exists(output_subdir):
        os.makedirs(output_subdir, exist_ok=True)
    
    # Define output files
    output_report = os.path.join(output_subdir, 'report.tsv')
    output_lib = os.path.join(output_subdir, 'generated_library.tsv')
    
    # Build and execute command
    print(f"\n{'='*60}")
    print(f"Processing library: {file_name}")
    print(f"{'='*60}")
    
    command = build_diann_command(
        diann_path, spectra_command, library_path,
        output_report, output_lib, threads, qvalue, unimod_tone
    )
    
    print("Running DIA-NN...")
    exit_code = os.system(command)
    
    if exit_code == 0:
        print(f"Successfully completed: {file_name}")
    else:
        print(f"Warning: DIA-NN exited with code {exit_code} for {file_name}")
    
    return output_subdir


def summarize_results(out_folder):
    """
    Generate summary report from all DIA-NN results.
    
    Parameters
    ----------
    out_folder : str
        Base output folder containing all result subdirectories
        
    Returns
    -------
    pd.DataFrame
        Summary dataframe with library names and ID counts
    """
    # Get list of result subdirectories
    subfolders = [f.name for f in os.scandir(out_folder) if f.is_dir()]
    
    lib_store = []
    total_ids_store = []
    
    print(f"\n{'='*60}")
    print("Summarizing Results")
    print(f"{'='*60}")
    
    for subfolder in subfolders:
        results_path = os.path.join(out_folder, subfolder, 'report.tsv')
        
        # Check if report exists
        if not os.path.exists(results_path):
            print(f"Warning: No report found for {subfolder}")
            continue
        
        # Load and count results
        results = pd.read_csv(results_path, sep='\t')
        
        # Remove duplicates based on modified sequence
        results_unique = results.drop_duplicates(subset='Modified.Sequence')
        
        lib_store.append(subfolder)
        total_ids_store.append(len(results_unique))
        
        print(f"{subfolder}: {len(results_unique)} unique identifications")
    
    # Create summary dataframe
    summary_df = pd.DataFrame({
        'Library': lib_store,
        'ID_count': total_ids_store
    })
    
    # Sort by ID count (descending)
    summary_df = summary_df.sort_values('ID_count', ascending=False).reset_index(drop=True)
    
    return summary_df


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    """
    Main execution block.
    
    Processes all libraries in the specified directory and generates
    a comparative summary.
    """
    
    print("="*60)
    print("DIA-NN Automated Library Evaluation")
    print("="*60)
    print(f"Library directory: {LIB_DIRECTORY}")
    print(f"Output directory: {OUT_FOLDER}")
    print(f"DIA-NN path: {DIANN_PATH}")
    print(f"Number of spectra files: {len(SPECTRA_FILES)}")
    print(f"Q-value threshold: {QVALUE_THRESHOLD}")
    print("="*60)
    
    # Build spectra command
    spectra_command = build_spectra_command(SPECTRA_FILES)
    
    # Get list of libraries
    lib_list = get_library_list(LIB_DIRECTORY)
    print(f"\nFound {len(lib_list)} libraries to process")
    
    # Create output folder if it doesn't exist
    if not os.path.exists(OUT_FOLDER):
        os.makedirs(OUT_FOLDER, exist_ok=True)
    
    # Process each library
    for library_path in lib_list:
        run_diann_for_library(
            library_path=library_path,
            lib_directory=LIB_DIRECTORY,
            diann_path=DIANN_PATH,
            spectra_command=spectra_command,
            out_folder=OUT_FOLDER,
            threads=THREADS,
            qvalue=QVALUE_THRESHOLD,
            unimod_tone=UNIMOD_TONE
        )
    
    # Generate summary report
    summary_df = summarize_results(OUT_FOLDER)
    
    # Save summary
    summary_output = os.path.join(OUT_FOLDER, 'summary_report.csv')
    summary_df.to_csv(summary_output, index=False)
    
    print(f"\n{'='*60}")
    print("Processing Complete!")
    print(f"{'='*60}")
    print(f"Summary saved to: {summary_output}")
    print(f"\nTop 5 libraries by identification count:")
    print(summary_df.head().to_string(index=False))
    print("="*60)
