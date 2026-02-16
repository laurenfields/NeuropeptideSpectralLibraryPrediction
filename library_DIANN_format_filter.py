#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
DIA-NN Library Formatter and Filter
====================================

This script combines intensity and retention time predictions from different models
and formats them into DIA-NN compatible spectral libraries. It also applies filtering
to optimize library size and quality.

The script:
1. Merges intensity predictions with RT predictions
2. Converts fragment annotations to DIA-NN format
3. Calculates precursor m/z values
4. Filters peptides by length
5. Selects top N most intense fragments per peptide

Author: lafields2
Created: Fri Apr 25 13:02:27 2025
"""

import pandas as pd
import pyopenms as oms

# ============================================================================
# Configuration
# ============================================================================

# Input paths
PREDICTED_LIBRARY_PATH = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_building\Round2_FullDB\Prosit_2024_intensity_PTMs_gl\predicted_library.csv"
PREDICTED_RT_DIR = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_building\RT_eval"

# Output directory
OUTPUT_DIRECTORY = r"D:\Manuscripts\2024_NN_SpectralLibrary_Benchmarking\Library_building\libraries_for_evaluation"

# Model names
INTENSITY_MODEL_NAME = 'Prosit_2024_intensity_PTMs_gl'

# List of RT models to combine with the intensity predictions
RT_MODEL_NAME_LIST = [
    'AlphaPeptDeep_rt_generic',
    'Chronologer_RT',
    'Deeplc_hela_hf',
    'Prosit_2019_irt',
    'Prosit_2024_irt_PTMs_gl'
]

# Filtering parameters
N_TOP_FRAGMENTS = 12  # Number of most intense fragments to keep per peptide
MAX_PEPTIDE_LENGTH = 45  # Maximum peptide length (in amino acids)


# ============================================================================
# Helper Functions
# ============================================================================

def calculate_precursor_mz(sequence, charge):
    """
    Calculate precursor m/z value using pyOpenMS.
    
    Parameters
    ----------
    sequence : str
        Peptide sequence with modifications in bracket notation
        (e.g., "PEPTIDE[UNIMOD:35]K")
    charge : int
        Precursor charge state
        
    Returns
    -------
    float
        Precursor m/z value
    """
    # Convert UNIMOD notation to OpenMS format
    # OpenMS uses parentheses instead of square brackets
    sequence_openms = sequence.replace('[', '(').replace(']', ')')
    
    # Parse sequence and calculate m/z
    seq = oms.AASequence.fromString(sequence_openms)
    mz = seq.getMZ(charge)
    
    return mz


def parse_fragment_annotation(annotation):
    """
    Parse fragment annotation into components.
    
    Fragment annotations follow the format: "fr[ion_type][number]^[charge]"
    Example: "frb3^1" = b-ion, position 3, charge 1+
    
    Parameters
    ----------
    annotation : str
        Fragment annotation string
        
    Returns
    -------
    tuple
        (fragment_type, fragment_number, fragment_charge, loss_type)
    """
    # Extract fragment type (b or y)
    fragment_type = annotation[2]  # Position 2 contains 'b' or 'y'
    
    # Extract fragment charge (last character before '^')
    fragment_charge = annotation[-2]
    
    # Extract fragment number (between type and charge)
    fragment_index_raw = annotation[3:]
    fragment_number = fragment_index_raw[:-3]
    
    # Loss type (currently only 'None' supported, but could be extended)
    loss_type = 'None'
    
    return fragment_type, fragment_number, fragment_charge, loss_type


def merge_intensity_and_rt(intensity_lib, rt_lib):
    """
    Merge intensity predictions with RT predictions.
    
    Both libraries must have matching peptide sequences, modifications,
    and charge states.
    
    Parameters
    ----------
    intensity_lib : pd.DataFrame
        Library with intensity predictions
    rt_lib : pd.DataFrame
        Library with RT/iRT predictions
        
    Returns
    -------
    pd.DataFrame
        Merged library with both intensity and RT information
    """
    merge_columns = [
        'peptide_sequences',
        'peptide',
        'n_mods',
        'precursor_charges',
        'collision_energies'
    ]
    
    merged_lib = pd.merge(
        intensity_lib,
        rt_lib,
        on=merge_columns,
        how='inner'
    )
    
    return merged_lib


def format_for_diann(predicted_lib):
    """
    Format predicted library into DIA-NN format.
    
    DIA-NN requires specific column names and data organization:
    - ModifiedPeptide: Peptide sequence with modifications
    - PrecursorMz: Precursor m/z value
    - iRT: Indexed retention time
    - FragmentMz, FragmentType, FragmentCharge, etc.: Fragment information
    
    Parameters
    ----------
    predicted_lib : pd.DataFrame
        Merged library with intensity and RT predictions
        
    Returns
    -------
    pd.DataFrame
        Library formatted for DIA-NN
    """
    # Calculate precursor m/z values
    mz_store = []
    for idx in range(len(predicted_lib)):
        sequence = predicted_lib['peptide_sequences'].iloc[idx]
        charge = int(predicted_lib['precursor_charges'].iloc[idx])
        mz = calculate_precursor_mz(sequence, charge)
        mz_store.append(mz)
    
    # Parse fragment annotations
    loss_type = []
    frag_charge_store = []
    frag_type_store = []
    frag_index_store = []
    
    for annotation in predicted_lib['annotation']:
        frag_type, frag_number, frag_charge, loss = parse_fragment_annotation(annotation)
        frag_type_store.append(frag_type)
        frag_index_store.append(frag_number)
        frag_charge_store.append(frag_charge)
        loss_type.append(loss)
    
    # Create formatted library
    formatted_lib = pd.DataFrame()
    formatted_lib['ModifiedPeptide'] = predicted_lib['peptide_sequences']
    formatted_lib['LabeledPeptide'] = predicted_lib['peptide_sequences']
    formatted_lib['StrippedPeptide'] = predicted_lib['peptide']
    formatted_lib['PrecursorCharge'] = predicted_lib['precursor_charges']
    formatted_lib['PrecursorMz'] = mz_store
    
    # Handle RT/iRT column (different models use different names)
    try:
        formatted_lib['iRT'] = predicted_lib['irt']
    except KeyError:
        formatted_lib['iRT'] = predicted_lib['rt']
    
    formatted_lib['CollisionEnergy'] = predicted_lib['collision_energies']
    formatted_lib['ProteinIds'] = 'unknown'  # Protein IDs would be added separately
    formatted_lib['RelativeFragmentIntensity'] = predicted_lib['intensities']
    formatted_lib['FragmentMz'] = predicted_lib['mz']
    formatted_lib['FragmentNumber'] = frag_index_store
    formatted_lib['FragmentType'] = frag_type_store
    formatted_lib['FragmentCharge'] = frag_charge_store
    formatted_lib['FragmentLossType'] = loss_type
    
    return formatted_lib


def filter_and_select_fragments(formatted_lib, max_length, n_fragments):
    """
    Apply filtering and fragment selection to optimize library.
    
    Steps:
    1. Filter peptides by maximum length
    2. Select top N most intense fragments per peptide
    
    Parameters
    ----------
    formatted_lib : pd.DataFrame
        Formatted DIA-NN library
    max_length : int
        Maximum peptide length (amino acids)
    n_fragments : int
        Number of top fragments to keep per peptide
        
    Returns
    -------
    pd.DataFrame
        Filtered and optimized library
    """
    # Filter by peptide length
    filtered_lib = formatted_lib[
        formatted_lib['StrippedPeptide'].str.len() <= max_length
    ].reset_index(drop=True)
    
    print(f"After length filter (≤{max_length} AA): {len(filtered_lib['ModifiedPeptide'].unique())} unique peptides")
    
    # Select top N most intense fragments per peptide
    filtered_lib = filtered_lib.groupby("ModifiedPeptide", group_keys=False).apply(
        lambda x: x.nlargest(n_fragments, 'RelativeFragmentIntensity')
    )
    
    filtered_lib.reset_index(drop=True, inplace=True)
    
    print(f"After selecting top {n_fragments} fragments: {len(filtered_lib)} total entries")
    
    return filtered_lib


# ============================================================================
# Main Processing
# ============================================================================

def process_library_combination(intensity_path, rt_path, rt_model_name,
                                 intensity_model_name, output_dir,
                                 max_length=45, n_fragments=12):
    """
    Process a single intensity + RT model combination.
    
    Parameters
    ----------
    intensity_path : str
        Path to intensity predictions CSV
    rt_path : str
        Path to RT predictions CSV
    rt_model_name : str
        Name of RT model
    intensity_model_name : str
        Name of intensity model
    output_dir : str
        Output directory
    max_length : int
        Maximum peptide length filter
    n_fragments : int
        Number of fragments to keep per peptide
        
    Returns
    -------
    str
        Path to output file
    """
    print(f"\n{'='*60}")
    print(f"Processing: {intensity_model_name} + {rt_model_name}")
    print(f"{'='*60}")
    
    # Load predictions
    print("Loading predictions...")
    predicted_lib_no_rt = pd.read_csv(intensity_path)
    predicted_rt = pd.read_csv(rt_path)
    
    print(f"Intensity library: {len(predicted_lib_no_rt)} entries")
    print(f"RT library: {len(predicted_rt)} entries")
    
    # Merge intensity and RT predictions
    print("Merging predictions...")
    predicted_lib = merge_intensity_and_rt(predicted_lib_no_rt, predicted_rt)
    print(f"Merged library: {len(predicted_lib)} entries")
    
    # Format for DIA-NN
    print("Formatting for DIA-NN...")
    formatted_lib = format_for_diann(predicted_lib)
    
    # Apply filters and fragment selection
    print("Applying filters...")
    final_lib = filter_and_select_fragments(formatted_lib, max_length, n_fragments)
    
    # Save output
    output_filename = f'unrefined_{intensity_model_name}_intensity_{rt_model_name}_RT_library.csv'
    output_path = os.path.join(output_dir, output_filename)
    final_lib.to_csv(output_path, index=False)
    
    print(f"Saved to: {output_path}")
    print(f"Final library: {len(final_lib['ModifiedPeptide'].unique())} unique peptides")
    
    return output_path


# ============================================================================
# Main Execution
# ============================================================================

if __name__ == "__main__":
    """
    Main execution block.
    
    Generates DIA-NN formatted libraries for all combinations of
    intensity and RT predictions.
    """
    import os
    
    print("="*60)
    print("DIA-NN Library Formatter and Filter")
    print("="*60)
    print(f"Intensity model: {INTENSITY_MODEL_NAME}")
    print(f"RT models: {', '.join(RT_MODEL_NAME_LIST)}")
    print(f"Maximum peptide length: {MAX_PEPTIDE_LENGTH} AA")
    print(f"Fragments per peptide: {N_TOP_FRAGMENTS}")
    print("="*60)
    
    # Process each RT model combination
    for rt_model_name in RT_MODEL_NAME_LIST:
        predicted_rt_path = os.path.join(
            PREDICTED_RT_DIR,
            rt_model_name,
            'predicted_library.csv'
        )
        
        process_library_combination(
            intensity_path=PREDICTED_LIBRARY_PATH,
            rt_path=predicted_rt_path,
            rt_model_name=rt_model_name,
            intensity_model_name=INTENSITY_MODEL_NAME,
            output_dir=OUTPUT_DIRECTORY,
            max_length=MAX_PEPTIDE_LENGTH,
            n_fragments=N_TOP_FRAGMENTS
        )
    
    print("\n" + "="*60)
    print("All library combinations processed successfully!")
    print("="*60)
