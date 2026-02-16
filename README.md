# Spectral Library Generation and Evaluation Toolkit

A comprehensive Python toolkit for generating, formatting, evaluating, and comparing spectral libraries for DIA (Data-Independent Acquisition) proteomics using neural network-based prediction models.

## Overview

This toolkit provides end-to-end functionality for:
- **Library Generation**: Create predicted spectral libraries using multiple prediction models via the Koina API
- **RT Prediction**: Generate retention time predictions using various RT models
- **Library Formatting**: Format libraries for DIA-NN compatibility with filtering and optimization
- **Library Evaluation**: Automated DIA-NN analysis to assess library performance
- **Library Comparison**: Compare predicted vs experimental libraries using Pearson correlation
- **Intensity Normalization**: Normalize fragment intensities for consistent analysis

## Table of Contents

- [Installation](#installation)
- [Dependencies](#dependencies)
- [Scripts Overview](#scripts-overview)
- [Workflow](#workflow)
- [Usage Examples](#usage-examples)
- [Configuration](#configuration)
- [Output Files](#output-files)
- [Citation](#citation)
- [License](#license)

## Installation

### Prerequisites

- Python 3.8 or higher
- DIA-NN 1.8.1 or higher (for library evaluation)

### Python Dependencies

Install required Python packages:

```bash
pip install pandas numpy scipy matplotlib seaborn pyopenms koinapy more-itertools
```

### Package Versions (Tested)

```
pandas>=1.5.0
numpy>=1.23.0
scipy>=1.9.0
matplotlib>=3.6.0
seaborn>=0.12.0
pyopenms>=2.8.0
koinapy>=0.1.0
more-itertools>=9.0.0
```

## Dependencies

### Core Libraries

- **pandas**: Data manipulation and analysis
- **numpy**: Numerical computing
- **scipy**: Scientific computing (Pearson correlation)
- **matplotlib/seaborn**: Visualization and plotting
- **pyopenms**: Mass spectrometry data processing and m/z calculations
- **koinapy**: Interface to Koina prediction server
- **more-itertools**: Advanced iteration tools (chunking)

### External Software

- **DIA-NN**: For spectral library searching and evaluation
  - Download from: https://github.com/vdemichev/DiaNN

## Scripts Overview

### 1. `koina_lib_building.py`

Generates predicted spectral libraries using various intensity prediction models available through the Koina API.

**Supported Models:**
- AlphaPeptDeep_ms2_generic (QE, LUMOS)
- Prosit variants (2019, 2020, 2024, 2025)
- UniSpec (QE, QEHFX, LUMOS, NONE)
- ms2pip variants (HCD2021, Immuno_HCD)

**Key Features:**
- Model-specific filtering (peptide length, modifications, charge states)
- Batch processing for large datasets
- Automatic modification reformatting
- Support for multiple instrument types

**Output:** Predicted spectral libraries with fragment m/z and intensity predictions

---

### 2. `koina_RT_retrieve.py`

Generates retention time (RT) or indexed retention time (iRT) predictions.

**Supported RT Models:**
- AlphaPeptDeep_rt_generic
- Chronologer_RT
- Deeplc_hela_hf
- Prosit_2019_irt
- Prosit_2024_irt_PTMs_gl

**Key Features:**
- Model-specific constraints handling
- Modification notation conversion
- Support for both RT and iRT predictions

**Output:** RT/iRT predictions for peptide-charge combinations

---

### 3. `library_DIANN_format_filter.py`

Combines intensity and RT predictions, formats them for DIA-NN, and applies optimization filters.

**Key Features:**
- Merges intensity and RT predictions from different models
- Calculates precursor m/z values using pyOpenMS
- Converts fragment annotations to DIA-NN format
- Filters by peptide length
- Selects top N most intense fragments per peptide

**Output:** DIA-NN compatible spectral libraries (CSV format)

---

### 4. `diann_auto.py`

Automates library evaluation using DIA-NN for batch processing of multiple libraries.

**Key Features:**
- Processes multiple libraries automatically
- Consistent DIA-NN parameters across all libraries
- Generates comparative summary statistics
- Supports custom modification specifications

**Output:** 
- Individual DIA-NN reports for each library
- Summary CSV comparing all libraries
- Generated libraries from DIA-NN refinement

---

### 5. `compare2libs.py`

Compares two spectral libraries by calculating Pearson correlation coefficients for matched peptides.

**Key Features:**
- Peak matching within specified m/z tolerance
- Pearson correlation calculation for matched fragments
- Distribution visualization
- Summary statistics

**Output:**
- CSV file with per-peptide correlation coefficients
- Histogram plots of PCC distribution
- Summary statistics

---

### 6. `library_normalize.py`

Normalizes fragment ion intensities in spectral libraries on a per-peptide basis.

**Key Features:**
- Scales intensities relative to maximum per peptide
- Validation of normalization results
- Summary statistics reporting

**Output:** Normalized spectral library with intensities in [0, 1] range

## Workflow

### Complete Pipeline

```
1. Peptide Database
        ↓
2. Generate Intensity Predictions (koina_lib_building.py)
        ↓
3. Generate RT Predictions (koina_RT_retrieve.py)
        ↓
4. Format and Filter Libraries (library_DIANN_format_filter.py)
        ↓
5. Evaluate with DIA-NN (diann_auto.py)
        ↓
6. Compare Libraries (compare2libs.py)
```

### Alternative Workflows

**Experimental Library Processing:**
```
Experimental Library → Normalize (library_normalize.py) → Compare with Predicted
```

**Quick Comparison:**
```
Predicted Library → Compare with Experimental (compare2libs.py)
```

## Usage Examples

### Example 1: Generate Predicted Libraries

```python
# 1. Generate intensity predictions
python koina_lib_building.py

# Edit script to uncomment desired models:
# predict_prosit_2024_ptms(LIBRARY_PATH, OUTPUT_DIR)
# predict_unispec(LIBRARY_PATH, OUTPUT_DIR, instrument_type='QE')
```

### Example 2: Add RT Predictions

```python
# 2. Generate RT predictions
python koina_RT_retrieve.py

# Edit script to uncomment desired models:
# predict_prosit_2024_irt_ptms(LIBRARY_PATH, OUTPUT_DIR)
# predict_alphapeptdeep_rt(LIBRARY_PATH, OUTPUT_DIR)
```

### Example 3: Format for DIA-NN

```python
# 3. Combine predictions and format
python library_DIANN_format_filter.py

# Configure in script:
# INTENSITY_MODEL_NAME = 'Prosit_2024_intensity_PTMs_gl'
# RT_MODEL_NAME_LIST = ['Prosit_2024_irt_PTMs_gl', 'AlphaPeptDeep_rt_generic']
# N_TOP_FRAGMENTS = 12  # Number of fragments to keep
# MAX_PEPTIDE_LENGTH = 45  # Maximum peptide length
```

### Example 4: Evaluate Libraries

```python
# 4. Run DIA-NN evaluation
python diann_auto.py

# Configure in script:
# LIB_DIRECTORY = "path/to/libraries"
# SPECTRA_FILES = ["path/to/data.mzML"]
# OUT_FOLDER = "path/to/results"
```

### Example 5: Compare Libraries

```python
# 5. Compare predicted vs experimental
python compare2libs.py

# Configure in script:
# LIB1_PATH = "path/to/predicted_library.csv"
# LIB2_PATH = "path/to/experimental_library.csv"
# TOLERANCE = 0.1  # m/z tolerance
```

### Example 6: Normalize Experimental Library

```python
# Normalize experimental library
python library_normalize.py

# Configure in script:
# NON_NORM_PATH = "path/to/raw_library.csv"
# OUTPUT_PATH = "path/to/normalized_library.csv"
```

## Configuration

### File Paths

All scripts use absolute file paths for I/O. Update the following variables in each script:

**koina_lib_building.py:**
```python
OUTPUT_DIR = r"path/to/output"
LIBRARY_PATH = r"path/to/input_peptides.csv"
```

**koina_RT_retrieve.py:**
```python
OUTPUT_DIR = r"path/to/RT_predictions"
LIBRARY_PATH = r"path/to/input_peptides.csv"
```

**library_DIANN_format_filter.py:**
```python
PREDICTED_LIBRARY_PATH = r"path/to/intensity_predictions.csv"
PREDICTED_RT_DIR = r"path/to/RT_predictions"
OUTPUT_DIRECTORY = r"path/to/formatted_libraries"
```

**diann_auto.py:**
```python
LIB_DIRECTORY = r"path/to/libraries"
DIANN_PATH = r"path/to/DiaNN.exe"
OUT_FOLDER = r"path/to/results"
SPECTRA_FILES = [r"path/to/data.mzML"]
```

### Input File Formats

**Peptide Library (CSV):**
Required columns for Koina predictions:
- `peptide_sequences`: Modified peptide sequence (e.g., "PEPTIDE[UNIMOD:35]K")
- `peptide`: Unmodified peptide sequence
- `precursor_charges`: Charge state (integer)
- `collision_energies`: Collision energy (float)
- `n_mods`: Number of modifications (integer)

**Experimental Library (CSV):**
For comparison and normalization:
- `peptide`: Peptide sequence
- `peptide_sequences`: Modified peptide sequence
- `precursor_charges`: Charge state
- `mz`: Fragment m/z values
- `intensities`: Fragment intensities

## Output Files

### Koina Predictions

**Directory structure:**
```
output_dir/
├── AlphaPeptDeep_ms2_generic_QE/
│   └── predicted_library.csv
├── Prosit_2024_intensity_PTMs_gl/
│   └── predicted_library.csv
└── ...
```

**File contents:**
- Predicted fragment m/z values
- Predicted fragment intensities
- Fragment annotations
- Original peptide information

### RT Predictions

**Directory structure:**
```
RT_dir/
├── AlphaPeptDeep_rt_generic/
│   └── predicted_library.csv
├── Prosit_2024_irt_PTMs_gl/
│   └── predicted_library.csv
└── ...
```

**File contents:**
- RT or iRT values
- Original peptide information

### DIA-NN Formatted Libraries

**Format (CSV):**
```
ModifiedPeptide,LabeledPeptide,StrippedPeptide,PrecursorCharge,PrecursorMz,iRT,CollisionEnergy,ProteinIds,RelativeFragmentIntensity,FragmentMz,FragmentNumber,FragmentType,FragmentCharge,FragmentLossType
```

### DIA-NN Results

**Directory structure per library:**
```
library_name/
├── report.tsv          # Main results
├── generated_library.tsv  # Refined library
├── matrix_*.tsv        # Quantification matrices
└── ...
```

**Summary report:**
- Library name
- Total identification count
- Sorted by performance

### Comparison Results

**CSV file:**
```
key,pcc
PEPTIDE_2,0.95
PEPTIDE_3,0.87
...
```

**Plot files:**
- `*_tol0.1.png`: Standard PCC distribution (0-1.2)
- `*_tol0.1_extended.png`: Extended PCC distribution (-1.2 to 1.2)

## Model Limitations

### Intensity Models

| Model | Max Length | Max Charge | Supported PTMs | Notes |
|-------|-----------|-----------|----------------|-------|
| AlphaPeptDeep_ms2_generic | 32 AA | - | No [UNIMOD:2] | Instrument-specific |
| Prosit_2019_intensity | 30 AA | 6 | None | No PTM support |
| Prosit_2020_intensity_HCD | 30 AA | 6 | None | HCD only |
| Prosit_2024_intensity_PTMs_gl | 30 AA | 6 | Limited | Excludes some PTMs |
| UniSpec | 40 AA | 5 (QE/QEHFX) | Limited | Instrument-specific |
| ms2pip_HCD2021 | 30 AA | - | Most | HCD fragmentation |

### RT Models

| Model | Length Range | Supported PTMs | Output Type |
|-------|-------------|----------------|-------------|
| AlphaPeptDeep_rt_generic | - | No [UNIMOD:2] | RT |
| Chronologer_RT | - | Most | RT |
| Deeplc_hela_hf | 4-60 AA | Most | RT |
| Prosit_2019_irt | ≤30 AA | None | iRT |
| Prosit_2024_irt_PTMs_gl | ≤30 AA | Limited | iRT |

## UNIMOD Modifications

Common modifications referenced in the scripts:

- `[UNIMOD:2]`: Acetylation (N-terminus)
- `[UNIMOD:27]`: Pyroglutamate from E
- `[UNIMOD:28]`: Pyroglutamate from Q
- `[UNIMOD:35]`: Oxidation (M)
- `[UNIMOD:40]`: Phosphorylation (STY)

## Troubleshooting

### Common Issues

**1. Koina Connection Errors**
```
Solution: Check internet connection and verify Koina server status
Server: koina.wilhelmlab.org:443
```

**2. pyOpenMS Import Errors**
```
Solution: Install pyOpenMS with conda:
conda install -c bioconda pyopenms
```

**3. DIA-NN Not Found**
```
Solution: Update DIANN_PATH variable with full path to DiaNN.exe
Verify: Run DiaNN.exe --version in terminal
```

**4. Memory Issues with Large Libraries**
```
Solution: Reduce CHUNK_SIZE in koina_lib_building.py
Or process libraries in smaller batches
```

**5. Modification Format Errors**
```
Solution: Ensure modifications use UNIMOD notation with square brackets
Example: PEPTIDE[UNIMOD:35]K (not PEPTIDE(ox)K)
```

## Performance Tips

1. **Use batch processing** for large datasets (adjust `CHUNK_SIZE`)
2. **Filter early** to reduce downstream processing time
3. **Use multiple threads** in DIA-NN (`THREADS` parameter)
4. **Pre-filter peptides** by length before Koina predictions
5. **Cache predictions** to avoid redundant API calls

## Best Practices

### Library Generation

1. **Model Selection**: Choose models appropriate for your instrument and sample type
2. **Modification Support**: Verify model supports your PTMs before prediction
3. **Length Filtering**: Apply length filters before prediction to save time
4. **Validation**: Always validate predictions against experimental data

### Library Evaluation

1. **Consistent Parameters**: Use same DIA-NN settings across comparisons
2. **Representative Data**: Test on data similar to your analysis samples
3. **Multiple Metrics**: Evaluate libraries on multiple datasets
4. **Document Settings**: Keep detailed records of all parameters used

### Comparison Analysis

1. **Appropriate Tolerance**: Choose m/z tolerance based on instrument resolution
2. **Normalization**: Normalize intensities before comparison
3. **Statistical Rigor**: Report both mean/median and distribution
4. **Visual Inspection**: Always plot distributions to check for outliers

## Citation

If you use this toolkit in your research, please cite:
(citation to be added)

Additionally, please cite the relevant prediction models:

- **Prosit**: Gessulat et al. (2019) Nature Methods
- **AlphaPeptDeep**: Zeng et al. (2022) Nature Communications
- **UniSpec**: Diaz-Gay et al. (2023) Cell Systems
- **DeepLC**: Bouwmeester et al. (2021) Nature Methods
- **Koina**: Wilhelm et al. (2021) bioRxiv
- **DIA-NN**: Demichev et al. (2020) Nature Methods


---

**Version:** 1.0.0  
**Last Updated:** February 2026  
**Author:** lafields2
