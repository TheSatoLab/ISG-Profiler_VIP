# ISG-VIP: Viral Infection Predictor

Machine learning-based viral infection prediction using ISG (Interferon-Stimulated Gene) expression profiles.

## Requirements

- Python 3.12
- OpenMP (Required for LightGBM)

## Quick Start

### 1. Preparation

First, get the source code.

**Get source**

```bash
git clone https://github.com/TheSatoLab/ISG-Profiler_VIP.git

cd ISG-Profiler_VIP/isg-vip/
```

### 2. Installation

#### Step 2-1: Prepare System Dependencies (LightGBM / OpenMP)

This application requires LightGBM. Please install the necessary libraries for your operating system.

**For Mac OS (Apple Silicon / Intel)**
Homebrew is required.

```bash
brew install lightgbm
```

**For Linux (Debian / Ubuntu)**
If you encounter `OSError: libgomp.so.1: cannot open shared object file`, install the OpenMP library:

```bash
apt-get update && apt-get install -y libgomp1
```

#### Step 2-2: Install Python Dependencies

Choose **one** of the following methods to install the Python environment.

- **Option A (Recommended):** Use the existing conda environment from ISG-Profiler.
- **Option B:** Create a new virtual environment manually.

**Option A: Install via conda (Recommended)**

If you have already set up the `isg-profiler` environment following the [ISG-Profiler README](../isg-profiler/README.md#option-a-install-via-conda-recommended), you can reuse it.

```bash
# Activate the existing environment
conda activate isg-profiler

# Install ISG-VIP into this environment
pip install .
```

**Option B: Install manually (pip)**

Use this method if you cannot use conda. Choose one of the following methods:

**Method 1: Create virtual env (Recommended)**

```bash
# Check your current directory.
pwd # Output: {YOUR_INSTALLTION_PATH}/ISG-Profiler_VIP/isg-vip

# Create virtual env
python3 -m venv .venv
source .venv/bin/activate

## Update venv pip
pip install --upgrade pip
## Install requirements
pip install .
```

**Method 2: Install directly**

```bash
# Check your current directory.
pwd # Output: {YOUR_INSTALLTION_PATH}/ISG-Profiler_VIP/isg-vip

pip install .
```

### 3. Tutorial: Run with example data

Before processing your actual samples, let's verify the installation using the **example dataset**.

**1. Check current directory**

Ensure you are in the project root directory (`isg-vip`).

```bash
# Check current directory path
pwd
# Example output: .../ISG-Profiler_VIP/isg-vip
```

**2. Prepare example data**

> [!TIP]
> **Already downloaded?**
> If you have already downloaded `example_files_20260126.zip` in the ISG-Profiler tutorial, you can skip the download. Copy or move the extracted directory to the current location, or simply ensure it is accessible.

If you don't have the example data yet, download it from GitHub Releases:

```bash
# Download example data zip
curl -f -L -O https://github.com/TheSatoLab/ISG-Profiler_VIP/releases/download/example_files_20260126/example_files_20260126.zip

# Unzip (This will create 'example_files_20260126' directory)
unzip example_files_20260126.zip
```

**3. Run ISG-VIP**

Run the prediction using the example files. We explicitly specify the input file paths using command line arguments.

```bash
python3 -m isg_vip \
  --gene_count_file example_files_20260126/isg-vip/input/per_gene_count.tsv \
  --metadata example_files_20260126/isg-vip/input/sample_metadata.tsv \
  --output example_output
```

**4. Check the results**

If the command finishes successfully, check the `example_output` directory.

```bash
ls -F example_output/
# Output: Infection_Prediction_0.csv ... Infection_Prediction_Stacking_final.csv
```

> [!TIP]
> **Verification Complete!** Now you are ready to process your actual samples. Proceed to the next step.

### 4. Prepare sample files

Place your input files in the `input` directory. The following files are required:

- `per_gene_count.tsv`: ISG expression data.
- `sample_metadata.tsv`: Sample metadata.

#### `per_gene_count.tsv`

This file corresponds to the `per_gene_count.tsv` output from ISG-Profiler (generated **without** the `--per_species` option).
For details on how to generate this file, please refer to the [ISG-Profiler README](../isg-profiler/README.md#per_gene_counttsv-without---per-spices).

| Column Header          | Data Type | Description                  |
| :--------------------- | :-------- | :--------------------------- |
| **sample_id**          | String    | Same as ISG-Profiler output. |
| **hum_symbol**         | float     | Same as ISG-Profiler output. |
| **raw_count**          | String    | Same as ISG-Profiler output. |
| **type**               | String    | Same as ISG-Profiler output. |
| **normalized_count**   | String    | Not used in ISG-VIP.         |
| **standardized_count** | String    | Not used in ISG-VIP.         |

#### `sample_metadata.tsv`

This file is the same as the ISG-Profiler input `sample_metadata.tsv`.

| Column Header    | Data Type | Description           |
| :--------------- | :-------- | :-------------------- |
| **sample_id**    | String    | Same as ISG-Profiler. |
| **species_host** | String    | Same as ISG-Profiler. |
| **order_host**   | String    | Same as ISG-Profiler. |
| **clade_host**   | String    | Not used in ISG-VIP.  |

> [!IMPORTANT]
> Always verify that your `per_gene_count.tsv` and `sample_metadata.tsv` share the same **sample_id**s.

### 5. Execute

Run the prediction tool using `python3`.

```bash
# Runs ISG-VIP using default paths
python3 -m isg_vip

# Runs ISG-VIP with specific file paths
python3 -m isg_vip \
  --gene_count_file input/per_gene_count.tsv \
  --metadata input/sample_metadata.tsv \
  --output output
```

#### **isg_vip** Command Line Options

| Option              | Description                    | Default Value               |
| :------------------ | :----------------------------- | :-------------------------- |
| `-h, --help`        | Show help message.             | -                           |
| `--gene_count_file` | Path to `per_gene_count.tsv`.  | `input/per_gene_count.tsv`  |
| `--metadata`        | Path to `sample_metadata.tsv`. | `input/sample_metadata.tsv` |
| `--output`          | Output directory.              | `output`                    |

## Outputs

Results are written to the directory specified by `--output` (default: `output`) in **CSV** format.

| File Name                                          | Description                                                                          |
| :------------------------------------------------- | :----------------------------------------------------------------------------------- |
| `Infection_Prediction_{0-4}.csv`                   | Individual fold predictions generated by base models (LightGBM, LogisticRegression). |
| `Infection_Prediction_Stacking_{0-4}_external.csv` | Stacking ensemble results for each specific fold.                                    |
| `Infection_Prediction_Stacking_all.csv`            | All stacking ensemble results from each specific fold results.                       |
| `Infection_Prediction_Stacking_final.csv`          | Final consolidated predictions derived from a 5-fold majority vote.                  |

## Model Architecture

5-fold stacking ensemble:

- Base models: LightGBM + Logistic Regression
- Meta model: LightGBM (trained on base model outputs)
- Final prediction: Majority vote across folds

## Security & Model Integrity

This project uses `pickle` format for trained models. Since loading pickle files can execute arbitrary code, we provide SHA-256 checksums to ensure the files have not been tampered with.

If you want to verify the integrity of the model files, run the following command in the root directory:

```bash
# Verify that the model files match the original checksums
sha256sum -c checksums.sha256
```

## Troubleshooting

### Error: `No module named isg_vip`

If you see the following error:

> `/path/to/python3: `No module named isg_vip`

This means the Python package was not installed correctly in your current environment. This often happens if the `pip install .` step was skipped or run in the wrong directory.

**Solution:**
Navigate to the repository root (`isg-vip`), activate your virtual environment, and reinstall:

```bash
# Ensure you are in the correct directory
cd ISG-Profiler_VIP/isg-vip/

# Activate your environment
conda activate isg-profiler  # or 'source .venv/bin/activate'

# Re-install the package
pip install .
```

### Installation Conflicts / Pip Errors

If `pip install .` fails or behaves unexpectedly (e.g., conflicts with global user packages installed in your home directory):

**Solution:**
Export `PYTHONNOUSERSITE=1` to ignore local user packages before installing.

```bash
# Prevent conflicts with user-site packages
export PYTHONNOUSERSITE=1

# Install dependencies
pip install .
```

## Developer Guide

### Install dependencies

- Python 3.12.x is required.
- For update Python, test by example files and rewrite pyproject.toml

```bash
# Create virtual env
python3 -m venv .venv
source .venv/bin/activate

# Update pip
pip install --upgrade pip

# Install dependencies
python3 -m pip install -e ".[dev]"
```

### Dependencies version

When upgrading dependencies, please also pin the ISG-Profiler dependencies to the SAME VERSION.

### Notice for Model Updates

When updating the model, please ensure consistency in the scikit-learn version used during training.

Refer to the pyproject.toml in this repository and pin the scikit-learn dependency to the specific version used for training. The model files are located in the `model_dir` directory.

Additionally, you must update the checksums whenever any files in model_dir/ are added or modified. This allows users to verify that the .pkl files have not been tampered with. To update the checksums.sha256 file, run the following command from the project root:

```bash
find src/isg_vip/model_dir/ -type f \( -name "*.pkl" -o -name "*.npy" \) -exec sha256sum {} + | sort -k 2 > checksums.sha256
```

Always commit the updated `checksums.sha256` alongside the new model files.
