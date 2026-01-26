# ISG-Profiler

Cross-species ISG profiling for hidden virus detection in animal RNA-Seq data.

## Requirements

- bash >= 3.2
- Python 3.12
- salmon 1.10.x (`https://github.com/COMBINE-lab/salmon`)

## Quick Start

### 1. Preparation (Common Steps)

First, get the source code and reference files. **This step is required for all users.**

**Get source**

```bash
git clone https://github.com/TheSatoLab/ISG-Profiler_VIP.git

cd ISG-Profiler_VIP/isg-profiler/
```

**Download salmon index**
Download salmon index and place it into `reference` directory.

```bash
# Download zip file
curl -f -L -O https://github.com/TheSatoLab/ISG-Profiler_VIP/releases/download/Isoform_241003_salmon/Isoform_241003_salmon.zip

# Unzip downloaded file
unzip Isoform_241003_salmon.zip
# and move it into `reference` directory
mv Isoform_241003_salmon reference/
```

Confirm the reference directory structure:

```bash
$ ls reference/Isoform_241003_salmon
complete_ref_lens.bin  duplicate_clusters.tsv  pos.bin           refAccumLengths.bin  refseq.bin
ctable.bin             info.json               pre_indexing.log  ref_indexing.log     seq.bin
ctg_offsets.bin        mphf.bin                rank.bin          reflengths.bin       versionInfo.json
```

### 2. Installation (Choose Method A or B)

Choose **one** of the following procedures to install Salmon and dependencies.

- **Option A (Recommended):** Install via conda (Supports x86_64 / ARM including Mac OS Apple Silicon)
- **Option B:** Install manually (Supports Linux x86_64 only)

> [!IMPORTANT]
> **For Apple Silicon users:** Direct salmon installation from the GitHub repository (Option B) is not supported due to build constraints. Please use **Option A (conda)**.

#### Option A: Install via conda (Recommended)

This method handles both Salmon and Python environment setup.

```bash
# Check your current directory.
pwd # Output: {YOUR_INSTALLTION_PATH}/ISG-Profiler_VIP/isg-profiler

# Configure conda channels:
conda config --add channels conda-forge
conda config --add channels bioconda

# Create virtual env and activate
conda create -n isg-profiler salmon python=3.12
conda activate isg-profiler

# After activation, install python dependencies
pip install .
```

[Optional] Remove defaults channels if required

```bash
# Please check conda / miniconda EULA BEFORE use default channels.
conda config --remove channels defaults

# Check default chanels is removed
conda config --show-sources
# If 'defaults' still appears above, remove it by specifying the target file:
conda config --file {TARGET_FILE} --remove channels defaults
```

> [!TIP]
> **Setup Complete!** If you chose Option A, please skip Option B and proceed to **"3. Tutorial"**.

#### Option B: Install manually (Linux x86_64)

Use this method if you cannot use conda or prefer to manage binaries manually.

**Step B-1: Install Salmon from GitHub**

> [!NOTE]
> Skip this step if you already have Salmon installed.
> If you prefer to build Salmon from source, please refer to the [official documentation](https://salmon.readthedocs.io/en/latest/building.html#installation).

```bash
# Set the installation directory (Please update this with your own absolute path)
YOUR_SALMON_DIR_INSTALL_PATH=""
# Set the bin path (Please update this with your own bin path)
SALMON_VERSION=1.10.0
mkdir -p "${YOUR_SALMON_DIR_INSTALL_PATH}"

curl --fail -L "https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-${SALMON_VERSION}_linux_x86_64.tar.gz" -o salmon.tar.gz

# Extract the downloaded archive
tar -xzf salmon.tar.gz -C "${YOUR_SALMON_DIR_INSTALL_PATH}" --strip-components=1

# Check your salmon install directory.
ls -F "${YOUR_SALMON_DIR_INSTALL_PATH}"
## Output: bin/  lib/  sample_data.tgz

# Create a symbolic link in your absolute PATH
YOUR_SALMON_BIN_PATH=""
mkdir -p "${YOUR_SALMON_BIN_PATH}"
ln -s "${YOUR_SALMON_DIR_INSTALL_PATH}"/bin/salmon "${YOUR_SALMON_BIN_PATH}"

# Check your installed directory
ls "${YOUR_SALMON_BIN_PATH}"
## Output: salmon@

# Clean up the downloaded file
rm salmon.tar.gz

"${YOUR_SALMON_BIN_PATH}/salmon" --version # Output: salmon 1.10.0
```

> [!NOTE]
> You can specify the absolute path to Salmon via a command-line option. See "**isg_profiler.sh** Command Line Options" for details.

**Step B-2: Install Python dependencies**

Choose one of the following methods:

**Method 1: Create virtual env (Recommended)**

```bash
# Check your current directory.
pwd # Output: {YOUR_INSTALLTION_PATH}/ISG-Profiler_VIP/isg-profiler

# Create virtual env
python3 -m venv .venv
source .venv/bin/activate

## Update venv pip
pip install --upgrade pip
## then install requirements
pip install .
```

**Method 2: Install directly**

```bash
# Check your current directory.
pwd # Output: {YOUR_INSTALLTION_PATH}/ISG-Profiler_VIP/isg-profiler

pip install .
```

---

### 3. Tutorial: Run with example data

Before processing your actual samples, let's verify the installation using a small **example dataset**.

**1. Check current directory**

Ensure you are in the project root directory (`isg-profiler`) where the main script is located.

```bash
# Check current directory path
pwd
# Example output: .../ISG-Profiler_VIP/isg-profiler

# Confirm the execution script exists
ls -F isg_profiler.sh
# Output: isg_profiler.sh*
```

**2. Download example data**

Get the example dataset from GitHub Releases.

```bash
# Download example data zip
curl -f -L -O https://github.com/TheSatoLab/ISG-Profiler_VIP/releases/download/example_files_20260126/example_files_20260126.zip

# Unzip (This will create 'example_files_20260126' directory)
unzip example_files_20260126.zip
```

**3. Run ISG-Profiler**

Run the script specifying the example directories.

```bash
./isg_profiler.sh \
  --fastq_dir example_files_20260126/isg-profiler/input/fastq \
  --metadata example_files_20260126/isg-profiler/input/sample_metadata.tsv \
  --out_dir example_output
```

**4. Check the results**

If the command finishes successfully, check the `example_output` directory.

```bash
ls -F example_output/
# Should output: ISG_score.tsv  per_gene_count.tsv  ...
```

> [!TIP]
> **Verification Complete!** Now you are ready to process your actual samples. Proceed to the next step.

### 4. Prepare sample files

Place your input files into the `input` directory.

Example: In case of two samples, ERR12917750 (paired-end) and ERR2012446 (single-end), the structure should look like this:

```bash
input
├── fastq
│   ├── ERR12917750_1.cleaned.fastq
│   ├── ERR12917750_2.cleaned.fastq
│   └── ERR2012446.cleaned.fastq
└── sample_metadata.tsv
```

#### 1. `sample_metadata.tsv`

Create a TSV file with the following columns:

| Column Header    | Data Type | Description                                                                                                    |
| :--------------- | :-------- | :------------------------------------------------------------------------------------------------------------- |
| **sample_id**    | String    | NCBI SRA RUN ID                                                                                                |
| **species_host** | String    | The scientific binomial name of the host organism. Spaces are replaced by underscores (e.g., `Gallus_gallus`). |
| **order_host**   | String    | The taxonomic order to which the host species belongs (e.g., `Galliformes`).                                   |
| **clade_host**   | String    | The broader taxonomic clade or class of the host (e.g., `Aves`, `Mammalia`).                                   |

> Check `sample_metadata.example.tsv` for a reference.

#### `fastq` directory

Store your **uncompressed** FASTQ files here.

> [!NOTE]
> Compressed files (e.g., `.fastq.gz`) are not supported. Please use uncompressed `.fastq` files.

| type            | file name                           | description                                                                                     |
| :-------------- | :---------------------------------- | :---------------------------------------------------------------------------------------------- |
| single-end read | `{NCBI_SRA_RUN_ID}.cleaned.fastq`   | RNA-Seq cleaned data (single-end read). Passed to `salmon quant` option `-r` [`--unmatedReads`] |
| paired-end read | `{NCBI_SRA_RUN_ID}_1.cleaned.fastq` | Read 1 FASTQ file, require adapter trimming. Passed to `salmon quant` option `-1` [`--mates1`]  |
| paired-end read | `{NCBI_SRA_RUN_ID}_2.cleaned.fastq` | Read 2 FASTQ file, require adapter trimming. Passed to `salmon quant` option `-2` [`--mates2`]  |

### Run script

```bash
# Run with default settings (threads=4, input=./input/fastq, output=./output)
./isg_profiler.sh

# To overwrite defaults, append arguments (see Command Line Options below):
./isg_profiler.sh \
  --thread 8 \
  --out_dir ./custom_output \
  --salmon_bin /my_salmon_path/bin/salmon
```

#### **isg_profiler.sh** Command Line Options

| Option                | Required | Default                     | Description                                    |
| :-------------------- | :------: | :-------------------------- | :--------------------------------------------- |
| `--thread <int>`      |    N     | `4`                         | Number of threads.                             |
| `--fastq_dir <path>`  |    N     | `input/fastq`               | Directory containing fastq files.              |
| `--out_dir <path>`    |    N     | `output`                    | Output directory for salmon.                   |
| `--ref_dir <path>`    |    N     | `reference`                 | Reference directory.                           |
| `--metadata <path>`   |    N     | `input/sample_metadata.tsv` | Sample metadata file.                          |
| `--salmon_bin <path>` |    N     | `salmon`                    | Path to salmon executable.                     |
| `--per_species`       |    N     | -                           | If set, group counts by hum_symbol and tax_id. |
| `--help`              |    N     | -                           | Show the help message and exit.                |

#### **quant_normalizer** Command Line Options

| Option              | Required | Description                                                                                                                                      |
| :------------------ | :------: | :----------------------------------------------------------------------------------------------------------------------------------------------- |
| `-h, --help`        |          | Show the help message and exit.                                                                                                                  |
| `--reference_dir`   |   Yes    | Directory containing reference files (e.g., gene2refseq list, Aves/Mars removal lists, mean & SD tables, and Amniota398_sp_id.list).             |
| `--sample_metadata` |   Yes    | Path to the sample metadata table (TSV). Required columns: `sample_id`, `species_host`, `order_host`, `clade_host`.                              |
| `--sf_dir`          |   Yes    | Directory containing Salmon `quant.sf` files. Files must be named as `<sample_id>_quant.sf`.                                                     |
| `--out_dir`         |   Yes    | Output directory path (e.g., `isg_profiler_out/result`).                                                                                         |
| `--per_species`     |          | If set, groups counts by both `hum_symbol` and `tax_id` and attaches species information. Note: **ISG_score will NOT be computed** in this mode. |
| `--log_level`       |          | Set the logging level. Choose from `info`, `debug`, or `warning`.                                                                                |

> [!IMPORTANT]
> Always verify that your `--sf_dir` and `--sample_metadata` share the same sample IDs to ensure mappings.

### Outputs

Salmon files:

- `{NCBI_SRA_RUN_ID}_quant.sf`: salmon quant result file

> [!NOTE]
> About salmon options, see [Salmon Documentation](https://salmon.readthedocs.io/en/latest/salmon.html).

Results written to `--out_dir` option (default: `output` directory):
default mode

- `ISG_score.tsv`
- `per_gene_count.tsv`

`--per-spices` mode

- `per_gene_per_spices_count.tsv`

> [!NOTE]
> Only default mode result can be used for ISG-VIP. See ISG-VIP documentation for detail.

#### default mode

##### `ISG_score.tsv`

| Column Header    | Data Type | Description                                        |
| :--------------- | :-------- | :------------------------------------------------- |
| **sample_id**    | String    | NCBI SRA RUN ID.                                   |
| **ISG score**    | Float     | Mean of total amount of normalized ISG expression. |
| **species_host** | String    | Added from `sample_metadata.tsv`.                  |
| **order_host**   | String    | Added from `sample_metadata.tsv`.                  |
| **clade_host**   | String    | Added from `sample_metadata.tsv`.                  |

##### `per_gene_count.tsv` (without `--per-spices`)

| Column Header          | Data Type | Description                         |
| :--------------------- | :-------- | :---------------------------------- |
| **sample_id**          | String    | NCBI SRA RUN ID.                    |
| **hum_symbol**         | float     | gene name symbol.                   |
| **raw_count**          | String    | Salmon raw count result.            |
| **type**               | String    | Label to check ISG or control gene. |
| **normalized_count**   | String    | Normalized count.                   |
| **standardized_count** | String    | standardized count.                 |

#### `--per-spices` mode

##### `per_gene_per_spices_count.tsv`

| Column Header          | Data Type | Description                                   |
| :--------------------- | :-------- | :-------------------------------------------- |
| **sample_id**          | String    | Same as default mode.                         |
| **hum_symbol**         | float     | Same as default mode.                         |
| **raw_count**          | String    | Same as default mode.                         |
| **type**               | String    | Same as default mode.                         |
| **normalized_count**   | String    | Same as default mode, but grouped by spicies. |
| **standardized_count** | String    | Same as default mode, but grouped by spicies. |
| **species**            | String    | Species genes to be mapped.                   |

### Execute normalizer directly:

You can excute normalize tool directly when you already have salmon quant.sf files:

```bash
python -m quant_normalizer \
  --sample_metadata input/sample_metadata.tsv \
  --sf_dir salmon_res \
  --out_dir isg_profiler_out/test_out
```

## Troubleshooting

### Salmon Execution Issues

If you encounter "command not found" for `salmon` or the script fails to run the quantification step:

1.  **Check Installation:**
    Ensure Salmon is installed and reachable:

    ```bash
    salmon --version
    ```

2.  **Specify Binary Path (Non-conda installs):**
    If you installed Salmon manually (e.g., from GitHub source) and it's not in your system `$PATH`, you must provide the path explicitly using the `--salmon_bin` option:

    ```bash
    ./isg_profiler.sh --salmon_bin /path/to/your/salmon
    ```

3.  **Recommendation:**
    We strongly recommend using **conda** for installation. Conda automatically manages the `$PATH` and dependencies, preventing these issues.

### Error: `No module named quant_normalizer`

If you see the following error:

> `/path/to/python3: No module named quant_normalizer`

This means the Python package was not installed correctly in your current environment. This often happens if the `pip install .` step was skipped or run in the wrong directory.

**Solution:**
Navigate to the repository root (`isg-profiler`), activate your virtual environment, and reinstall:

```bash
# Ensure you are in the correct directory
cd ISG-Profiler_VIP/isg-profiler/

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
python -m venv .venv
source .venv/bin/activate

# Update pip
pip install --upgrade pip

# Install dependencies
python -m pip install -e ".[dev]"
```

### Dependencis version

Pin the `scikit-learn` and `numpy` versions to match those used in ISG-VIP.

### Update Index file

To update the salmon index, please follow the steps below:

1. Create the salmon index.
1. Compress the index into a zip file named `Isoform_{YYMMDD}_salmon.zip`.  
   Example: For January 02, 2026, use `Isoform_260102_salmon.zip`.
1. Create a new release on the GitHub `main` branch. Set both the release title and the tag name to `Isoform_{YYMMDD}_salmon`.
1. Update the download URL and salmon index path in the Installation section of the README.
1. Update `SALMON_INDEX_DIR_NAME` in `isg-profiler/isg_profiler.sh` to `Isoform_YYMMDD_salmon`.
