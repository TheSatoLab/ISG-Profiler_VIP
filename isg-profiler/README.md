# ISG-Profiler

Cross-species ISG profiling for hidden virus detection in animal RNA-Seq data.

## Requirements

- bash >= 3.2
- Python 3.12
- salmon 1.10.x (`https://github.com/COMBINE-lab/salmon`)

## Quick Start

### Installation

Get source

```bash
git clone https://github.com/TheSatoLab/ISG-Profiler_VIP.git

cd isg-profiler
```

Download salmon index and place it into `reference` directory.

```bash
# Download zip file
curl -f -L -O https://github.com/TheSatoLab/ISG-Profiler_VIP/releases/download/Isoform_241003_salmon/Isoform_241003_salmon.zip

# Unzip downloaded file
unzip Isoform_241003_salmon.zip
# and move it into `reference` directory
mv Isoform_241003_salmon reference/
```

```bash
$ ls reference/Isoform_241003_salmon
complete_ref_lens.bin  duplicate_clusters.tsv  pos.bin           refAccumLengths.bin  refseq.bin
ctable.bin             info.json               pre_indexing.log  ref_indexing.log     seq.bin
ctg_offsets.bin        mphf.bin                rank.bin          reflengths.bin       versionInfo.json
```

To install salmon and dependency, choose one of the following procedures:

- Install via conda (for x86_64 / ARM users)
- Install salmon from github repository (for x86_64 users only)

#### Install via conda

```bash
# Configure conda channels:
conda config --add channels conda-forge
conda config --add channels bioconda
# crate virtual env and activate
conda create -n isg-profiler salmon python=3.12
conda activate isg-profiler

# After activation, install dependency
pip install .
```

```bash
# [Optional] Remove defaults channels if required.
# Please check conda / miniconda EULA BEFORE use default channels.
conda config --remove channels defaults

# Check default chanels is removed
conda config --show-sources
# If 'defaults' still appears above, remove it by specifying the target file:
conda config --file {TARGET_FILE} --remove channels defaults
```

#### Install salmon from github repository

Install salmon from github repository:

```bash
SALMON_VERSION=1.10.0

mkdir -p /usr/local/salmon
curl --fail -L "https://github.com/COMBINE-lab/salmon/releases/download/v${SALMON_VERSION}/salmon-${SALMON_VERSION}_linux_x86_64.tar.gz" -o salmon.tar.gz

tar -xzf salmon.tar.gz -C /usr/local/salmon/ --strip-components=1 && \
    ln -s /usr/local/salmon/bin/salmon /usr/local/bin/salmon && \
    rm salmon.tar.gz

salmon --version # Output: salmon 1.10.0
```

Install requirements:

A) [Recommend] Create virtual env and install requirements:

```bash
# Create virtual env
python -m venv .venv
source .venv/bin/activate

## Update venv pip
pip install --upgrade pip
## then install requirements
pip install .
```

B) Install requirements directly:

```bash
pip install .
```

### Prepare sample files

Put samples to `input` directory:

- `fastq` direcotory (check following)
- `sample_metadata.tsv`

#### `sample_metadata.tsv`

TSV format, with folloing columns and contents:

| Column Header    | Data Type | Description                                                                                                    |
| :--------------- | :-------- | :------------------------------------------------------------------------------------------------------------- |
| **sample_id**    | String    | NCBI SRA BioProject_ID.                                                                                        |
| **species_host** | String    | The scientific binomial name of the host organism. Spaces are replaced by underscores (e.g., `Gallus_gallus`). |
| **order_host**   | String    | The taxonomic order to which the host species belongs (e.g., `Galliformes`).                                   |
| **clade_host**   | String    | The broader taxonomic clade or class of the host (e.g., `Aves`, `Mammalia`).                                   |

Check `sample_metadata.example.tsv` example file.

#### `fastq` directory

| type            | file name                         | description                                                                                     |
| :-------------- | :-------------------------------- | :---------------------------------------------------------------------------------------------- |
| single-end read | `{BIOPROJECT_ID}.cleaned.fastq`   | RNA-Seq cleaned data (single-end read). Passed to `salmon quant` option `-r` [`--unmatedReads`] |
| paired-end read | `{BIOPROJECT_ID}_1.cleaned.fastq` | Read 1 FASTQ file, require adapter trimming. Passed to `salmon quant` option `-1` [`--mates1`]  |
| paired-end read | `{BIOPROJECT_ID}_2.cleaned.fastq` | Read 2 FASTQ file, require adapter trimming. Passed to `salmon quant` option `-2` [`--mates2`]  |

Example: In case of two sample, ERR12917750 (paired end) and ERR2012446 (single end)

```bash
.
├── fastq
│   ├── ERR12917750_1.cleaned.fastq
│   ├── ERR12917750_2.cleaned.fastq
│   └── ERR2012446.cleaned.fastq
└── sample_metadata.example.tsv
```

### Run script

```bash
./isg_profiler.sh
```

#### **isg_profiler.sh** Command Line Options

| Option               | Required | Description                                                |
| :------------------- | :------: | :--------------------------------------------------------- |
| `--thread <int>`     |          | Number of threads (default: 4).                            |
| `--fastq_dir <path>` |          | Directory containing fastq files (default: input/fastq).   |
| `--out_dir <path>`   |          | Output directory for salmon (default: output).             |
| `--ref_dir <path>`   |          | Reference directory (default: reference).                  |
| `--metadata <path>`  |          | Sample metadata file (default: input/sample_metadata.tsv). |
| `--per_species`      |          | (Optional) If set, group counts by hum_symbol and tax_id.  |
| `--help`             |          | Show the help message and exit.                            |

### Execute normalizer directly:

You can excute normalize tool directly when you already have salmon quant.sf files:

```bash
python -m quant_normalizer \
  --sample_metadata input/sample_metadata.tsv \
  --sf_dir salmon_res \
  --out_prefix isg_profiler_out/test_out
```

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

- `{BIOPROJECT_ID}_quant.sf`: salmon quant result file

> [!NOTE]
> About salmon options, see [Salmon Documentation](https://salmon.readthedocs.io/en/latest/salmon.html).

Results written to `output-prefix` option:

| option \ File name  | `ISG_score.tsv` | `per_gene_count.tsv` |
| :------------------ | :-------------- | :------------------- |
| default mode        | ✓               | ✓                    |
| `--per-spices` mode | x               | ✓                    |

> [!NOTE]
> Only default mode result can be used for ISG-VIP. See ISG-VIP documentation for detail.

#### default mode

##### `ISG_score.tsv`

| Column Header    | Data Type | Description                                        |
| :--------------- | :-------- | :------------------------------------------------- |
| **sample_id**    | String    | NCBI SRA BioProject_ID.                            |
| **ISG score**    | Float     | Mean of total amount of normalized ISG expression. |
| **species_host** | String    | Added from `sample_metadata.tsv`.                  |
| **order_host**   | String    | Added from `sample_metadata.tsv`.                  |
| **clade_host**   | String    | Added from `sample_metadata.tsv`.                  |

##### `per_gene_count.tsv` (without `--per`)

| Column Header          | Data Type | Description                         |
| :--------------------- | :-------- | :---------------------------------- |
| **sample_id**          | String    | NCBI SRA BioProject_ID.             |
| **hum_symbol**         | float     | gene name symbol.                   |
| **raw_count**          | String    | Salmon raw count result.            |
| **type**               | String    | Label to check ISG or control gene. |
| **normalized_count**   | String    | Normalized count.                   |
| **standardized_count** | String    | standardized count.                 |

#### `--per-spices` mode

##### `per_gene_count.tsv`

| Column Header          | Data Type | Description                                   |
| :--------------------- | :-------- | :-------------------------------------------- |
| **sample_id**          | String    | Same as default mode.                         |
| **hum_symbol**         | float     | Same as default mode.                         |
| **raw_count**          | String    | Same as default mode.                         |
| **type**               | String    | Same as default mode.                         |
| **normalized_count**   | String    | Same as default mode, but grouped by spicies. |
| **standardized_count** | String    | Same as default mode, but grouped by spicies. |
| **species**            | String    | Species genes to be mapped.                   |

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
