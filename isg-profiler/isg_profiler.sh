#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Luca Nishimura & Jumpei Ito

set -eu

# === DEFAULT CONFIG ===
THREAD=4
FASTQ_DIR="input/fastq"
OUTPUT_DIR="output"
REF_DIR="reference"
## NOTE: Salmon index directory name is hard-coded because index `info.json` contains this string
SALMON_INDEX_DIR_NAME="Isoform_241003_salmon"
SAMPLE_METADATA="input/sample_metadata.tsv"
PER_SPECIES_OPT="" # Default is empty (= disabled)
# ======================

# Help messages
usage() {
  echo "Usage: $0 [OPTIONS]"
  echo "Options:"
  echo "  --thread <int>          Number of threads (default: ${THREAD})"
  echo "  --fastq_dir <path>      Directory containing fastq files (default: ${FASTQ_DIR})"
  echo "  --out_dir <path>        Output directory for salmon (default: ${OUTPUT_DIR})"
  echo "  --ref_dir <path>        Reference directory (default: ${REF_DIR})"
  echo "  --metadata <path>       Sample metadata file (default: ${SAMPLE_METADATA})"
  echo "  --per_species           (Optional) If set, group counts by hum_symbol and tax_id"
  echo "  --help                  Show this help message"
  exit 0
}

# === Parse arguments ===
while [[ $# -gt 0 ]]; do
  key="$1"
  # Error if argument has no value, except help and boolean flags
  case "$key" in
  --thread | --fastq_dir | --out_dir | --ref_dir | --metadata)
    if [[ -z "${2:-}" ]] || [[ "${2:-}" == --* ]]; then
      echo "Error: Argument for $key is missing"
      usage
    fi

    case "$key" in
    --thread) THREAD="$2" ;;
    --fastq_dir) FASTQ_DIR="$2" ;;
    --out_dir) OUTPUT_DIR="$2" ;;
    --ref_dir) REF_DIR="$2" ;;
    --metadata) SAMPLE_METADATA="$2" ;;
    esac
    shift 2
    ;;
  # bool flag (without value)
  --per_species)
    PER_SPECIES_OPT="--per_species"
    shift 1
    ;;
  --help)
    usage
    ;;
  *)
    echo "Error: Unknown option $1"
    usage
    ;;
  esac
done

# ==== macOS/Bash 3.2 Compatible Path Normalization ====
get_abs_path() {
  local target="$1"
  local abs_path

  if [ -d "$target" ]; then
    abs_path=$(cd "$target" >/dev/null 2>&1 && pwd -P)
  elif [ -f "$target" ]; then
    local dir_name=$(dirname "$target")
    local file_name=$(basename "$target")
    abs_path="$(cd "$dir_name" >/dev/null 2>&1 && pwd)/${file_name}"
  else
    echo "Error: '${target}' does not exist, cannot resolve absolute path." >&2
    exit 1
  fi
  echo "$abs_path"
}

# Rewrite all path
if [[ ! -d "${FASTQ_DIR}" ]]; then
  echo "Error: Input directory '${FASTQ_DIR}' does not exist." >&2
  exit 1
fi
FASTQ_DIR=$(get_abs_path "${FASTQ_DIR}")

if [[ ! -f "${SAMPLE_METADATA}" ]]; then
  echo "Error: Metadata file '${SAMPLE_METADATA}' does not exist." >&2
  exit 1
fi
SAMPLE_METADATA=$(get_abs_path "${SAMPLE_METADATA}")

if [[ ! -d "${REF_DIR}" ]]; then
  echo "Error: Reference directory '${REF_DIR}' does not exist." >&2
  exit 1
fi
REF_DIR=$(get_abs_path "${REF_DIR}")

mkdir -p "${OUTPUT_DIR}"
OUTPUT_DIR=$(get_abs_path "${OUTPUT_DIR}")

# =======================================================

# Define script directory
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
BIN_DIR="${SCRIPT_DIR}/bin"

echo "====== Configuration ======"
echo "Threads:              ${THREAD}"
echo "Input FASTQ Dir:      ${FASTQ_DIR}"
echo "Input metadata file:  ${SAMPLE_METADATA}"
echo "Reference Dir:        ${REF_DIR}"
echo "Output Dir:           ${OUTPUT_DIR}"

# Check PER_SPECIES_OPT variable and echo Enabled / Disabled
if [[ -n "${PER_SPECIES_OPT}" ]]; then
  echo "Per Species Mode:     Enabled"
else
  echo "Per Species Mode:     Disabled"
fi
echo "====================="

# WARNING
SALMON_INDEX_PATH="${REF_DIR}/${SALMON_INDEX_DIR_NAME}"

"${BIN_DIR}/run_salmon.sh" \
  --thread "${THREAD}" \
  --fastq_dir "${FASTQ_DIR}" \
  --out_dir "${OUTPUT_DIR}" \
  --salmon_index "${SALMON_INDEX_PATH}"

# isg_profiler output prefix setup
# 出力先ディレクトリ配下にプレフィックスを作成するように変更
PROFILER_OUT_DIR="${OUTPUT_DIR}/isg_profiler_res"

# NOTE: ${PER_SPECIES_OPT} is intentionally unquoted to allow it to be empty
python3 -m quant_normalizer \
  --reference_dir "${REF_DIR}" \
  --sample_metadata "${SAMPLE_METADATA}" \
  --sf_dir "${OUTPUT_DIR}" \
  --out_dir "${PROFILER_OUT_DIR}" \
  ${PER_SPECIES_OPT}
