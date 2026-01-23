#!/usr/bin/env bash
# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Luca Nishimura & Jumpei Ito

set -eu

# === ARGUMENTS ===
# Default arguments:
THREAD=""
FASTQ_DIR=""
OUTPUT_DIR=""
SALMON_INDEX_PATH=""
SALMON_BIN="salmon" # Default to using command in PATH

usage() {
  echo "Usage: $0 --thread <int> --fastq_dir <path> --out_dir <path> --index <path> [--salmon_bin <path>]"
  exit 1
}

# parse arguments
while [[ $# -gt 0 ]]; do
  key="$1"
  case "$key" in
  --thread | --fastq_dir | --out_dir | --salmon_index | --salmon_bin)
    if [[ -z "${2:-}" ]] || [[ "${2:-}" == --* ]]; then
      echo "Error: Argument for $key is missing"
      usage
    fi

    case "$key" in
    --thread) THREAD="$2" ;;
    --fastq_dir) FASTQ_DIR="$2" ;;
    --out_dir) OUTPUT_DIR="$2" ;;
    --salmon_index) SALMON_INDEX_PATH="$2" ;;
    --salmon_bin) SALMON_BIN="$2" ;;
    esac
    shift 2
    ;;
  *)
    echo "Error: Unknown argument $1"
    usage
    ;;
  esac
done

# Check required arguments:
if [[ -z "$THREAD" ]] || [[ -z "$FASTQ_DIR" ]] || [[ -z "$OUTPUT_DIR" ]] || [[ -z "$SALMON_INDEX_PATH" ]]; then
  echo "Error: Missing required arguments."
  usage
fi
# ==================

TMP_DIR="${OUTPUT_DIR}_tmp"
REF="${SALMON_INDEX_PATH}"

mkdir -p "${OUTPUT_DIR}" "${TMP_DIR}"

run_salmon() {
  local ID=$1
  local SAMPLE_TMP="${TMP_DIR}/${ID}"
  local FQ1="${FASTQ_DIR}/${ID}_1.cleaned.fastq"
  local FQ2="${FASTQ_DIR}/${ID}_2.cleaned.fastq"
  local FQ_SINGLE="${FASTQ_DIR}/${ID}.cleaned.fastq"

  echo ">>> Processing ID: ${ID}"
  mkdir -p "${SAMPLE_TMP}"

  # common options to be passed to `salmon quant`
  local salmon_opts=(
    "-i" "${REF}"
    "-l" "A"
    "-p" "${THREAD}"
    "--validateMappings"
    "-o" "${SAMPLE_TMP}"
  )

  # different options by fastq file
  if [[ -f "$FQ1" && -f "$FQ2" ]]; then
    echo "Mode: Paired-end"
    salmon_opts+=("-1" "$FQ1" "-2" "$FQ2")

  elif [[ -f "$FQ_SINGLE" ]]; then
    echo "Mode: Single-end"
    salmon_opts+=("-r" "$FQ_SINGLE")

  else
    echo "Warning: Files for ${ID} not found. Skipping."
    return 1
  fi

  # Call salmon using the variable
  "${SALMON_BIN}" quant "${salmon_opts[@]}"

  if [ -f "${SAMPLE_TMP}/quant.sf" ]; then
    cp "${SAMPLE_TMP}/quant.sf" "${OUTPUT_DIR}/${ID}_quant.sf"
  fi
  rm -r "${SAMPLE_TMP}"
}

IFS=$'\n'
IDS=($(find "${FASTQ_DIR}" -maxdepth 1 -name "*.cleaned.fastq" -exec basename {} \; |
  sed -E 's/(_[12])?\.cleaned\.fastq//' |
  sort -u))
unset IFS

# Check array is empty (＝fail find or no file) or not
if [ ${#IDS[@]} -eq 0 ]; then
  echo "Error: No valid fastq files found in '${FASTQ_DIR}'." >&2
  exit 1
fi

for ID in "${IDS[@]}"; do
  run_salmon "${ID}"
done
