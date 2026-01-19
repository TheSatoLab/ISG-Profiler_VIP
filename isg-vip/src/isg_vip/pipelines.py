# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Hiroaki Unno & Jumpei Ito

import argparse
import logging
import os
import sys
import warnings
from pathlib import Path

import pandas as pd
from pandas.errors import PerformanceWarning
from sklearn.exceptions import (
    # NOTE: This class contains, but cause IDE error
    InconsistentVersionWarning,  # type: ignore
)

from isg_vip.io.constants import (
    LoadedPerGeneCountTsvCols,
    MetadataTsvCols,
    PerGeneCountTsvCols,
)

# HACK: Add to __main__ namespace to import .pkl model
from isg_vip.prediction.ensemble import (
    CustomNormalizer,
)

setattr(sys.modules["__main__"], "CustomNormalizer", CustomNormalizer)
from isg_vip import __version__  # noqa: E402
from isg_vip.io.data_loader import load_per_gene_count  # noqa: E402
from isg_vip.io.model_loader import ISGModelArtifacts  # noqa: E402
from isg_vip.io.output_writer import export_final_prediction  # noqa: E402
from isg_vip.prediction.ensemble import (  # noqa: E402
    execute_prediction,
)
from isg_vip.prediction.run_inference import predict_each_fold  # noqa: E402
from isg_vip.utils.logger import setup_logger  # noqa: E402

PACKAGE_ROOT = Path(__file__).resolve().parent
print(PACKAGE_ROOT)
MODEL_DIR = PACKAGE_ROOT / "model_dir"
GENE_LIST_PATH = PACKAGE_ROOT / "reference" / "gene_list.txt"

CWD = Path.cwd()
INPUT_DIR = CWD / "input"
OUTPUT_DIR = CWD / "output"


def parse_args():
    parser = argparse.ArgumentParser(
        prog="isg_vip",
        description=(
            "ISG VIP: Machine learning-based viral infection prediction using "
            "ISG (Interferon-Stimulated Gene) expression profiles with ensemble methods."
        ),
    )

    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    parser.add_argument(
        "--gene_count_file",
        required=False,
        default=INPUT_DIR / "per_gene_count.tsv",
        help="per_gene_count.tsv. See README.md file.",
    )

    parser.add_argument(
        "--metadata",
        required=False,
        default=INPUT_DIR / "sample_metadata.tsv",
        help="sample_metadata.tsv. See README.md file.",
    )

    parser.add_argument(
        "--output",
        required=False,
        default=OUTPUT_DIR,
        help="Output directory.",
    )

    args = parser.parse_args()
    gene_count_file = Path(args.gene_count_file).resolve()
    metadata_file = Path(args.metadata).resolve()
    output_dir = Path(args.output)
    # Output directory
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    return (gene_count_file, metadata_file, output_dir)


def main():
    (gene_count_file, metadata_file, output_dir) = parse_args()

    logger = setup_logger(None, level=logging.INFO)
    # NOTE: ignore, this code do not concern about performance
    warnings.simplefilter("ignore", category=PerformanceWarning)
    # NOTE: this cause future warnings:
    # if you want to update pandas,
    warnings.simplefilter("ignore", category=FutureWarning)
    # NOTE: This project model not match scikit-learn version.
    # Already confirmed output value is same in 1.5.1.
    # TODO: When updating model, remove this supress.
    warnings.simplefilter("ignore", category=InconsistentVersionWarning)

    # Files load
    try:
        artifacts = ISGModelArtifacts.from_directory(MODEL_DIR)
    except FileNotFoundError as e:
        logger.critical(e)
        exit(1)

    # Load data and filter
    info = load_per_gene_count(gene_count_file, GENE_LIST_PATH)

    # Create pivot table
    info_ = info.groupby(PerGeneCountTsvCols.ID)[[LoadedPerGeneCountTsvCols.ALL_SUM]].mean()
    gene_list = sorted(info[PerGeneCountTsvCols.HUM_SYMBOL].unique())
    info_ = pd.concat([info_, pd.DataFrame({g: 0 for g in gene_list}, index=info_.index)], axis=1)
    pivot_data = info.pivot_table(
        index=MetadataTsvCols.ID,
        columns=PerGeneCountTsvCols.HUM_SYMBOL,
        values=LoadedPerGeneCountTsvCols.NORM_CNTL_LOG,
        aggfunc="sum",
    )
    info_.update(pivot_data)

    # Merge with metadata
    meta_info = pd.read_csv(metadata_file, sep="\t")[
        [MetadataTsvCols.SAMPLE_ID, MetadataTsvCols.SPECIES_HOST]
    ]
    meta_info = meta_info.rename(
        columns={
            MetadataTsvCols.SAMPLE_ID: MetadataTsvCols.ID,
            MetadataTsvCols.SPECIES_HOST: MetadataTsvCols.H_SPECIES,
        }
    )

    df2 = pd.merge(info_, meta_info, on=PerGeneCountTsvCols.ID)

    # Species filter
    mbio_all = df2.copy()

    # Additional metadata
    group_info = pd.read_csv(metadata_file, sep="\t")
    meta_info_ = (
        group_info[
            [MetadataTsvCols.SAMPLE_ID, MetadataTsvCols.SPECIES_HOST, MetadataTsvCols.ORDER_HOST]
        ]
        .copy()
        .rename(
            columns={
                MetadataTsvCols.SAMPLE_ID: MetadataTsvCols.ID,
                MetadataTsvCols.SPECIES_HOST: MetadataTsvCols.H_SPECIES,
                MetadataTsvCols.ORDER_HOST: MetadataTsvCols.ORDER,
            }
        )
    )

    # Features: Use all columns of imported records.
    X = mbio_all.copy()

    columns_to_process = [MetadataTsvCols.H_SPECIES, MetadataTsvCols.ORDER]  # Target columns
    exclude_columns = [
        MetadataTsvCols.ID,
        MetadataTsvCols.H_SPECIES,
        MetadataTsvCols.ORDER,
    ]  # Columns to exclude
    predict_each_fold(
        artifacts,
        output_dir,
        columns_to_process,
        exclude_columns,
        meta_info_=meta_info_,
        copiedX=X,
    )

    all_dfs = execute_prediction(
        artifacts,
        output_dir,
        columns_to_process,
        exclude_columns,
        X.copy(),
        meta_info_,
    )
    export_final_prediction(output_dir, meta_info_.copy(), all_dfs)


if __name__ == "__main__":
    main()
