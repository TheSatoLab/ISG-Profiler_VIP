# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Hiroaki Unno & Jumpei Ito

from pathlib import Path

import numpy as np
import pandas as pd

from isg_vip.io.constants import (
    GeneType,
    LoadedPerGeneCountTsvCols,
    PerGeneCountTsvCols,
)


def load_per_gene_count(info_file_path: Path, gene_list_path: Path):
    """Load data and filter, then normalize"""
    info = pd.read_csv(info_file_path, sep="\t")[
        [
            PerGeneCountTsvCols.SAMPLE_ID,
            PerGeneCountTsvCols.HUM_SYMBOL,
            PerGeneCountTsvCols.RAW_COUNT,
            PerGeneCountTsvCols.TYPE,
        ]
    ]
    info = info.rename(columns={PerGeneCountTsvCols.SAMPLE_ID: PerGeneCountTsvCols.ID})
    info_filled = _zero_filling_missing_genes(gene_list_path, info)

    # Calcurate following columns:
    # all_sum: use raw_count. group by sample_id, then sum grouped values
    # cntl_sum: use raw_count. group by sample_id and filtered by type == 'cntl' then sum values
    info_filled = info_filled.assign(
        **{
            LoadedPerGeneCountTsvCols.ALL_SUM: info_filled.groupby(PerGeneCountTsvCols.ID)[
                PerGeneCountTsvCols.RAW_COUNT
            ].transform("sum"),
            LoadedPerGeneCountTsvCols.CNTL_SUM: info_filled.groupby(PerGeneCountTsvCols.ID)[
                PerGeneCountTsvCols.RAW_COUNT
            ].transform(
                lambda x: x[
                    info_filled.loc[x.index, PerGeneCountTsvCols.TYPE] == GeneType.CNTL
                ].sum()
            ),
        }
    )
    # WARNING: DO NOT CHANGE THRESHOLD as it was used during the model creation.
    info_filled = info_filled[info_filled[LoadedPerGeneCountTsvCols.CNTL_SUM] > 10000]

    # Normalization
    info_filled = info_filled.assign(
        **{
            LoadedPerGeneCountTsvCols.NORM_CNTL: lambda x: x[PerGeneCountTsvCols.RAW_COUNT]
            / x[LoadedPerGeneCountTsvCols.CNTL_SUM],
            LoadedPerGeneCountTsvCols.NORM_CNTL_LOG: lambda x: np.log2(
                x[LoadedPerGeneCountTsvCols.NORM_CNTL] * (10e5) + 1
            ),
            LoadedPerGeneCountTsvCols.NORM_ALL: lambda x: x[PerGeneCountTsvCols.RAW_COUNT]
            / x[LoadedPerGeneCountTsvCols.ALL_SUM],
            LoadedPerGeneCountTsvCols.NORM_ALL_LOG: lambda x: np.log2(
                x[LoadedPerGeneCountTsvCols.NORM_ALL] * (10e5) + 1
            ),
        }
    )

    return info_filled


def _zero_filling_missing_genes(gene_list_path: Path, info: pd.DataFrame):
    """
    Zero padding; ISG-Profiler result does not have all genes. Some non-exist genes are missings.
    Check for removed gene list:
    Aves:         reference/Aves_rem.txt
    Marsupialia:  reference/Mars_rem.txt
    :param gene_list_path: all gene list
    :type gene_list_path: Path
    :param info: zero padding target
    :type info: pd.DataFrame
    """
    if info.empty:
        raise ValueError("'info' is empty. No zero padding targets")

    gene_df = pd.read_csv(gene_list_path, header=None)
    if gene_df.empty:
        raise ValueError(f"Invalid gene list {gene_list_path}.")

    # Create combination of sample_ID and all genes
    target_genes = gene_df[0].unique().tolist()
    samples = info[PerGeneCountTsvCols.ID].unique().tolist()
    # Create a Cartesian product of all samples and target genes.
    index = pd.MultiIndex.from_product(
        [samples, target_genes],
        names=[PerGeneCountTsvCols.ID, PerGeneCountTsvCols.HUM_SYMBOL],
    )
    # Convert the MultiIndex into a DataFrame with regular columns to facilitate merging.
    all_combinations = pd.DataFrame(index=index).reset_index()
    info_filled = pd.merge(
        all_combinations,
        info,
        on=[PerGeneCountTsvCols.ID, PerGeneCountTsvCols.HUM_SYMBOL],
        how="left",
    )

    fill_values = {
        PerGeneCountTsvCols.RAW_COUNT: 0,  # Zero filling NaN value.
        PerGeneCountTsvCols.TYPE: GeneType.ISG,  # Zero-filled gene type is ISG.
    }
    info_filled = info_filled.fillna(value=fill_values)

    return info_filled
