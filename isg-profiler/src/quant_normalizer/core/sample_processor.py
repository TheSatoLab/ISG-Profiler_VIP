# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Luca Nishimura & Jumpei Ito

import logging
from pathlib import Path

import numpy as np
from pandas import DataFrame, read_csv

logger = logging.getLogger(__name__)


def summarize_sf_for_sample(
    sample_id: str,
    clade_host: str,
    gene_info: DataFrame,
    sf_dir: Path,
    aves_neg_genes: set,
    mars_neg_genes: set,
    gene_mean_sd_list: DataFrame,
    per_species: bool = False,
) -> DataFrame:
    """
    Process a single Salmon quant.sf file and return per-gene statistics per sample.

    If per_species is False:
        group by hum_symbol
    If per_species is True:
        group by (hum_symbol, tax_id)
    """

    sf_path = sf_dir / f"{sample_id}_quant.sf"

    # If the sample does not exist, return an empty dataframe silently
    if not sf_path.exists():
        logger.warning(
            f"{sf_path} was not found. check sample_metadata.example.tsv file and fastaq files."
        )
        base_cols = [
            "sample_id",
            "hum_symbol",
            "raw_count",
            "type",
            "normalized_count",
            "standardized_count",
        ]
        if per_species:
            base_cols.insert(2, "tax_id")  # sample_id, hum_symbol, tax_id, ...
        return DataFrame(columns=base_cols)

    # Load Salmon quant file and filter NumReads > 0
    sf_data_temp = read_csv(sf_path, sep="\t")
    # sf_data_temp = sf_data_temp.loc[sf_data_temp["NumReads"] > 0].copy()

    # Join with reference annotation (gene_info)
    merged = sf_data_temp.merge(
        gene_info,
        how="left",
        left_on="Name",
        right_on="Isoform",
    )

    # Decide grouping columns depending on per_species flag
    if per_species:
        group_cols = ["hum_symbol", "tax_id"]
    else:
        group_cols = ["hum_symbol"]

    # Aggregate raw read counts
    grouped = (
        merged.dropna(subset=group_cols)
        .groupby(group_cols, as_index=False)
        .agg(
            raw_count=("NumReads", "sum"),
            type=("type", "first"),
        )
    )

    # Remove negative-control genes depending on host clade
    if clade_host == "Aves":
        grouped = grouped.loc[~grouped["hum_symbol"].isin(aves_neg_genes)].copy()
    elif clade_host == "Marsupialia":
        grouped = grouped.loc[~grouped["hum_symbol"].isin(mars_neg_genes)].copy()

    # Count normalization using control genes
    cntl_total = grouped.loc[grouped["type"] == "cntl", "raw_count"].sum()

    if cntl_total > 0:
        scaled = grouped["raw_count"] / cntl_total
        grouped["normalized_count"] = np.log2(scaled * 1000000 + 1)
    else:
        grouped["normalized_count"] = np.nan

    # Merge mean/sd reference table for standardization (by hum_symbol)
    grouped = grouped.merge(
        gene_mean_sd_list[["hum_symbol", "mean_norm_genes_log", "sd_norm_genes_log"]],
        on="hum_symbol",
        how="left",
    )

    # Compute standardized value: (normalized - mean) / sd
    grouped["standardized_count"] = (
        grouped["normalized_count"] - grouped["mean_norm_genes_log"]
    ) / grouped["sd_norm_genes_log"]

    grouped = grouped.drop(columns=["mean_norm_genes_log", "sd_norm_genes_log"])

    # Insert sample_id column at the beginning
    grouped.insert(0, "sample_id", sample_id)

    return grouped


def process_samples(
    sample_metadata: DataFrame,
    gene_info: DataFrame,
    sf_dir: Path,
    aves_neg_genes: set[str],
    mars_neg_genes: set[str],
    gene_mean_sd_list: DataFrame,
    per_species: bool,
):
    """Process each sample"""
    all_sample_gene_count_list = []
    for _, row in sample_metadata.iterrows():
        sample_id = row["sample_id"]
        clade_host = row["clade_host"]
        logger.debug(f"processing sample_id: {sample_id} clade_host: {clade_host}")

        sample_df = summarize_sf_for_sample(
            sample_id=sample_id,
            clade_host=clade_host,
            gene_info=gene_info,
            sf_dir=sf_dir,
            aves_neg_genes=aves_neg_genes,
            mars_neg_genes=mars_neg_genes,
            gene_mean_sd_list=gene_mean_sd_list,
            per_species=per_species,
        )

        if sample_df.empty:
            logger.warning(f"Empty sample: {sample_id} (per_species={per_species})")
        else:
            logger.info(f"Processed sample: {sample_id} (per_species={per_species})")
            all_sample_gene_count_list.append(sample_df)

    return all_sample_gene_count_list
