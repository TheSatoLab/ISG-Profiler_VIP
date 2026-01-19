#!/usr/bin/env python
# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Luca Nishimura & Jumpei Ito


import argparse
from pathlib import Path

import pandas as pd

from quant_normalizer import __version__
from quant_normalizer.core.isg_scorer import calculate_isg_scores
from quant_normalizer.core.sample_processor import process_samples
from quant_normalizer.io.output_writer import write_to_tsv
from quant_normalizer.io.reference_loader import load_reference_data
from quant_normalizer.utils.logger import parse_args_as_log_level, setup_logger


def parse_args():
    """Parse arguments"""
    parser = argparse.ArgumentParser(
        prog="quant_normalizer",
        description="ISG Profiler: "
        "Compute normalized and standardized ISG profiles from Salmon quant.sf files."
    )

    parser.add_argument("--version", action="version", version=f"%(prog)s {__version__}")

    parser.add_argument(
        "--reference_dir",
        required=True,
        help="Directory containing reference files "
        "(gene2refseq list, Aves/Mars removal lists, mean & SD tables, Amniota398_sp_id.list).",
    )
    parser.add_argument(
        "--sample_metadata",
        required=True,
        help="Path to sample metadata table (TSV). "
        "Must include: sample_id, species_host, order_host, clade_host.",
    )
    parser.add_argument(
        "--sf_dir",
        required=True,
        help="Directory containing Salmon quant.sf files (named as <sample_id>_quant.sf).",
    )
    parser.add_argument(
        "--out_dir",
        required=True,
        help="Output directory (example: isg_profiler_out/result).",
    )
    parser.add_argument(
        "--per_species",
        action="store_true",
        help="If set, group counts by both hum_symbol and tax_id and attach species. "
        "ISG_score will NOT be computed in this mode.",
    )

    parser.add_argument(
        "--log_level",
        choices=["info", "debug", "warning"],
        default="info",
        help="Log level",
    )

    args = parser.parse_args()

    reference_dir = Path(args.reference_dir)
    sample_metadata_path = Path(args.sample_metadata)
    sf_dir = Path(args.sf_dir)
    out_dir = Path(args.out_dir)
    per_species = args.per_species

    log_level = args.log_level

    return (reference_dir, sample_metadata_path, sf_dir, out_dir, per_species, log_level)


def main():
    (reference_dir, sample_metadata_path, sf_dir, out_dir, per_species, log_level) = parse_args()
    logger = setup_logger(None, level=parse_args_as_log_level(log_level))
    out_dir.mkdir(parents=True, exist_ok=True)
    ref = load_reference_data(reference_dir, sample_metadata_path, per_species)

    all_sample_gene_count_list = process_samples(
        sample_metadata=ref.sample_metadata,
        gene_info=ref.gene_info,
        sf_dir=sf_dir,
        aves_neg_genes=ref.aves_neg_genes,
        mars_neg_genes=ref.mars_neg_genes,
        gene_mean_sd_list=ref.gene_mean_sd_list,
        per_species=per_species,
    )

    # Combine outputs
    if len(all_sample_gene_count_list) > 0:
        all_sample_gene_count_df = pd.concat(all_sample_gene_count_list, ignore_index=True)
    else:
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
        all_sample_gene_count_df = pd.DataFrame(columns=base_cols)

    # If per_species, attach species information via tax_id
    if (
        per_species
        and "tax_id" in all_sample_gene_count_df.columns
        and ref.species_map_df is not None
    ):
        all_sample_gene_count_df = all_sample_gene_count_df.merge(
            ref.species_map_df, on="tax_id", how="left"
        )  # adds "species" column

    # Save per-gene results
    write_to_tsv(
        all_sample_gene_count_df,
        out_dir / "per_gene_count.tsv",
    )

    # --------------------------------------------------------
    # Compute ISG score (mean of standardized ISGs)
    # NOTE: per_species mode -> do NOT compute ISG_score
    # --------------------------------------------------------
    if per_species:
        logger.debug("per_species mode: skipped to generate ISG_score.tsv")
    else:
        isg_score_df = calculate_isg_scores(all_sample_gene_count_df, ref.sample_metadata)

        # Save ISG scores
        write_to_tsv(
            isg_score_df,
            out_dir / "ISG_score.tsv",
        )


if __name__ == "__main__":
    main()
