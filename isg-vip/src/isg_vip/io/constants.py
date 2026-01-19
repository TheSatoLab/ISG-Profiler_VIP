# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Hiroaki Unno & Jumpei Ito


class PerGeneCountTsvCols:
    """
    per_gene_count.tsv column names
    """

    SAMPLE_ID = "sample_id"
    HUM_SYMBOL = "hum_symbol"
    RAW_COUNT = "raw_count"
    TYPE = "type"
    NORMALIZED_COUNT = "normalized_count"
    STANDARDIZED_COUNT = "standardized_count"

    ID = "ID"
    """`sample_id` is renamed to this string."""

    ALL_COLUMNS = [
        SAMPLE_ID,
        HUM_SYMBOL,
        RAW_COUNT,
        TYPE,
        NORMALIZED_COUNT,
        STANDARDIZED_COUNT,
    ]


class LoadedPerGeneCountTsvCols:
    ALL_SUM = "all_sum"
    CNTL_SUM = "cntl_sum"

    NORM_ALL = "norm_all"
    NORM_CNTL = "norm_cntl"
    NORM_CNTL_LOG = "norm_cntl_log"
    NORM_ALL_LOG = "norm_all_log"


class MetadataTsvCols:
    """
    sample_metadata.tsv column names
    """

    SAMPLE_ID = "sample_id"
    SPECIES_HOST = "species_host"
    ORDER_HOST = "order_host"
    CLADE_HOST = "clade_host"

    ID = "ID"
    """`sample_id` is renamed to this string."""
    ORDER = "order"
    """`order_host` is renamed to this string."""
    H_SPECIES = "h_species"
    """`Host_species` is renamed to this string."""

    ALL_COLUMNS = [
        SAMPLE_ID,
        SPECIES_HOST,
        ORDER_HOST,
        CLADE_HOST,
    ]


class GeneType:
    ISG = "ISG"
    CNTL = "cntl"
