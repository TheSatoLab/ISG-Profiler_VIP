# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Luca Nishimura & Jumpei Ito

from dataclasses import dataclass
from pathlib import Path
from typing import Final, Optional

import pandas as pd


class ReferenceFiles:
    """
    Define reference file names
    """

    GENE2REFSEQ: Final[str] = "gene2refseq_Amniota_ISGcntl_241003.list"
    AVES_REM: Final[str] = "Aves_rem.txt"
    MARS_REM: Final[str] = "Mars_rem.txt"
    AVES_MAM_MBIO_ISG_CNTL_MNSD: Final[str] = "Aves_Mam_mbio_ISGcntl_mnsd.txt"
    SP_ID_LIST: Final[str] = "Amniota398_sp_id.list"


@dataclass(frozen=True)
class ReferenceData:
    """Hold loaded reference files data"""

    gene_info: pd.DataFrame
    gene_mean_sd_list: pd.DataFrame
    sample_metadata: pd.DataFrame
    aves_neg_genes: set[str]
    mars_neg_genes: set[str]
    species_map_df: Optional[pd.DataFrame] = None


def _read_tsv(path: Path) -> pd.DataFrame:
    """
    Read TSV (\t) separated data

    :param path: tsv file path
    :type path: Path
    :return: tsv file DataFrame
    :rtype: DataFrame
    """
    return pd.read_csv(path, sep="\t")


def _load_gene_set(path: Path, column: str) -> set[str]:
    """
    Read unique gene set from tsv file.

    :param path: target file path to read
    :type path: Path
    :param column: which one to read
    :type column: str
    :return: loaded gene set
    :rtype: set[str]
    """
    df = _read_tsv(path)
    return set(df[column].dropna())


def load_reference_data(
    reference_dir: Path, sample_metadata_path: Path, per_species: bool
) -> ReferenceData:
    """load reference directory data files"""
    # --------------------------------------------------------
    # Load reference files
    # --------------------------------------------------------
    gene_info_path = reference_dir / ReferenceFiles.GENE2REFSEQ
    gene_mean_sd_list_path = reference_dir / ReferenceFiles.AVES_MAM_MBIO_ISG_CNTL_MNSD

    gene_info = _read_tsv(gene_info_path)
    gene_mean_sd_list = _read_tsv(gene_mean_sd_list_path)

    # Load sample metadata
    sample_metadata = _read_tsv(sample_metadata_path)

    # Load negative gene lists
    aves_neg_genes = _load_gene_set(reference_dir / ReferenceFiles.AVES_REM, column="hum_symbol")
    mars_neg_genes = _load_gene_set(reference_dir / ReferenceFiles.MARS_REM, column="hum_symbol")

    # If per_species mode, load tax_id -> species mapping
    species_map_df = None
    if per_species:
        species_map_path = reference_dir / ReferenceFiles.SP_ID_LIST
        species_map_df = _read_tsv(species_map_path)  # columns: tax_id, species

    return ReferenceData(
        gene_info=gene_info,
        gene_mean_sd_list=gene_mean_sd_list,
        sample_metadata=sample_metadata,
        aves_neg_genes=aves_neg_genes,
        mars_neg_genes=mars_neg_genes,
        species_map_df=species_map_df,
    )
