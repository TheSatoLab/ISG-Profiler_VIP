# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Luca Nishimura & Jumpei Ito

import logging
from pathlib import Path

from pandas import DataFrame

logger = logging.getLogger(__name__)


def write_to_tsv(df: DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"TSV file was exported to {output_path}")
    df.to_csv(
        output_path,
        sep="\t",
        index=False,
    )
