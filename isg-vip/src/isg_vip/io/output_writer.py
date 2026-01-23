# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Hiroaki Unno & Jumpei Ito

import logging
from pathlib import Path

import pandas as pd
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


def write_to_csv(df: DataFrame, output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    logger.info(f"CSV file was exported to {output_path}")
    df.to_csv(
        output_path,
        index=False,
    )


def export_final_prediction(
    dir_name: Path,
    merged_df: DataFrame,
    all_dfs: list,
):
    for df in all_dfs:
        merged_df = pd.merge(merged_df, df, on="ID", how="inner")

    label_cols = [f"Prediction_score_fold{n}" for n in range(5)]
    merged_df["Final_Prediction_score(mean)"] = merged_df[label_cols].mean(axis=1)
    label_cols = [f"Prediction_Label_fold{n}" for n in range(5)]
    merged_df["Final_Prediction_Label"] = merged_df[label_cols].mode(axis=1)[0]

    # Save final results
    merged_df = merged_df.sort_values(by="ID")
    base_name = "Infection_Prediction_Stacking_all"
    write_to_csv(merged_df, dir_name / f"{base_name}.csv")

    final_df = merged_df[["ID", "Final_Prediction_score(mean)", "Final_Prediction_Label"]]
    final_df = final_df.rename(
        columns={
            "Final_Prediction_score(mean)": "Prediction_score(mean)",
            "Final_Prediction_Label": "Prediction_Label",
        }
    )

    base_name_final = "Infection_Prediction_Stacking_final"
    write_to_csv(final_df, dir_name / f"{base_name_final}.csv")
    write_to_tsv(final_df, dir_name / f"{base_name_final}.txt")
