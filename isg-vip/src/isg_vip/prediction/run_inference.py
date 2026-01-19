# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Hiroaki Unno & Jumpei Ito

from pathlib import Path

import pandas as pd

from isg_vip.io.model_loader import ISGModelArtifacts, ModelType
from isg_vip.io.output_writer import write_to_csv
from isg_vip.prediction.ensemble import (
    get_train_categories_from_encoder,
    prediction,
    replace_unseen_categories,
)


def predict_each_fold(
    artifacts: ISGModelArtifacts,
    dir_name: Path,
    columns_to_process: list[str],
    exclude_columns: list[str],
    meta_info_: pd.DataFrame,
    copiedX: pd.DataFrame,
):
    """Predict specific times"""
    for n in range(5):
        # Normalize
        normalizer = artifacts.get_normalizer(n)
        X_test = normalizer.transform(copiedX)

        # LightGBM prediction
        final_model_lgb = artifacts.get_model(ModelType.LGB, n)
        y_test_pred_label_lgb, y_test_pred_prob_lgb = prediction(
            artifacts,
            final_model_lgb,
            columns_to_process,
            exclude_columns,
            X_test,
            ModelType.LGB,
            n,
        )

        # Logistic Regression prediction
        final_model_lr = artifacts.get_model(ModelType.LR, n)
        y_test_pred_label_lr, y_test_pred_prob_lr = prediction(
            artifacts,
            final_model_lr,
            columns_to_process,
            exclude_columns,
            X_test,
            ModelType.LR,
            n,
        )

        # Ensure consistency by mapping out-of-vocabulary labels to 'unknown'
        encoder = artifacts.get_encoder(ModelType.LGB, n)
        train_categories = get_train_categories_from_encoder(encoder, columns_to_process)
        X_test = replace_unseen_categories(
            X_test, columns=columns_to_process, train_categories=train_categories
        )
        missing_values = X_test[X_test.isnull().any(axis=1)].index
        X_test_final = X_test.drop(index=missing_values)

        # Build results DataFrame
        df_long = pd.DataFrame(
            {
                "ID": list(X_test_final["ID"]) * 2,
                "Model": ["LightGBM"] * len(X_test_final)
                + ["LogisticRegression"] * len(X_test_final),
                "Prediction_score": list(y_test_pred_prob_lgb) + list(y_test_pred_prob_lr),
                "Prediction_Label": list(y_test_pred_label_lgb) + list(y_test_pred_label_lr),
            }
        )
        df_long["Prediction_Label"] = df_long["Prediction_Label"].replace(
            {0: "Negative", 1: "Positive"}
        )

        # Merge & save
        merged_df = pd.merge(meta_info_, df_long, on="ID", how="inner")
        write_to_csv(merged_df, dir_name / f"Infection_Prediction_{n}.csv")
