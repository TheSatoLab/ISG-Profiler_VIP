# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Hiroaki Unno & Jumpei Ito

from pathlib import Path
from typing import Any

import lightgbm as lgb
import numpy as np
import pandas as pd
from sklearn.base import BaseEstimator, TransformerMixin

from isg_vip.io.model_loader import ISGModelArtifacts, ModelType
from isg_vip.io.output_writer import write_to_csv
from isg_vip.preprocessing.normalizer import cal_z


class CustomNormalizer(BaseEstimator, TransformerMixin):
    def __init__(self, isg_list, group_info):
        self.isg_list = isg_list
        self.group_info = group_info
        self.mean_ = None
        self.std_ = None
        self.gene_means_ = {}
        self.gene_stds_ = {}

    def transform(self, X):
        X = X.copy()
        # normalize total expression
        X["sum_all_log"] = X["all_sum"].apply(lambda x: np.log2(x + 1))
        X["sum_all_log_z"] = (X["sum_all_log"] - self.mean_) / self.std_
        X.drop(columns=["all_sum", "sum_all_log"], inplace=True)

        # standardize each gene
        for col in self.gene_means_:
            X[col] = (X[col] - self.gene_means_[col]) / self.gene_stds_[col]

        # map species to order
        group_info_ = self.group_info[self.group_info["species"].isin(X["h_species"])]
        mapper = group_info_.drop_duplicates("species").set_index("species")["order"]
        X["order"] = X["h_species"].map(mapper)

        # calculate mean ISG score
        X["mean_ISGscore"] = X[self.isg_list].mean(axis=1)
        return X


def get_train_categories_from_encoder(encoder, columns):
    """Extract learned categories from the encoder and map them to column names."""
    categories = {}
    for col, cats in zip(columns, encoder.categories_):
        categories[col] = list(cats)
    return categories


def replace_unseen_categories(data, columns: list[str], train_categories):
    data = data.copy()
    for col in columns:
        valid = set(train_categories[col])
        data[col] = data[col].where(data[col].isin(valid), "unknown")
    return data


def prediction(
    artifacts: ISGModelArtifacts,
    final_model: Any,
    columns_to_process: list[str],
    exclude_columns: list[str],
    X_test: pd.DataFrame,
    m_type: ModelType,
    fold: int,
):
    encoder = artifacts.get_encoder(m_type, fold)
    train_categories = get_train_categories_from_encoder(encoder, columns_to_process)
    X_test = replace_unseen_categories(
        X_test, columns=columns_to_process, train_categories=train_categories
    )

    # Encode categorical columns
    X_test_encoded = encoder.transform(X_test[columns_to_process])
    X_test_encoded = pd.DataFrame(
        X_test_encoded,
        index=X_test.index,
        columns=encoder.get_feature_names_out(columns_to_process),
    )

    # Combine with numeric columns
    X_test_numeric = X_test.drop(columns=exclude_columns)
    X_test_final = pd.concat([X_test_numeric, X_test_encoded], axis=1)

    # Drop rows with missing values
    missing_values = X_test_final[X_test_final.isnull().any(axis=1)].index
    X_test_final = X_test_final.drop(index=missing_values)

    # Predict probabilities
    if isinstance(final_model, lgb.Booster):
        y_test_pred_prob = final_model.predict(X_test_final)
    else:
        y_test_pred_prob = final_model.predict_proba(X_test_final)[:, 1]

    # Apply threshold to get labels
    best_threshold_train = artifacts.get_threshold(m_type)[fold]
    y_test_pred_label = (y_test_pred_prob >= best_threshold_train).astype(int)

    return y_test_pred_label, y_test_pred_prob


def train_stacking_model(
    artifacts: ISGModelArtifacts,
    dir_name: Path,
    columns_to_process: list[str],
    exclude_columns: list[str],
    X_test: pd.DataFrame,
    n: int,
):
    # Create meta features for each fold
    X_test_meta = (
        pd.read_csv(dir_name / f"Infection_Prediction_{n}.csv")
        .pivot(index="ID", columns="Model", values="Prediction_score")
        .fillna(0)
    )
    for model in ["LightGBM", "LogisticRegression"]:
        X_test_meta = cal_z(X_test_meta.astype(float), model)

    encoder = artifacts.get_encoder(ModelType.META, n)
    train_categories = get_train_categories_from_encoder(encoder, columns_to_process)
    X_test = replace_unseen_categories(
        X_test, columns=columns_to_process, train_categories=train_categories
    )
    missing_values = X_test[X_test.isnull().any(axis=1)].index
    X_test_final = X_test.drop(index=missing_values)

    # Merge with additional features
    X_test_meta = X_test_meta.merge(X_test_final[["ID", "h_species", "order"]], on="ID", how="left")

    # Run prediction with trained meta model
    final_model_meta = artifacts.get_model(ModelType.META, n)
    y_test_pred_label_meta, y_test_pred_prob_meta = prediction(
        artifacts,
        final_model_meta,
        columns_to_process,
        exclude_columns,
        X_test_meta,
        ModelType.META,
        n,
    )

    return final_model_meta, y_test_pred_label_meta, y_test_pred_prob_meta, X_test_meta


def execute_prediction(
    artifacts: ISGModelArtifacts,
    dir_name: Path,
    columns_to_process: list[str],
    exclude_columns: list[str],
    copiedX: pd.DataFrame,
    meta_info_: pd.DataFrame,
):
    # Save results for each fold
    for n in range(5):
        normalizer = artifacts.get_normalizer(n)
        X_test = normalizer.transform(copiedX)
        meta_model, y_test_pred_meta, y_test_pred_prob_meta, X_test_meta = train_stacking_model(
            artifacts,
            dir_name,
            columns_to_process,
            exclude_columns,
            X_test,
            n,
        )

        df_long = pd.DataFrame(
            {
                "ID": X_test_meta["ID"],
                "Model": "Stacking",
                "Prediction_score": y_test_pred_prob_meta,
                "Prediction_Label": y_test_pred_meta,
            }
        )
        df_long["Prediction_Label"] = df_long["Prediction_Label"].replace(
            {0: "Negative", 1: "Positive"}
        )

        merged_df = pd.merge(meta_info_, df_long, on="ID", how="inner").sort_values(by="ID")
        write_to_csv(merged_df, dir_name / f"Infection_Prediction_Stacking_{n}_external.csv")

    # TODO: ひとつ前のfor loopと内容が同じなので統合可能なはず。要確認。
    # Combine all folds and create final label by majority vote
    all_dfs = []
    for n in range(5):
        normalizer = artifacts.get_normalizer(n)
        X_test = normalizer.transform(copiedX)
        meta_model, y_test_pred_meta, y_test_pred_prob_meta, X_test_meta = train_stacking_model(
            artifacts,
            dir_name,
            columns_to_process,
            exclude_columns,
            X_test,
            n,
        )

        df_long = pd.DataFrame(
            {
                "ID": X_test_meta["ID"],
                f"Prediction_score_fold{n}": y_test_pred_prob_meta,
                f"Prediction_Label_fold{n}": y_test_pred_meta,
            }
        )
        df_long[f"Prediction_Label_fold{n}"] = df_long[f"Prediction_Label_fold{n}"].replace(
            {0: "Negative", 1: "Positive"}
        )

        all_dfs.append(df_long)

    return all_dfs
