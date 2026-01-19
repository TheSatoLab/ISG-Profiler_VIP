# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Hiroaki Unno & Jumpei Ito

from pandas import DataFrame


def cal_z(df: DataFrame, col_name: str):
    """
    convert `col_name` z-score normalization

    :param df: modified dataframe
    :param col_name: target column for convert z-score
    """
    df[col_name] = (df[col_name] - df[col_name].mean()) / df[col_name].std()
    return df
