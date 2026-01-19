# SPDX-License-Identifier: GPL-3.0-only
# SPDX-FileCopyrightText: Copyright 2026 Luca Nishimura & Jumpei Ito

from pandas import DataFrame


def calculate_isg_scores(
    all_sample_gene_count_df: DataFrame, sample_metadata: DataFrame
) -> DataFrame:
    """
    Calculates ISG scores for each sample and merges them with host metadata.

    The ISG score is calculated as the mean of the 'standardized_count' for genes
    labeled as 'ISG' in the 'type' column. This score is then merged with
    taxonomic information from the metadata.

    :param all_sample_gene_count_df: DataFrame containing gene counts and types for all samples.
        The DataFrame is expected to have columns such as 'sample_id', 'hum_symbol',
        'raw_count', 'type', 'normalized_count', and 'standardized_count'.
        **Required columns:** 'sample_id', 'type', 'standardized_count'.
    :type all_sample_gene_count_df: DataFrame
    :param sample_metadata: DataFrame containing host metadata for the samples.
        **Required columns:** 'sample_id', 'species_host', 'order_host', 'clade_host'.
    :type sample_metadata: DataFrame
    :return: A DataFrame containing the calculated ISG scores and host taxonomy.
        The resulting columns are:
        - sample_id (str)
        - ISG_score (float): Mean standardized count of ISG genes.
        - species_host (str)
        - order_host (str)
        - clade_host (str)
    :rtype: DataFrame

    Examples:
        >>> # Example input structure based on provided data
        >>> df_counts = pd.DataFrame({
        ...     'sample_id': ['ERR12917750', 'ERR12917750', 'ERR2012446', 'ERR2012446'],
        ...     'hum_symbol': ['ACTR2', 'ADAR', 'ZCCHC2', 'ZNFX1'],
        ...     'type': ['cntl', 'ISG', 'ISG', 'ISG'],
        ...     'standardized_count': [-1.76, -1.81, -0.10, 0.45]
        ... })
        >>> df_meta = pd.DataFrame({
        ...     'sample_id': ['ERR12917750', 'ERR2012446'],
        ...     'species_host': ['Gallus_gallus', 'Gallus_gallus'],
        ...     'order_host': ['Galliformes', 'Galliformes'],
        ...     'clade_host': ['Aves', 'Aves']
        ... })
        >>> calculate_isg_scores(df_counts, df_meta)
             sample_id  ISG_score   species_host   order_host clade_host
        0  ERR12917750      -1.81  Gallus_gallus  Galliformes       Aves
        1   ERR2012446       0.175 Gallus_gallus  Galliformes       Aves
    """
    isg_score_df = (
        all_sample_gene_count_df.loc[all_sample_gene_count_df["type"] == "ISG"]
        .groupby("sample_id", as_index=False)
        .agg(ISG_score=("standardized_count", "mean"))
    )

    isg_score_df = isg_score_df.merge(
        sample_metadata[["sample_id", "species_host", "order_host", "clade_host"]],
        on="sample_id",
        how="left",
    )

    return isg_score_df
