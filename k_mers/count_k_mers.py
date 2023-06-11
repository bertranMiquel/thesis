# This script contains a function to compute the number of k-mers per scaffold or species
# It will create two feather files:
# - count_k_mers_scaffold.feather: containing the number of k-mers per scaffold
# - count_k_mers_species.feather: containing the number of k-mers per species
# This information is useful to compute the similarity metric between scaffolds or species in k_mers_metric.py

import pandas as pd

def count_k_mers(df_k_mers: pd.DataFrame, output_folder='../Data/Intermediate/k_mers/') -> None:
    # From df_k_mers, create a dataframe with the number of k_mers per scaffold / species
    df_k_mers_scaffold = df_k_mers.groupby(['species', 'replicon_accession', 'k'], as_index=False).count()
    df_k_mers_species = df_k_mers_scaffold.groupby(['species', 'k'], as_index=False).sum()

    df_k_mers_scaffold = df_k_mers_scaffold[['species', 'replicon_accession', 'k', 'k_mer']].rename(columns = {'k_mer': 'total_k_mers'})
    df_k_mers_species = df_k_mers_species[['species', 'k', 'k_mer']].rename(columns = {'k_mer': 'total_k_mers'})

    # Save them into a feather file
    df_k_mers_scaffold.to_feather(output_folder + 'count_k_mers_replicon_accession.feather')
    df_k_mers_species.to_feather(output_folder + 'count_k_mers_species.feather')