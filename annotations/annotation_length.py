# This script takes the annotations dataframe
# and returns the length in the desired level

import pandas as pd
from pathlib import Path
import os

def get_length(df_annotations: pd.DataFrame, level: str) -> None:
    # Get the length of the annotations in the desired level
    # Length: number of genes
    df_length = (df_annotations[[level, 'strict_locus', 'start', 'stop']]
                .groupby(level, as_index=False)
                .agg(
                    strict_locus = ('strict_locus', 'count'),
                    start = ('start', 'min'),
                    stop = ('stop', 'max')
                ))
    df_length.rename(columns = {'strict_locus': 'length'}, inplace = True)

    df_length.to_feather('../Data/Intermediate/' + level + '_length.feather')

# Update both length files, aggregating them by the species and replicon_accession
def both_length(df_annotations: pd.DataFrame) -> None:

    get_length(df_annotations, 'species')
    get_length(df_annotations, 'replicon_accession')


# Read the files with the annotations length
def read_length(level: str) -> pd.DataFrame:
    # Read the files with the annotations length
    return pd.read_feather('../Data/Intermediate/' + level + '_length.feather')

# Take only the annotations with more than threshold genes
def remove_short_replicons(df_annotations: pd.DataFrame, df_replicon_accession_length: pd.DataFrame, threshold: int) -> pd.DataFrame:
    # Remove the replicons that are less than 10 genes from df_annotations
    large_replicons = df_replicon_accession_length[df_replicon_accession_length['length'] > threshold]['replicon_accession']
    df_annotations = df_annotations[df_annotations['replicon_accession'].isin(large_replicons)]
    
    return df_annotations