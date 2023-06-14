# This script will update all the files depending on:
# - annotations (after feature_enriched step)
# - k-mers (after k_mers step)

# It will update:
# - annotations file
# - length files
# - k-mers feather file
# - k-mers count files

import pandas as pd
import logging
import argparse

import annotations.annotation_length as al
import k_mers.read_k_mers as rkm
import k_mers.count_k_mers as ckm
import annotations.read_annotations as ra
# from test import Path
from pathlib import Path

def update_all_files(input_folder, output_folder) -> None:
    # Update annotations file
    input_path_folder = Path(input_folder + 'interim/feature_enriched/')
    df_annotations = ra.read_annotations(input_path_folder)
    df_annotations.reset_index(inplace=True)
    df_annotations.to_feather(output_folder + 'df_annotations.feather')
    logging.info('Annotations file updated')

    # Update length file
    al.both_length(df_annotations)
    logging.info('Length files updated')

    # Update k-mers feather file
    rkm.update_k_mers(input_folder + 'processed/k_mers/', mode='strict')
    rkm.update_k_mers(input_folder + 'processed/k_mers/', mode='relaxed')
    logging.info('K-mers files updated')

    # Now, read the k-mers, as they are now updated
    df_k_mers_relaxed = pd.read_feather(output_folder + 'df_k_mers_relaxed.feather')
    df_k_mers_strict = pd.read_feather(output_folder + 'df_k_mers_strict.feather')

    # Update k-mers counter feather file
    ckm.count_k_mers(df_k_mers_relaxed, output_folder=output_folder + 'relaxed_')
    ckm.count_k_mers(df_k_mers_strict, output_folder=output_folder + 'strict_')
    logging.info('K-mers counter files updated')
    logging.info('All files updated')

if __name__ == '__main__':
    # Get the arguments
    parser = argparse.ArgumentParser()
    parser.add_argument('--input', type=str, default='../Data/Intermediate/Dmitrii/', help='Input file path')
    parser.add_argument('--log', type=str, default='INFO', help='Log level') # Not using it actually
    parser.add_argument('--output', type=str, default='../Data/Intermediate/accumulated/', help='Output file path')
    args = parser.parse_args()

    logging.basicConfig(
        format='%(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.INFO,
    )

    update_all_files(input_folder=args.input, output_folder=args.output)