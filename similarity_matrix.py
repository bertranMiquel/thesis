# Aim: get a similarity matrix from all the species in the database

import pandas as pd
import logging
import argparse

import k_mers.k_mers_metric as kmm
import annotations.annotation_length as al



def similarity_matrix(mode = 'relaxed', output_file_path = '../Data/Intermediate/metric_matrix_') -> None:

    k_mers_counter = pd.read_feather('../Data/Intermediate/accumulated/' + mode + '_count_k_mers_species.feather')

    df_k_mers = pd.read_feather('../Data/Intermediate/accumulated/df_k_mers_' + mode + '.feather')

    df_length = al.read_length('species')

    # Compute the similarity measure for each species versus all the others
    # Do this using the defined metric
    measure = pd.DataFrame()
    for species in k_mers_counter['species'].unique():
        similarity = kmm.k_mers_to_metric(origin=species, level='species', k_mers_counter=k_mers_counter, df_length=df_length, df_k_mers=df_k_mers)
        similarity.rename(columns={'species': 'species2'}, inplace=True)
        similarity['species1'] = species

        measure = pd.concat([measure, similarity])

    # Convert the dataframe to a matrix
    matrix = measure.pivot(index='species1', columns='species2', values='metric')

    #  Save the metric matrix
    matrix.to_csv(output_file_path + '/metric_matrix_' + mode + '.csv', index=False)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='Create a similarity matrix from all the species in the database')
    parser.add_argument('-m', '--mode', type=str, default='strict', help='Mode to create the similarity matrix. Options: strict, relaxed')
    parser.add_argument('-o', '--output', type=str, default='../Data/Intermediate/metric/', help='Output file path')
    args = parser.parse_args()
    
    logging.basicConfig(
        format='%(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.INFO,
    )

    similarity_matrix(mode=args.mode, output_file_path=args.output)