# The aim of this script is to read k-mers from all generated files
# We can find those in the folder Data/Intermediate/k_mers/processed/
# where we have one folder per different species

import os
import pandas as pd
from sympy import isprime
import yaml
from pathlib import Path
import species.map_species_name as msn

# Define a function to read a single feather file and return a dictionary
def read_file(filepath, mode = 'strict'):
    if mode in filepath:
        # Read the feather file
        df = pd.read_feather(filepath)

        folder, filename = os.path.split(filepath)
        
        # Get the species name from the filename
        species_gcf_id = folder.split('/')[-1]

        # Get the k-mer length from the filename
        k = filename.split('.')[0]

        # Create a dictionary with the species name, the k-mer length and the dataframe
        data_dict = {'species_gcf_id': species_gcf_id, 'k': k, 'df': df}

        return data_dict
    else:
        return None


def process_folder(path: Path, mode: str):
    # Read all the feather files and store the data in a list of dictionaries
    data_list = []
    for folder in os.scandir(path):
        # There are not only folders, but also other files in this directory
        # Checck it is a folder
        if folder.is_dir():
            folder_path = os.path.join(path, folder.name)
            for filename in os.listdir(folder_path):
                filepath = os.path.join(folder_path, filename)
                k = filename.split('.')[0]
                # Do not read sorted files
                if 'sorted' in filepath:
                    continue
                # Check if it is a number and if it is prime
                if k.isdigit():
                    if isprime(int(k)):
                        data_dict = read_file(filepath, mode=mode)
                        if data_dict is not None:
                            data_list.append(data_dict)
    return data_list


def update_k_mers(input_file_path: Path, mode: str) -> None:
    # Read all the feather files of the specified folder and store them into a list of dictionaries
    data_list = process_folder(input_file_path, mode=mode)

    # Convert the list of dictionaries into a single DataFrame, including the species and k columns
    df_list = pd.concat([d['df'].assign(species_gcf_id=d['species_gcf_id'], k=d['k']) for d in data_list], ignore_index=True)
    df_list = df_list[['species_gcf_id', 'k', 'k_mer', 'replicon_accession']]           #.rename(columns={'replicon_accession': 'scaffold'})

    # Create a column in the dataframe with the species name
    file_species = msn.map_species_to_file()
    df_list['species'] = df_list['species_gcf_id'].map(file_species)

    # Save the dataframe as a feather file
    df_list.to_feather('../Data/Intermediate/accumulated/df_k_mers_' + mode + '.feather')
