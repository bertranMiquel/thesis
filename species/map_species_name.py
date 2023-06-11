# Function to map those folder names ( = gcf_id) to the species name

import os
import yaml

def map_species_to_file(input_path='./k_mers_generation/gff_to_k_mer_pipeline/config/config.yaml'):
    # Open config.yaml file
    with open(input_path, 'r') as f:
        config = yaml.load(f, Loader=yaml.FullLoader)

    # Create a dictionary with the filename and the species name from the config file
    # First, take only the organisms with a filename
    organisms = {k: v for k, v in config['organisms'].items() if v['filename'] is not None}

    # Create the dictionary with the filename as key and the species as value
    file_species = {v['filename']: v['name'] for k, v in organisms.items()}

    return file_species