import pandas as pd
from pathlib import Path
import os
import sys
import species.map_species_name as msn

# Read all the annotations feather files of one folder and store them into a dataframe  
def read_annotations(input_file_path: Path, species=[]) -> pd.DataFrame:

    # Read all the feather files of the specified folder and store them into a dataframe
    df_annotations = pd.DataFrame()

    # If there are species specified, read only those
    if not species:
        for file in os.listdir(input_file_path):
            data = pd.read_feather(os.path.join(input_file_path, file))            
            # Add the species GCF id to the dataframe
            data['species_gcf_id'] = file.split('.')[0]            
            df_annotations = pd.concat([df_annotations, data])
    
    elif species:
        # Map the species name to the filename with the help of the config file
        dict_species = msn.map_species_to_file()
        species_filename = [k for k, v in dict_species.items() if v in species]
        for file in os.listdir(input_file_path):
            # Remove '.features.feather' part of the string and check if it is in the list of species
            if file.split('.features')[0] in species_filename:
                data = pd.read_feather(os.path.join(input_file_path, file))
                # Add the species GCF id to the dataframe
                data['species_gcf_id'] = file.split('.')[0]
                df_annotations = pd.concat([df_annotations, data])


    return df_annotations