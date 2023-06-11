
import pandas as pd
from pathlib import Path
import os


# Read all the annotations feather files of one folder and store them into a dataframe  
def read_annotations(input_file_path: Path) -> pd.DataFrame:
    # Read all the feather files of the specified folder and store them into a dataframe
    df_annotations = pd.DataFrame()
    for file in os.listdir(input_file_path):
        data = pd.read_feather(os.path.join(input_file_path, file))
        data['species_gcf_id'] = file.split('.')[0]
        
        df_annotations = pd.concat([df_annotations, data])
    
    return df_annotations


# Create the ouput folder if it does not exist
def create_output_folder(output_path = '../Data/Intermediate/alignment/') -> Path:
    # Create the output folder
    timestamp = dt.now().strftime('%Y%m%d_%H%M%S')
    output_path = output_path + timestamp

    # Create the output folder
    os.mkdir(output_path)
    logging.info(f'Output folder created: {output_path}')

    # Create the output folder for the plots
    os.mkdir(output_path + 'plots')
    logging.info(f'Output folder for the plots created: {output_path + "plots"}')

    return output_path