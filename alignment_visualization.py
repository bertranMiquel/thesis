# In this script the alignment visualization is implemented
# The alignment visualization is a 2D plot of the alignment
# We read the output file of the alignment called blocks.csv under the specific folder of the alignment
# Each block means a part of the genome aligned with another part of the genome of a different species

# Import packages
from typing import List
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as n
from pathlib import Path
import logging
import argparse

# Read the output file of the alignment
# Give format to some columns
def read_blocks(input_file_path) -> pd.DataFrame:
    blocks = pd.read_csv(input_file_path)

    # Redefine the start and stop columns as integers
    blocks.start = blocks.start.astype('int64')
    blocks.stop = blocks.stop.astype('int64')

    # Set match_perc as percentage format with 2 decimals
    blocks.match_perc = blocks.match_perc.astype('float64')
    blocks.match_perc = blocks.match_perc.apply(lambda x: '{:.2%}'.format(x))

    return blocks

# Merge the length file with the blocks file to obtain the longitudes of the scaffolds
def merge_length(input_file, blocks: pd.DataFrame) -> pd.DataFrame:
    # Read the length file
    length = pd.read_feather(input_file)

    # Merge the length file with the blocks file
    blocks = blocks.merge(length, how='left', on='replicon_accession', suffixes=('', '_replicon_accession'))
    # drop the length column of the replicon_accession, since represents the number of genes, might led to confussion
    blocks = blocks.drop(columns=['length'])

    return blocks

# Create ids for the species and the scaffolds
def assign_id(order: List, blocks: pd.DataFrame) -> None:
    # Create the id for the target scaffold
    blocks['target_replicon_accession_id'] = blocks.sort_values('alignment_id').reset_index().groupby('target_replicon_accession').ngroup()

    # Create the id for the species
    # Sort blocks table by the desired order of the species
    # Create a new column with the desired order of the species
    if not order:
        blocks['species_ordered'] = pd.Categorical(blocks['species'])
    else:    
        # First let us flat the list of lists
        order = [item for sublist in order for item in sublist]
        order = list(dict.fromkeys(order))
        blocks['species_ordered'] = pd.Categorical(blocks['species'], categories=order, ordered=True)

    # get the category codes (i.e. IDs) for each entry in the ordered column
    blocks['species_id'] = blocks['species_ordered'].cat.codes
    blocks['species_id'] = blocks['species_id'].astype('int64')

# Recompute the start and stop position of each scaffold to order them one after the other
# So, in the alignment visualization they will be represented one after the other and not superimposed
def recalculate_start_stop(blocks: pd.DataFrame) -> None:
    # Sort values for proper iteration
    blocks.sort_values(['species', 'replicon_accession'], inplace=True) # or alignment_id, depend on order priority
    blocks.reset_index(inplace=True, drop=True)

    # A new variable is needed to keep the original stop value and 
    # update the stop and start values by the same amount
    prev_stop = 0

    for i in range(1,len(blocks)):

        if blocks.loc[i, 'species'] == blocks.loc[i-1, 'species']:
            if blocks.loc[i, 'replicon_accession'] != blocks.loc[i-1, 'replicon_accession']:
                prev_stop =  blocks.loc[i-1, 'stop_replicon_accession']

            blocks.loc[i, 'stop'] += prev_stop
            blocks.loc[i, 'start'] += prev_stop
            blocks.loc[i, 'start_replicon_accession'] = prev_stop
            blocks.loc[i, 'stop_replicon_accession'] += prev_stop

    # Substitute nan values to 0
    blocks.start_replicon_accession.fillna(0, inplace=True)


# Create a color mapping per each alignment_id
# Create a dictionary where the key is the alignment_id and the value is the color
# The color should be a tuple of RGB values
# I want to create a color mapping that is different for each target_replicon_accession
# If a target_replicon_accession has multiple alignment_ids, then the color mapping should be the same for all of them
# We need to map each target_replicon_accession to a integer value
# For doing so, first I choose a color map and then I normalize the values of the target_replicon_accession
def color_mapping(blocks: pd.DataFrame) -> dict:
    
    # source for different color mappings: https://matplotlib.org/stable/tutorials/colors/colormaps.html#qualitative
    cmap = mpl.cm.tab20 # mpl.cm.viridis
    norm = mpl.colors.Normalize(vmin=blocks['alignment_id'].min(), vmax=blocks['alignment_id'].max())

    # Create a dictionary with the color mapping
    return dict(zip(blocks.alignment_id.unique(), cmap(norm(blocks['alignment_id']))))

# We would like to represent the scaffold information to directly plot its information. 
# Therefore, we will create a new dictionary with the scaffold information as entry
# containing the start and end of the scaffold, and the species_id.
def scaffold_data(blocks: pd.DataFrame) -> dict:
    # Create a dictionary with the scaffold information
    # key: replicon_accession
    # value: tuple of start and stop position of the scaffold
    # and the species_id of the same row as the key
    
    # Create a unique dataframe with the scaffold name, species and the maximum start and stop values
    scaffold_info = (blocks
                        .groupby(['replicon_accession', 'species_id'], as_index=False)
                        [['start_replicon_accession', 'stop_replicon_accession']]
                        .agg(
                            {'start_replicon_accession': 'min', 
                            'stop_replicon_accession': 'max'}
                        )
    )

    # Convert it to a dictionary
    scaffold_info_dict = scaffold_info.set_index('replicon_accession').T.to_dict('list')

    return scaffold_info_dict

# Plot the alignment
def plot_alignment(blocks: pd.DataFrame, color_mapping: dict, replicon_accession_length: dict, output_file_path: str, timestamp: str) -> None:
    
    # Create a dictionary with the blocks information it is easy to handle the information
    # Create a dictionary where the key is the alignment_id and the value is the array of tuples, including the species_id
    blocks.sort_values(['alignment_id', 'species_id'], ascending=True, inplace=True)
    blocks_dict = (blocks
                    .groupby('alignment_id')[
                                                ['start_length', 
                                                'species_id', 
                                                'replicon_accession', 
                                                'match_perc', 
                                                'stop_replicon_accession']
                                            ]
                    .apply(lambda x: x.values.tolist())
                    .to_dict())

    # Start plotting
    fig, ax = plt.subplots()

    # Define plot size
    fig.set_size_inches(30, 10)

    # Each key is an alignment_id, which determines the color of the bar
    # Each block is represented per species, so we need a y-level per each species
    # The y-level is determined by the species_id
    # The x-level is determined by the start and length of the block
    for key, value in blocks_dict.items():
            for i in value:
                    # Draw the block
                    ax.broken_barh([i[0]],  # start and length x-axis
                                    (i[1], 0.5),            # start and height y-axis
                                    facecolors=color_mapping[key], # = target_replicon_accession; key, if we want per alignment
                                    zorder=2)


            # Add a line between the blocks for each alignment_id
            # The start positions of each block are connected by a line
            ax.plot([value[1][0][0], value[0][0][0]],
                    [value[1][1], value[0][1] + 0.5],
                    color=color_mapping[key], 
                    linestyle='-', 
                    linewidth=1)

            # The end (= start + length) positions of each block are connected by a line
            ax.plot([value[1][0][0] + value[1][0][1], value[0][0][0] + value[0][0][1]],
                    [1, 0.5], #[i[1] for i in value], 
                    color=color_mapping[key], 
                    linestyle='-', 
                    linewidth=1)
            
            # Set match_perc KPI between the two lines
            # First, compute the middle of the two lines for the x-axis
            mid_bottom_block_pos = (value[0][0][0] + (value[0][0][0] + value[0][0][1])) / 2 
            mid_top_block_pos = (value[1][0][0] + (value[1][0][0] + value[1][0][1])) / 2
            
            in_between = (mid_top_block_pos + mid_bottom_block_pos) / 2

            # Insert the KPI
            ax.text(x=in_between, # + (value[1][0][0] + (value[1][0][0] + value[1][0][1]) / 2)) / 2,
                    y=0.75,
                    s=value[1][3],
                    ha='center',
                    va='center',
                    fontsize=10,
                    weight='bold',
                    rotation='horizontal')

    ax.set_ylim(0, blocks.species.nunique())
    ax.set_xlim(min(blocks.start) - 1000, max(blocks.stop) + 1000)
    #ax.set_xlabel('Genomic position')

    # Make grid lines visible
    ax.grid(True, axis='y', linestyle='-', color='black')

    # Hide axes
    for key, spine in ax.spines.items():
        spine.set_visible(False)

    # Set tick labels
    ax.set_yticks([i + 0.25 for i in blocks.species_id.unique()], labels=blocks.species.unique())
    ax.set_xticks([])

    # With the help of scaffold_data, we can draw a vertical line for each scaffold and its name
    for key, value in replicon_accession_length.items():
            # Draw the line
            plt.vlines(x=value[2], ymin=value[0], ymax=value[0]+0.5, color='black', linestyle='-', linewidth=1)
            
            # Write the name of the scaffold in each aligned block between a box
            # Create the box object
            box = dict(boxstyle='round', facecolor='white', alpha=1)

            ax.text(x=(value[1] + value[2])/2,
                    y=value[0]+0.25,
                    s=key,
                    ha='center',
                    va='center',
                    rotation='horizontal',
                    bbox=box)
            
    # Set the x-axis and y-axis scales to 'auto'
    #ax.autoscale(enable=True, axis='both', tight=True)

    # Add a title
    ax.set_title('Alignment of Gorilla vs. Human specifying % /of matching genes in aligned blocks', fontsize=20)

    # Save the plot
    fig.savefig(output_file_path + '/alignment_viz_' + timestamp +'.svg', format='svg', dpi=1200)


def alignment_vis(path: str, timestamp: str, order, output_file_path: str) -> None:

    # Create the input folder path
    folder_path = path + '/' + timestamp

    # Read the blocks
    blocks = read_blocks(input_file_path=folder_path + '/blocks.csv')
    logging.info('Blocks read')

    # Keep only the ones we are interested in
    # So, the pairs of comparing_species that appear in order
    if order:
        blocks = blocks[blocks['comparing_species'].isin(order)]
        logging.info('Blocks filtered')
    else:
        logging.info('No order specified, all blocks are kept')

    # Add length information
    blocks = merge_length('../Data/Intermediate/replicon_accession_length.feather', blocks=blocks)

    assign_id(blocks=blocks, order=order)

    recalculate_start_stop(blocks)

    # Create length column
    blocks['length'] = blocks['stop'] - blocks['start']

    # Create a tuple of the start and length of each block
    blocks['start_length'] = blocks.apply(lambda x: (x.start, x.length), axis=1)

    # Create a dictionary where the key is the alignment_id and the value is the color
    color_mapping_dict = color_mapping(blocks)

    # Create a dictionary where the key is the scaffold and the value is the length of the scaffold
    replicon_accession_length = scaffold_data(blocks)

    # Create and save the plot
    logging.info('Creating plot')
    plot_alignment(blocks, color_mapping_dict, replicon_accession_length, output_file_path, timestamp)


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description='Visualize the alignment of two species')
    parser.add_argument('-p', '--path', type=str, default='../Data/Intermediate/alignment/', help='Path to the folder containing the alignment folders')
    parser.add_argument('--timestamp', type=str, help='Timestamp of the folder containing the blocks.csv file')
    parser.add_argument('--order', nargs='+', default=[], help='Order of the species to be visualized')
    parser.add_argument('-o', '--output_file_path', type=str, default='../Data/Results/alignment_visualization', help='Path to the output')
    args = parser.parse_args()

    logging.basicConfig(
        format='%(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.INFO,
    )

    # Run the function
    alignment_vis(path=args.path, timestamp=args.timestamp, order=args.order, output_file_path=args.output_file_path)
    