# Phylogenetic tree
# The aim of this script is to create a phylogenic tree based on top of the measure results
# The tree is created using the biopython library and ete3. 

# Importing libraries
import pandas as pd
import numpy as np
from Bio import Phylo
from Bio.Phylo.TreeConstruction import DistanceMatrix, DistanceTreeConstructor
from scipy.spatial.distance import squareform
from ete3 import Tree, TreeStyle, faces
from datetime import datetime as dt
import logging
from pathlib import Path
import argparse

# The aim of the function is to detect if a column has more than 1 nan-value, then:
# 1. Get the index
# 2. Get the column name
# 3. Get the row name
# 4. Remove the columns
# 5. Remove the rows
def remove_nan(matrix: pd.DataFrame) -> pd.DataFrame:
    # Get the columns with more than 1 nan
    columns = matrix.columns[matrix.isna().sum() > 1]
    # Get the rows with more than 1 nan
    rows = matrix.index[matrix.isna().sum(axis=1) > 1]

    if not columns.empty:
        logging.info('The removed species are: ' + str(columns))

    # Remove the columns
    matrix = matrix.drop(columns=columns)
    # Remove the rows
    matrix = matrix.drop(index=rows)
    return matrix

# Transform the similarity measure to a distance
def measure_to_dist(z, a):
    return np.exp(1/z**a) - np.exp(1)

# Read the similarity matrix. Return:
# - species: list of species
# - lower_triangular_list_of_lists: list of lists with the lower triangular matrix transformed to a distance
def read_measure(input_folder_path, mode, scale= 1/4) -> tuple:

    # Read the similarity matrix

    # Convert mode into a string if it is not so
    if not isinstance(mode, str):
        mode = str(mode)
    
    input_file_path = input_folder_path + '/metric_matrix_' + mode + '.csv'

    # Read the file of the input path as a csv
    matrix = pd.read_csv(input_file_path)

    matrix = remove_nan(matrix)

    species = matrix.columns.to_list() 
    matrix = matrix.to_numpy()

    # Convert similarity metric to distance with the scale parameter
    lower_triangular_list_of_lists = [[measure_to_dist(matrix[i][j], scale) if matrix[i][j] > 0 else matrix[i][j] for j in range(i+1)] for i in range(len(matrix))]

    return species, lower_triangular_list_of_lists

def create_tree(input_file_path, method: str, mode:str, scale=1/4) -> Phylo.BaseTree.Tree:
    # Read the similarity matrix
    species, lower_triangular_list_of_lists = read_measure(input_file_path, mode, scale)

    # Create a distance matrix object from your similarity matrix
    matrix = DistanceMatrix(species, lower_triangular_list_of_lists)

    # Create a DistanceTreeConstructor object
    constructor = DistanceTreeConstructor()

    # Build the tree using the neighbor joining method
    if method == 'nj':
        return constructor.nj(matrix)
    elif method == 'upgma':
        return constructor.upgma(matrix)
    else:
        raise ValueError('Method not recognized')
    

# Create and save the phylogenetic tree with ete3 library
def philo_tree(tree, output_file_path, mode, method) -> None:

    t = Tree(tree.format("newick"), format=1, quoted_node_names=True) #tree_string

    # Different customizations options forr the tree
    '''
    # Draws a circular tree using a semi-circumference: 
    ts = TreeStyle()
    ts.mode = "c"
    ts.arc_start = -180 # 0 degrees = 3 o'clock
    ts.arc_span = 180

    ts = TreeStyle()
    ts.rotation = 90

    # Add title
    title_face = faces.TextFace("My Phylogenetic Tree")
    t.add_face(title_face, column=0, position="branch-top")

    # Add a description to the entire tree in the bottom-left corner
    desc_face = faces.TextFace("This is a description of the entire tree")
    desc_face.margin_top = 5
    desc_face.margin_left = 5
    desc_face.margin_right = 5
    desc_face.border.width = 1
    desc_face.border.color = "gray"

    t.add_face(desc_face, column=0, position="branch-bottom")

    # Add a blank face to the bottom-right corner to align the description
    # with the bottom-left corner
    blank_face = faces.TextFace("")
    blank_face.margin_top = 5
    blank_face.margin_right = 5

    t.add_face(blank_face, column=1, position="branch-bottom")
    '''
    t.show() #tree_style=ts

    # Save the tree into svg file
    timestamp = dt.now().strftime("%Y%m%d-%H%M%S")
    t.render(output_file_path + timestamp +'_' + mode + '_' + method + "_tree.svg")


if __name__ == '__main__':
    # Parse arguments
    parser = argparse.ArgumentParser(description='Create a phylogenetic tree based on the similarity measure')
    parser.add_argument('-i', '--input', type=str, help='Input folder path', required=True, default='../Data/Intermediate/metric/')
    parser.add_argument('-o', '--output', type=str, help='Output file path', required=True, default='../Data/Results/tree/')
    parser.add_argument('-m', '--method', type=str, help='Method to create the tree', required=True, choices=['nj', 'upgma'])
    parser.add_argument('-mode', '--mode', type=str, help='Mode to create the tree', required=True, default='strict')
    args = parser.parse_args()

    logging.basicConfig(
        format='%(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.INFO,
    )

    tree = create_tree(input_file_path=args.input, mode=args.mode, method=args.method)
    logging.info('Tree created')
    philo_tree(tree=tree, output_file_path=args.output, mode=args.mode, method=args.method)
    logging.info('Tree saved')
    