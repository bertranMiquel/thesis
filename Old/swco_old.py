'''
All functions defined for applying the Smith Watermann algorithm & other useful functions
'''
import pandas as pd
import numpy as np
from enum import IntEnum
import os

# Function to read all the excel tabs of an excel file 
# where each tab contains the species genome
def read_all(filepath):

    # Load Excel file using Pandas with `sheet_name=None`
    df_dict = pd.read_excel(filepath, sheet_name=None)

    df_species = pd.DataFrame()

    species = df_dict.keys()

    # Flag each specie
    for s in species:
        aux = df_dict.get(s)
        aux['Specie'] = s
        df_species = pd.concat([df_species, aux])

    return preprocessing(df_species)

# Data cleaning and creation of new fields
# - Gene_non_or: Gene name without number part
# - Gene: Gene name with orientation
# - Specie_Scaffold: Specie + Scaffold
# Clean the genes that are not LOC or - or empty
# Reset index
def preprocessing(specie):
    
    specie['Gene_non_or'] = specie['Locus'].str.split('(\d+)').str[0].str.upper()
    specie['Gene'] = specie['Gene_non_or'] + specie['Strand']
    specie['Specie_Scaffold'] = specie['Specie'] + '_' + specie['Replicon Accession']

    specie = specie[(specie['Gene_non_or'].str.contains('LOC') == False) & (specie['Gene_non_or'].str.isspace() == False) & (specie['Gene_non_or'] != '-') ]  
    
    specie.reset_index(inplace= True)
    specie.drop(columns=['index'], inplace=True)

    return specie


# Clean consecutive duplicates
# Consecutive duplicates are those that have:
# - Same Replicon Name
# - Same Accession
# - Same Gene
# - Same orientation
# - Start prev <= Start < End prev
def cons_duplicates_pre(data):
    data.sort_values(by = [data.columns[1], data.columns[2]], inplace=True)
    # Same Replicon Name, same Accession, same Gene, same orientation, Start prev <= Start < End prev
    data['New'] = np.where((data['#Replicon Name'].shift() == data['#Replicon Name']) & (data['Replicon Accession'].shift() == data['Replicon Accession']) & (data['Gene'].shift() == data['Gene']) & (data['Start'].shift() <= data['Start']) & (data['Start'] < data.iloc['Stop'].shift()), 0, 1)
    data['Change'] = data['New'].cumsum()
    return data.groupby('Change', as_index=False).first()


# Remove consecutive duplicates, here are considered as:
# - Same Replicon Name
# - Same Accession
# - Same Gene
# For creating the kmers
def cons_duplicates_kmers(data):
    data.sort_values(by = ['Specie', '#Replicon Name', 'Replicon Accession', 'Start'], inplace=True)
    data['New'] = np.where((data['#Replicon Name'].shift() == data['#Replicon Name']) & (data.iloc['Replicon Accession'].shift() == data.iloc['Replicon Accession']) & (data['Gene_non_or'].shift() == data['Gene_non_or']), 0, 1)
    data['Change'] = data['New'].cumsum()
    return data.groupby('Change', as_index=False).first().rename(columns={'Change':'Index'})


def match_scaffold(k, specie, scaffold, min, df_k_mers):
    
    sp_sca = specie + '_' + scaffold

    # Data clean
    # Filter
    df_k_mers = df_k_mers.loc[(df_k_mers['k'] == k), ['k_mers', 'Specie_Scaffold']]
    df_query = df_k_mers.loc[(df_k_mers['Specie_Scaffold'].str.contains(sp_sca, case=False)), ['k_mers', 'Specie_Scaffold']]

    df_query['Specie_Scaffold'] = df_query['Specie_Scaffold'].str.replace("\['", "", regex=True)
    df_query['Specie_Scaffold'] = df_query['Specie_Scaffold'].str.replace("'\]", "", regex=True)
    
    # List all scaffolds per k_mer
    sep = df_query['Specie_Scaffold'].str.split("', '", expand=True)

    df_query = sep.merge(df_query['k_mers'], left_index=True, right_index=True, how='right')

    # Make table: k_mer | Specie_Scaffold
    melt = pd.melt(df_query, id_vars=['k_mers']).dropna().drop(labels='variable', axis=1).rename(columns={'value': 'Specie_Scaffold'})

    match = (melt.groupby('Specie_Scaffold', as_index=False)
                .count()
                .rename(columns={'k_mers':'# Appearances'})
                .sort_values('# Appearances', ascending=False))
    
    match[['Specie', 'Scaffold']] = match['Specie_Scaffold'].astype(str).str.split("_", n = 1, expand=True)

    # Apply minimum threshold & do not return same specie apperances
    match = match.loc[(match['# Appearances'] > min) & (match['Specie'].str.contains(specie, case=False, regex=True)==False)]

    return match  


def match_scaffold_new(k, specie, scaffold, min, df_k_mers):

    # Take only the k_mers of interest
    df_k_mers = df_k_mers.loc[(df_k_mers['k'] == k), ['k_mers', 'Specie_Scaffold']]

    # filter dataframe to only include rows where is the desired specie and scaffold
    sp_sca = specie + '_' + scaffold
    df_filtered = df_k_mers.loc[(df_k_mers['Specie_Scaffold'].str.contains(sp_sca, case=False)), ['k_mers', 'Specie_Scaffold']]

    # group by k_mers, count occurrences of each k_mers, and sort in descending order by count
    grouped = df_filtered.groupby('k_mers')['Specie_Scaffold'].count().reset_index(name='Count')
    sorted_grouped = grouped.sort_values(by='Count', ascending=False)

    # Keep only the categories which have more than min occurrences
    sorted_grouped = sorted_grouped.loc[sorted_grouped['Count'] > min]

    return sorted_grouped  


# Create a data frame with the length of each scaffold for all the species
# Different lengths:
# - Number of genes
# - Start of base pairs per scaffold
# - Maximum value for stop of base pairs to know how many base pairs are per scaffold
def df_scaffold_length(data):
    # Create a data frame with the length of each scaffold for all the species
    data['Start'] = data['Start'].astype('int64')
    data['Stop'] = data['Stop'].astype('int64')

    df_scaffold_length = (data[['Specie', '#Replicon Name', 'Replicon Accession', 'Gene', 'Start', 'Stop']]
                            .groupby(['Specie', '#Replicon Name', 'Replicon Accession'], as_index=False)
                            .agg({'Gene': 'count', 'Start': 'min', 'Stop': 'max'}))

    df_scaffold_length.rename(columns={'Gene':'Number of genes'}, inplace=True)
    df_scaffold_length.to_csv('../Data/Intermediate/df_scaffold_length.csv', index=False)

    return df_scaffold_length

# Return the number of genes in a specific scaffold, given the dataframe with all the scaffold lengths
def scaffold_length(specie, scaffold, df_scaffold_length):
    return df_scaffold_length.loc[(df_scaffold_length['Specie'].str.contains(specie, case=False, regex=True)) & (df_scaffold_length['Replicon Accession'] == scaffold), 'Number of genes'].values[0]

# Return the number of base pairs in a specific scaffold, given the dataframe with all the scaffold lengths
def scaffold_base_pairs(specie, scaffold, df_scaffold_length):
    return df_scaffold_length.loc[(df_scaffold_length['Specie'].str.contains(specie, case=False, regex=True)) & (df_scaffold_length['Replicon Accession'] == scaffold), 'Stop'].values[0]

def excel_reader(species):
    data = pd.read_excel('../Data/Raw/Tables_Filtered_IK.xlsx', species)
    
    data['Gene'] = data['Locus'].str.split('(\d+)').str[0] + data['Strand']
    data['Gene_non_or'] = data['Locus'].str.split('(\d+)').str[0]
    data.reset_index(inplace=True)

    data = data[data['Gene'].str.contains('LOC')==False]

    return data


"""

LIFE97011 - Computing
Python Programming - Assessed Exercise No. 3
Task: Smith Waterman local alignment
@Author: Slaviana Pavlovich

"""

# Assigning the constants for the scores
class Score(IntEnum):
    # order by value it should be
    MATCH = 10      
    CONS_GAP = -2   # Want to penalize worst a first gap than a secong gap
    FIRST_GAP = -4
    MISMATCH = -3
    #INV = 0
    #DUP = 1

# Match: 10 , First_gap: -4, cons_gap: -1, mismatch: -5

# Assigning the constant values for the traceback
class Trace(IntEnum):
    STOP = 0
    LEFT = 1 
    UP = 2
    DIAGONAL = 3


def smith_waterman(seq_1, seq_2, mode):

    # Create arrays for making the comparison faster
    seq_1_list = seq_1[mode + '_locus'].tolist()
    seq_2_list = seq_2[mode + '_locus'].tolist()

    seq_1_index = seq_1.index.tolist()
    seq_2_index = seq_2.index.tolist()

    row = len(seq_1) + 1
    col = len(seq_2) + 1
    matrix = np.zeros(shape=(row, col), dtype=int)
    tracing_matrix = np.zeros(shape=(row, col), dtype=int) 

    # Initialising the variables to find the highest scoring cell
    max_score = -1
    max_index = (-1, -1)

    diagonal_score = 0
    vertical_score = 0
    horizontal_score = 0

    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):
            
            # Calculating the diagonal score (match score)
            match_value = Score.MATCH if seq_1_list[i - 1] == seq_2_list[j - 1] else Score.MISMATCH
                
            diagonal_score = matrix[i - 1, j - 1] + match_value

            # Calculating the vertical gap score, penalizing less consecutive gaps and identifying duplicates
            if matrix[i - 1, j] == vertical_score:
                vertical_score = matrix[i - 1, j] + Score.CONS_GAP
            else: vertical_score = matrix[i - 1, j] + Score.FIRST_GAP
            
            # Calculating the vertical gap score, penalizing consecutive gaps and identifying duplicates
            if matrix[i, j - 1] == horizontal_score:
                horizontal_score = matrix[i, j - 1] + Score.CONS_GAP
            else: horizontal_score = matrix[i, j - 1] + Score.FIRST_GAP
            
            # Taking the highest score 
            score = max(0, diagonal_score, vertical_score, horizontal_score)
            matrix[i, j] = score

            # Tracking where the cell's value is coming from    
            if score == 0: 
                tracing_matrix[i, j] = Trace.STOP
                
            elif score == diagonal_score: 
                tracing_matrix[i, j] = Trace.DIAGONAL 
                
            elif score == horizontal_score: 
                tracing_matrix[i, j] = Trace.LEFT
                
            elif score == vertical_score: 
                tracing_matrix[i, j] = Trace.UP
                
            # Tracking the cell with the maximum score
            # If we want different strings, here we can define a threshold and keep track of all the higher values
            if score > max_score:
                max_index = (i,j)
                max_score = matrix[i, j]

    # Tracing back the path to find the alignment
    # Initialising the variables for tracing
    aligned_seq_1 = []
    aligned_seq_2 = []
    index_seq_1 = []
    index_seq_2 = []
    (max_i, max_j) = max_index

    # Tracing and computing the pathway with the local alignment
    while tracing_matrix[max_i, max_j] != Trace.STOP:
        if tracing_matrix[max_i, max_j] == Trace.DIAGONAL:
            aligned_seq_1.insert(0, seq_1_list[max_i - 1])
            aligned_seq_2.insert(0, seq_2_list[max_j - 1])

            index_seq_1.insert(0, seq_1_index[max_i - 1])
            index_seq_2.insert(0, seq_2_index[max_j - 1])

            max_i = max_i - 1
            max_j = max_j - 1

        elif tracing_matrix[max_i, max_j] == Trace.UP:
            aligned_seq_1.insert(0, seq_1_list[max_i - 1])
            aligned_seq_2.insert(0, '-')

            index_seq_1.insert(0, seq_1_index[max_i - 1])
            index_seq_2.insert(0, '-')

            max_i = max_i - 1    

        elif tracing_matrix[max_i, max_j] == Trace.LEFT:
            aligned_seq_1.insert(0, '-')
            aligned_seq_2.insert(0, seq_2_list[max_j - 1])

            index_seq_1.insert(0, '-')
            index_seq_2.insert(0, seq_2_index[max_j - 1])
            
            max_j = max_j - 1

    return aligned_seq_1, aligned_seq_2, index_seq_1, index_seq_2, max_score, max_index, max_i, max_j, matrix
      

# Implementing the Smith Waterman local alignment
def smith_waterman_old(seq_1, seq_2):

    seq_1 = seq_1[['index', 'Gene', 'Gene_non_or']]
    seq_2 = seq_2[['index', 'Gene', 'Gene_non_or']]

    row = len(seq_1) + 1
    col = len(seq_2) + 1
    matrix = np.zeros(shape=(row, col), dtype=int)
    tracing_matrix = np.zeros(shape=(row, col), dtype=int) 

    # Initialising the variables to find the highest scoring cell
    max_score = -1
    max_index = (-1, -1)

    diagonal_score = 0
    vertical_score = 0
    horizontal_score = 0
    
    # Calculating the scores for all cells in the matrix
    for i in range(1, row):
        for j in range(1, col):
            
            # Calculating the diagonal score (match score)
            # If there is a match, always win? And then, we can skip all further comparations?
            match_value = Score.MATCH if seq_1.iloc[i - 1, 2] == seq_2.iloc[j - 1, 2] else Score.MISMATCH
            #Score.INV if seq_1.iloc[i - 1, 2] == seq_2.iloc[j - 1, 2] else  --> No sense, it will make sense if the overall is inversed and then one of them is turned around. So, first, check the overall waz and then compare them.
                
            diagonal_score = matrix[i - 1, j - 1] + match_value

            # Calculating the vertical gap score, penalizing less consecutive gaps and identifying duplicates
            #if seq_2.iloc[j-1, 2] == seq_2.iloc[j - 2, 2]:
                #vertical_score = matrix[i - 1, j] + Score.DUP elif
            if matrix[i - 1, j] == vertical_score:
                vertical_score = matrix[i - 1, j] + Score.CONS_GAP
            else: vertical_score = matrix[i - 1, j] + Score.FIRST_GAP
            
            # Calculating the vertical gap score, penalizing consecutive gaps and identifying duplicates
            #if seq_1.iloc[i-1, 2] == seq_1.iloc[i-2, 2]:
                #horizontal_score = matrix[i, j - 1] + Score.DUP elif
            if matrix[i, j - 1] == horizontal_score:
                horizontal_score = matrix[i, j - 1] + Score.CONS_GAP
            else: horizontal_score = matrix[i, j - 1] + Score.FIRST_GAP
            
            # Taking the highest score 
            matrix[i, j] = max(0, diagonal_score, vertical_score, horizontal_score)
            
            # Tracking where the cell's value is coming from    
            if matrix[i, j] == 0: 
                tracing_matrix[i, j] = Trace.STOP
                
            elif matrix[i, j] == diagonal_score: 
                tracing_matrix[i, j] = Trace.DIAGONAL 
                
            elif matrix[i, j] == horizontal_score: 
                tracing_matrix[i, j] = Trace.LEFT
                
            elif matrix[i, j] == vertical_score: 
                tracing_matrix[i, j] = Trace.UP
                
            # Tracking the cell with the maximum score
            # If we want different strings, here we can define a threshold and keep track of all the higher values
            if matrix[i, j] >= max_score:
                max_index = (i,j)
                max_score = matrix[i, j]


    # Initialising the variables for tracing
    aligned_seq_1 = []
    aligned_seq_2 = []
    index_seq_1 = []
    index_seq_2 = []
    (max_i, max_j) = max_index

    # Tracing and computing the pathway with the local alignment
    # if max_score > 3 * Score.MATCH
    while tracing_matrix[max_i, max_j] != Trace.STOP:
        if tracing_matrix[max_i, max_j] == Trace.DIAGONAL:
            aligned_seq_1.insert(0, seq_1.iloc[max_i - 1, 2])
            aligned_seq_2.insert(0, seq_2.iloc[max_j - 1, 2])

            index_seq_1.insert(0, seq_1.iloc[max_i - 1, 0])
            index_seq_2.insert(0, seq_2.iloc[max_j - 1, 0])

            max_i = max_i - 1
            max_j = max_j - 1

        elif tracing_matrix[max_i, max_j] == Trace.UP:
            aligned_seq_1.insert(0, seq_1.iloc[max_i - 1, 2])
            aligned_seq_2.insert(0, '-')

            index_seq_1.insert(0, seq_1.iloc[max_i - 1, 0])
            index_seq_2.insert(0, '-')

            max_i = max_i - 1    

        elif tracing_matrix[max_i, max_j] == Trace.LEFT:
            aligned_seq_1.insert(0, '-')
            aligned_seq_2.insert(0, seq_2.iloc[max_j - 1, 2])

            index_seq_1.insert(0, '-')
            index_seq_2.insert(0, seq_2.iloc[max_j - 1, 0])
            
            max_j = max_j - 1

    return aligned_seq_1, aligned_seq_2, index_seq_1, index_seq_2, max_score, max_index, max_i, max_j, matrix

# Optimise smith-waterman function
def optimise_sw(seq_1, seq_2, matrix, max_index, max_i, max_j):
    '''This function optimises the S-W algorithm by removing the gaps at the end of the alignment. 
    It is based on the fact that the S-W algorithm always starts from the cell with the highest score, 
    hence, the gaps at the end of the alignment are not relevant. 
    '''
    
    # Initialising the variables for tracing
    aligned_seq_1 = []
    aligned_seq_2 = []
    index_seq_1 = []
    index_seq_2 = []
    (max_i, max_j) = max_index
    
    # Tracing and computing the pathway with the local alignment
    while matrix[max_i, max_j] != 0:
        if matrix[max_i, max_j] == matrix[max_i - 1, max_j - 1] + Score.MATCH:
            aligned_seq_1.insert(0, seq_1.iloc[max_i - 1, 2])
            aligned_seq_2.insert(0, seq_2.iloc[max_j - 1, 2])

            index_seq_1.insert(0, seq_1.iloc[max_i - 1, 0])
            index_seq_2.insert(0, seq_2.iloc[max_j - 1, 0])

            max_i = max_i - 1
            max_j = max_j - 1

        elif matrix[max_i, max_j] == matrix[max_i - 1, max_j] + Score.FIRST_GAP:
            aligned_seq_1.insert(0, seq_1.iloc[max_i - 1, 2])
            aligned_seq_2.insert(0, '-')

            index_seq_1.insert(0, seq_1.iloc[max_i - 1, 0])
            index_seq_2.insert(0, '-')

            max_i = max_i - 1    

        elif matrix[max_i, max_j] == matrix[max_i, max_j - 1] + Score.FIRST_GAP:
            aligned_seq_1.insert(0, '-')
            aligned_seq_2.insert(0, seq_2.iloc[max_j - 1, 2])

            index_seq_1.insert(0, '-')
            index_seq_2.insert(0, seq_2.iloc[max_j - 1, 0])
            
            max_j = max_j - 1

    return aligned_seq_1, aligned_seq_2, index_seq_1, index_seq_2

# Detect if a gene has been duplicated
def duplicate(all):
    '''Evolutioary wise, a gene might be duplicated from one species to another.
    This function detects if a gene has been duplicated in the alignment done 
    by the S-W algorithm. It is based on the fact that if a gene is duplicated,
    it shows when there is a match and then a gap. This not happens in inverse
    order due to the scoring system, a match will be assigned before a gap,
    so only Match -> Gap with the same Gene as before is considered.
    
    Moreover, a gene might be duplicated not just once, but several times. 
    Hence, consecutive duplicates are checked when a duplicate is found.
    '''

    for i in range(1, len(all)):
        # Search for gaps and look before or after
        if all.loc[i, 'Result'] == 'Gap':
            if (all.loc[i-1, 'Result'] == 'Match') | (all.loc[i-1, 'Result'] == 'Inversion'):
                if (all.loc[i-1, 'Query Sequence Gene'] == all.loc[i, 'Query Sequence Gene']) | (all.loc[i-1,'Target Sequence Gene'] == all.loc[i,'Target Sequence Gene']):
                    all.loc[i, 'Result'] = 'Duplicate'

                    # Check for consecutives duplicates
                    j = i+1
                    while(
                        (j < len(all)) &
                        (all.loc[j, 'Result'] == 'Gap') & 
                        (all.loc[i,'Query Sequence Gene'] == all.loc[j,'Query Sequence Gene']) & 
                        (all.loc[i,'Target Sequence Gene'] == all.loc[j ,'Target Sequence Gene'])):
                            all.loc[j, 'Result'] = 'Duplicate'
                            j += 1
                    i = j-1
    return all


def merge(seq_1, seq_2, aligned_seq_1, index_seq_1, aligned_seq_2, index_seq_2) -> pd.DataFrame:
    '''S-W returns the aligned sequences and the indexes of the genes in the
    original sequences. However, all its information is lost. This function
    merges the aligned sequences with the original information of the genes,
    such as the replicon name, the replicon accession, the strand. Those
    are merged by the indexes of the genes in the original sequences, which
    are kept in the S-W function.
    '''

    # Merge the aligned sequences result from S-W with the genes information
    align_seq_1 = (pd.DataFrame({
                    'Original Position': index_seq_1,
                    'Gene': aligned_seq_1})
                    .merge(
                        seq_1[['index', '#Replicon Name', 'Replicon Accession', 'Strand', 'Start', 'Stop']],
                        left_on='Original Position',
                        right_on='index',
                        how='left'))
   
    align_seq_2 = (pd.DataFrame({
                    'Original Position': index_seq_2,
                    'Gene': aligned_seq_2})
                    .merge(
                        seq_2[['index', '#Replicon Name', 'Replicon Accession', 'Strand', 'Start', 'Stop']], 
                        left_on='Original Position',
                        right_on='index',
                        how='left'))

    align_seq_1.rename(columns={
                    'Original Position':'Target Sequence Original Position',
                    'Gene':'Target Sequence Gene',
                    '#Replicon Name':'Target Sequence Replicon Name',
                    'Replicon Accession':'Target Sequence Replicon Accession',
                    'Start':'Target Start', 
                    'Stop':'Target Stop',
                    'Strand': 'Target Sequence Orientation'
                    }, inplace=True)

    align_seq_2.rename(columns={
                    'Original Position':'Query Sequence Original Position',
                    'Gene':'Query Sequence Gene',
                    '#Replicon Name':'Query Sequence Replicon Name',
                    'Replicon Accession':'Query Sequence Replicon Accession',
                    'Start':'Query Start',
                    'Stop':'Query Stop',
                    'Strand': 'Query Sequence Orientation'
                    }, inplace=True)

    all = (align_seq_1.loc[:, align_seq_1.columns != 'index']
            .merge(
                align_seq_2.loc[:, align_seq_2.columns != 'index'],
                right_index=True, 
                left_index=True))           

    # Classify the alignment into Match, Mismatch, Inversion, Gap and Duplicate
    all['Result'] = np.where(all['Query Sequence Gene'] == all['Target Sequence Gene'],  
                        np.where(all['Query Sequence Orientation'] == all['Target Sequence Orientation'], 
                            'Match', 
                            'Inversion'),
                        np.where(((all['Query Sequence Gene'] == '-') | (all['Target Sequence Gene'] == '-')),
                            'Gap', 
                            'Mismatch'))

    # Reorder the columns
    all = all[[
            'Target Sequence Original Position',
            'Target Sequence Replicon Name',
            'Target Sequence Replicon Accession',
            'Target Start',
            'Target Stop',
            'Target Sequence Orientation',
            'Target Sequence Gene',
            'Result',
            'Query Sequence Gene',
            'Query Sequence Orientation',
            'Query Start',
            'Query Stop',
            'Query Sequence Replicon Accession', 
            'Query Sequence Replicon Name',
            'Query Sequence Original Position']]

    return all


# Check if a folder exists, if not, create it
def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Folder created at {path}")
    else:
        print(f"Folder already exists at {path}")

# Write the alignment results of scaffold alignment
def write_alignment(align_species, query_sp, query_sca, writer):

    # Excel sheet name limit is 31 characters
    # Check if we are over the limit and shorten the name if needed
    # But 3 characters are already used for the prefix
    query_name = query_sp + '_' + query_sca    
    if len(query_name) > 28:
        query_name = query_sp[0:3] + '_' + query_sca
    
    # Write the alignment results
    align_species.to_excel(writer, sheet_name = 'Al_' + query_name, index = False)

    # Restart the variable
    align_species = pd.DataFrame()