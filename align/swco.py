'''
All functions defined for applying the Smith Watermann algorithm & writting its results in the created output folder
'''
import pandas as pd
import numpy as np
from enum import IntEnum
import os
import matplotlib.pyplot as plt
import seaborn as sns
import time
from datetime import datetime as dt
import logging

logging.basicConfig(
    format='%(asctime)s - %(message)s',
    datefmt='%Y-%m-%d %H:%M:%S',
    level=logging.INFO,
)


# Create the ouput folder if it does not exist
def create_output_folder(output_path = '../Data/Intermediate/alignment/') -> str: # Actually is a path
    # Create the output folder
    timestamp = dt.now().strftime('%Y%m%d_%H%M%S')
    output_path = output_path + timestamp + '/'

    # Create the output folder
    os.mkdir(output_path)
    logging.info(f'Output folder created: {output_path}')

    # Create the output folder for the plots
    os.mkdir(output_path + '/plots/')
    logging.info(f'Output folder for the plots created: {output_path + "/plots/"}')

    return output_path

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
      

def merge(seq_1, seq_2, aligned_seq_1, index_seq_1, aligned_seq_2, index_seq_2) -> pd.DataFrame:
    '''S-W returns the aligned sequences and the indexes of the locuss in the
    original sequences. However, all its information is lost. This function
    merges the aligned sequences with the original information of the locuss,
    such as the replicon_name, the replicon_accession, the strand. Those
    are merged by the indexes of the locuss in the original sequences, which
    are kept in the S-W function.
    '''

    # Merge the aligned sequences result from S-W with the locuss information
    align_seq_1 = (pd.DataFrame({
                    'original_position': index_seq_1,
                    'locus': aligned_seq_1})
                    .merge(
                        seq_1[['index', 'replicon_name', 'replicon_accession', 'strand', 'start', 'stop']],
                        left_on='original_position',
                        right_on='index',
                        how='left'))
   
    align_seq_2 = (pd.DataFrame({
                    'original_position': index_seq_2,
                    'locus': aligned_seq_2})
                    .merge(
                        seq_2[['index', 'replicon_name', 'replicon_accession', 'strand', 'start', 'stop']], 
                        left_on='original_position',
                        right_on='index',
                        how='left'))

    align_seq_1.rename(columns={
                    'original_position':'Target Sequence original_position',
                    'locus':'Target Sequence locus',
                    'replicon_name':'Target Sequence replicon_name',
                    'replicon_accession':'Target Sequence replicon_accession',
                    'start':'Target start', 
                    'stop':'Target stop',
                    'strand': 'Target Sequence strand'
                    }, inplace=True)

    align_seq_2.rename(columns={
                    'original_position':'Query Sequence original_position',
                    'locus':'Query Sequence locus',
                    'replicon_name':'Query Sequence replicon_name',
                    'replicon_accession':'Query Sequence replicon_accession',
                    'start':'Query start',
                    'stop':'Query stop',
                    'strand': 'Query Sequence strand'
                    }, inplace=True)

    all = (align_seq_1.loc[:, align_seq_1.columns != 'index']
            .merge(
                align_seq_2.loc[:, align_seq_2.columns != 'index'],
                right_index=True, 
                left_index=True))           

    # Classify the alignment into Match, Mismatch, Gap and Duplicate
    all['Result'] = np.where(all['Query Sequence locus'] == all['Target Sequence locus'],
                        np.where(((all['Query Sequence locus'] == '-') | (all['Target Sequence locus'] == '-')),
                            'Gap', 'Match'),
                            'Mismatch')    #np.where(all['Query Sequence strand'] == all['Target Sequence strand'], 'Match','Inversion'),

    # Reorder the columns
    all = all[[
            'Target Sequence original_position',
            'Target Sequence replicon_name',
            'Target Sequence replicon_accession',
            'Target start',
            'Target stop',
            'Target Sequence strand',
            'Target Sequence locus',
            'Result',
            'Query Sequence locus',
            'Query Sequence strand',
            'Query start',
            'Query stop',
            'Query Sequence replicon_accession', 
            'Query Sequence replicon_name',
            'Query Sequence original_position']]

    return all

# Detect if a gene has been duplicated, in an scenario like this
# | target locus | Result | query locus |
# |--------------|--------|-------------|
# | ABBB         |  Match | ABBB        |
# | ABBB        |  Gap    | -           |
# | ABBB        |  Gap    | -           |

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
                if (all.loc[i-1, 'Query Sequence locus'] == all.loc[i, 'Query Sequence locus']) | (all.loc[i-1,'Target Sequence locus'] == all.loc[i,'Target Sequence locus']):
                    all.loc[i, 'Result'] = 'Duplicate'

                    # Check for consecutives duplicates
                    j = i+1
                    while(
                        (j < len(all)) &
                        (all.loc[j, 'Result'] == 'Gap') & 
                        (all.loc[i,'Query Sequence locus'] == all.loc[j,'Query Sequence locus']) & 
                        (all.loc[i,'Target Sequence locus'] == all.loc[j ,'Target Sequence locus'])):
                            all.loc[j, 'Result'] = 'Duplicate'
                            j += 1
                    i = j-1
    return all

# Check if a folder exists, if not, create it
def create_folder(path):
    if not os.path.exists(path):
        os.makedirs(path)
        print(f"Folder created at {path}")
    else:
        print(f"Folder already exists at {path}")

# Write the alignment results of scaffold alignment
def write_alignment(align_species, query_sp, query_sc, writer):
    # Excel sheet name limit is 31 characters
    # Check if we are over the limit and shorten the name if needed
    # But 3 characters are already used for the prefix
    query_name = query_sp + '_' + query_sc 
    if len(query_name) > 27:
        query_name = query_sp[0:3] + '_' + query_sc
    
    # Write the alignment results
    align_species.to_excel(writer, sheet_name = 'Al_' + query_name, index = False)

    # Restart the variable
    align_species = pd.DataFrame()


# Create a summary with some main KPIs of the alignment
# Store it into the excel workbook
def summary_results(writer, aligned_query, aligned_target, index_query, index_target, align, 
                    target_sp, target_sc, query_sp, query_sc, 
                    len_target, len_query, score, params, blocks, align_counter):

    # Keep the values of the aligned indeces
    index_query_no_gaps = []
    index_target_no_gaps = []

    # Remove gaps
    index_query_no_gaps = list(filter(lambda c: c!= '-', index_query))
    index_target_no_gaps = list(filter(lambda c: c!= '-', index_target))

    # Convert their elements to integers
    index_query_no_gaps = [int(i) for i in index_query_no_gaps]
    index_target_no_gaps = [int(i) for i in index_target_no_gaps]

    target_name = target_sp + '_' + target_sc
    query_name = query_sp + '_' + query_sc

    # Summary results of the alignement: include some information + parameters + KPIs
    summary = []
    summary = pd.DataFrame({
        'Description' : ['Target Sequence', 'Length Target Sequence', 'Query Sequence', 'Length Query Sequence', 'Similarity Ratio'],     #, '', 
                        #'Aligned Target Sequence length', 'Aligned Query Sequence length'],
        'Value': [target_name, len_target, query_name, len_query, score]        #, '',
                    #len(index_target_no_gaps), len(index_query_no_gaps)]
    })

    summary = pd.concat([summary, params])

    aligned_query_no_gaps = []
    aligned_target_no_gaps = []
    aligned_query_no_gaps = list(filter(lambda c: c!= '-', aligned_query))
    aligned_target_no_gaps = list(filter(lambda c: c!= '-', aligned_target)) 

    #KPIs
    perc_align_query = len(aligned_query_no_gaps) / len_query
    perc_align_target = len(aligned_target_no_gaps) / len_target

    kpi = []
    kpi = pd.DataFrame({
        'Description': ['','KPIs','Length Aligned Target Sequence', '% Covered Target sequence', 'Length Aligned Query Sequence', '% Covered Query sequence', 'Length Aligned Sequence', ''], #include where duplicates occur
        'Value': ['','', len(aligned_query_no_gaps), perc_align_target, len(aligned_target_no_gaps), perc_align_query, align.shape[0], '']
    })

    result_perc = align.groupby('Result').count() / align.shape[0]
    result_perc = result_perc.reset_index().iloc[:,[0,1]]
    result_perc.rename(columns = {
                        'Result':'Description', 
                        'Target Sequence original_position':'Value'},
                        inplace=True)

    summary = pd.concat([summary, kpi, result_perc])

    # Write
    # Excel sheet name limit is 31 characters
    # Check if we are over the limit and shorten the name if needed
    # But 3 characters are already used for the prefix
    query_name_excel = query_name
    if len(query_name_excel) > 27:
        query_name_excel = query_sp[0:3] + '_' + query_sc
    
    # Write the alignment results
    summary.to_excel(writer, sheet_name = 'Sum_' + query_name_excel, index = False)

    # Keep track of the results among all the scaffolds
    actual_results = pd.DataFrame(data = {target_sp, target_sc, query_sp, query_sc, perc_align_target, perc_align_query, align.shape[0]})
    
    # Keep the alignment results for alignment visualization
    # We want to keep the starting and stopping positions of the alignment in the target and query sequences
    # target_sca is introduced in both of them in order to aggrupate them in the plot
    # To report the match percentage, the value of the match is taken from the result_perc dataframe
    # Not always we have a match, so we need to check if the value is present and keep it in a variable
    # If not, we assign a 0 value
    if result_perc.loc[result_perc['Description'] == 'Match', 'Value'].empty:
        match_perc = 0
    else:
        match_perc = result_perc.loc[result_perc['Description'] == 'Match', 'Value'].values[0]

    # One with the target and one with the query ifnormation
    blocks.append([target_sp, target_sc, target_sc, [target_sp, query_sp], align_counter, min(align['Target start'].astype('float')), max(align['Target stop'].astype('float')), match_perc])
    blocks.append([query_sp, query_sc, target_sc, [target_sp, query_sp], align_counter, min(align['Query start'].astype('float')), max(align['Query stop'].astype('float')), match_perc])


    return actual_results, blocks, index_target_no_gaps


# Create a heatmap with the matrix information to represent which part of the aligned sequences are more similar
def heatmap(writer, matrix, target_sp, target_sc, query_sp, query_sc, output_path) -> None:
    # Heatmap of the alignment
    heatmap = sns.heatmap(matrix)
    plt.title('Heatmap of ' + target_sp + target_sc + ' vs. ' + query_sp + query_sc, fontsize = 12)
    plt.savefig(output_path + '/plots/' + query_sp + query_sc + '.png')
    plt.clf()

    # Define the query name for the sheet name, as it is definied above when saving summary information
    query_name = query_sp + '_' + query_sc
    if len(query_name) > 27:
        query_name = query_sp[0:3] + '_' + query_sc

    workbook  = writer.book
    worksheet = writer.sheets['Sum_' + query_name]
    worksheet.insert_image('H2', output_path + '/plots/' + query_sp + query_sc + '.png')

# When computing the alignment, we first compute the target scaffold with the most similar scaffold per species
# Then, if there are more than one similar scaffold per species, we will compare the not aligned part of the target scaffold
# versus the second most similar scaffold per species
# Hence, we have to remove the part of the target scaffold that has been already aligned
# We do it by using this function and the index of the target scaffold without gaps
def remove_match_parts(df_target, index_target_no_gaps):
    prev = []
    after = []
    prev = df_target.loc[df_target['index'].astype('int') < index_target_no_gaps[0] + 1]
    after = df_target.loc[index_target_no_gaps[-1] - 1 < df_target['index'].astype('int')]

    return pd.concat([prev, after])