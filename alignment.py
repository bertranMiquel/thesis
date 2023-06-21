# This python script will use some functions to perform the alignment of the annotations.

# Import the libraries
import pandas as pd
import time
import numpy as np
# import xlsxwriter
from matplotlib import pyplot as plt
import logging
import argparse

# Import the functions from the other scripts
import k_mers.k_mers_metric as kmm
import align.swco as swco
import annotations.annotation_length as al

# This function will perform the alignment of the annotations
# We will align the target species or scaffold vs. all the query species or scaffold
# Always in this direction: target vs. query
# We will use the Smith-Waterman algorithm
# If those variables are empty, we will align all the species/scaffolds vs. all the species/scaffolds given in the dataframe
def alignment(mode, df_k_mers: pd.DataFrame, k_mers_counter: pd.DataFrame, 
                df_replicon_accession_length: pd.DataFrame, df_annotations: pd.DataFrame,
                target_species=[], target_scaffold=[], query_species=[], query_scaffold=[], # relevant parameters
                output_path='../Data/Intermediate/alignment/'):

    # Create the output_folder where the alignment results and plots will be stored
    # Add the timestamp to the output folder path
    output_path = swco.create_output_folder(output_path)

    # If the target_species is empty, we will align all the species vs. all the species
    if not target_species:
        target_species = df_annotations['species'].unique()
    else:
        print('Target species: ', target_species)
        # Check if the target species are in the dataframe
        if not target_species in df_annotations['species'].unique():
            raise Exception('The target species introduced do not belong to the dataframe')

    # If the target_scaffold is empty, we will align all the scaffolds of the target_species vs. all the scaffolds of the query_species
    if not target_scaffold:
        target_scaffold = df_annotations.loc[df_annotations['species'].isin(target_species), 'replicon_accession'].unique()
    else:
        # Check if the target scaffold belongs to the target species
        all_scaffold = df_annotations.loc[df_annotations['species'].isin(target_species), 'replicon_accession'].unique().tolist()
        all_present = all(item in all_scaffold for item in target_scaffold)
        if not all_present:
            raise Exception('The target scaffolds introduced do not belong to one of the target species')

    # For visualization purposes, initialise some variables
    blocks = []
    align_counter = 0

    # Smith-Waterman algorithm has some parameters
    # Let's store them to write them in the output file
    # The parameters can be modified in the swco.py script
    params = pd.DataFrame({
        'Description': ['', 'Parameters', 'Match', 'Mismatch', 'First gap', 'Consecutive gap'], 
        'Value': ['', '', swco.Score.MATCH, swco.Score.MISMATCH, swco.Score.FIRST_GAP, swco.Score.CONS_GAP]})

    # When we do the alignment, we would like to keep the position of the genes in the table
    # So, we will create a new column with the index of the row
    df_annotations['index'] = df_annotations.index

    for target_sp in target_species:

        # Define the query_species
        if not query_species:
            query_species = df_annotations.loc[df_annotations['species'] != target_sp, 'species'].unique()
        # Check that target_sp is not in query_species
        elif target_sp in query_species:
            query_species.remove(target_sp)
        
        # Define the query_scaffold
        if not query_scaffold:
            query_scaffold = df_annotations.loc[df_annotations['species'].isin(query_species), 'replicon_accession'].unique()
        else:
            # Sanity check: verify that the query_scaffold belong to the query_species
            query_scaffold = df_annotations.loc[df_annotations['species'].isin(query_species) & df_annotations['replicon_accession'].isin(query_scaffold), 'replicon_accession'].unique()

        # For each target_scaffold
        for target_sc in target_scaffold:
            
            logging.info(f'Aligning {target_sp} {target_sc}')

            # Define target as the annotations of this target scaffold
            df_target = df_annotations.loc[(df_annotations['species'] == target_sp) & (df_annotations['replicon_accession'] == target_sc)]

            # Now, we want to order the query scaffolds per species by similarity to the target scaffold
            # So, we will first compare the target scaffold with the most similar query scaffolds
            # Then, we will compare the target scaffold non-matching part with the second most similar query scaffolds and so on
            # We have implemented a metric in order to compare them based in the k-mers
            # This is implemented in the k_mers_metric.py script
            
            # Filter to keep only the relevant query scaffolds for this alignment
            k_mers_counter_query = k_mers_counter.loc[(k_mers_counter['replicon_accession'].isin(query_scaffold)) | (k_mers_counter['replicon_accession'] == target_sc), ['species', 'replicon_accession', 'k', 'total_k_mers']]


            # Return the similar scaffolds ordered by species and similarity descending
            df_scaffold_cand = pd.DataFrame()
            logging.info(f'Calculating similar scaffolds for {target_sp} {target_sc}')
            df_scaffold_cand = kmm.k_mers_to_metric(target_sc, level='replicon_accession', k_mers_counter=k_mers_counter_query, df_length=df_replicon_accession_length, df_k_mers=df_k_mers)
            if df_scaffold_cand.empty:
                logging.warning(f'No similar scaffolds for {target_sp} {target_sc}')
                continue

            # Target scaffold is compared versus the similarest scaffold of one specie
            # Then, the non-matching part of the target scaffold is compared versus the second similarest scaffold of the same specie
            # So, our df_target is going to be updated after each alignment, if we compare to the same species
            # Then, we need a variable to keep the original df_target and to check if we change of query species
            df_target_original = df_target.copy()
            query_sp_old = ''

            # Create an excel file to store the alignment results per each target scaffold
            writer = pd.ExcelWriter(f'{output_path}/{target_sp}_{target_sc}.xlsx', engine='xlsxwriter')

            # Write the results and the align species variables are needed
            headers = ['Target Species', 'Target Scaffold', 'Query Species',
             'Query Scaffold', '% Covered Target sequence', '% Covered Query sequence', 'Number of aligned genes']
            
            align_species = pd.DataFrame()

            results = []

            # Go through the list of similar scaffolds to apply Smith-Waterman algorithm
            for query_sp, query_sc in zip(df_scaffold_cand['species'], df_scaffold_cand['replicon_accession']):
                
                # Keep the similarity ratio with us
                sim_ratio = df_scaffold_cand.loc[(df_scaffold_cand['species'] == query_sp) & (df_scaffold_cand['replicon_accession'] == query_sc), 'metric'].values[0]
                logging.info(f'Aligning {target_sp} {target_sc} vs. {query_sp} {query_sc} with similarity ratio: {sim_ratio:0.4f}')
                # I do not get why this is not printing correctly sim_ratio. If I tried appart, it works.
                
                # If we change of query species, we need to update the df_target
                if query_sp != query_sp_old:
                    df_target = df_target_original.copy()
                    query_sp_old = query_sp

                    # Write the alignment results of the previous species
                    if not align_species.empty:
                        swco.write_alignment(align_species, query_sp, query_sc, writer)

                # Define query as the annotations of this query scaffold
                df_query = df_annotations.loc[(df_annotations['species'] == query_sp) & (df_annotations['replicon_accession'] == query_sc)]
                
                # Increase the align counter
                align_counter += 1

                # Initialise variables for the alignment
                # Reset them, so no values from the past alignment are kept
                aligned_query = []
                aligned_target = []
                index_query = []
                index_target = []
                matrix = np.zeros((len(df_target), len(df_query)), dtype=int)

                ## Apply Smith-Waterman algorithm
                # The results are stored in a dataframe
                tic = time.perf_counter()
                aligned_target, aligned_query, index_target, index_query, max_score, max_index, max_i, max_j, matrix = swco.smith_waterman(df_target, df_query, mode=mode)
                toc = time.perf_counter()

                # Time
                logging.info(f'Alignment done in {toc - tic:0.4f} seconds')

                tic = time.perf_counter()

                ## Post processing
                # From here on, we will process the results of the alignment, transforming them into some format to report the main information.
                # That said, everything can be changed, except the blocks csv file which is used afterwards for the alignment visualization.

                # If the alignment is empty, we will not store the results
                if not aligned_target:
                    logging.info(f'No alignment found for {target_sp} {target_sc} vs. {query_sp} {query_sc}')
                    continue
                
                # Smith Waterman algorithm returns a target and a query aligned sequences, getting back the index and the locus. 
                # Merge this sequences to be able to represent them in a table
                # Moreover, merge the aligned sequences with the master information from the original sequences, such as, replicon_accession, start, end, strand, etc.
                align = pd.DataFrame()
                align = swco.merge(df_target, df_query, aligned_target, index_target, aligned_query, index_query)

                # | target locus | Result | query locus |
                # |--------------|--------|-------------|
                # | ABBB         |  Match | ABBB        |
                # | ABBB         |  Gap   | -           |
                # 
                # Scenarios like this should be tagged as a Duplicate genes, instead of a Gap
                # There is a function to do this
                # Look for the genes which have been duplicated and tag them as so, instead of tagging it as a Gap
                align = swco.duplicate(align)

                # Keep the results of the alignement per species
                align_species = pd.concat([align_species, align])

                toc = time.perf_counter()

                d = toc-tic
                logging.info(f'Postprocessed in {d:0.4f} seconds')

                ## Write the alignment results
                tic = time.perf_counter()

                # Create an aggregation of the results information and write it in the excel workbook
                # For this, we need all the information regarding the alignment and the two aligned sequences
                # Gives us back the index information for later usage
                # Also, another dataframe called blocks for doing the alignment visualization
                len_target = len(df_target)
                len_query = len(df_query)
                
                actual_results, blocks, index_target_no_gaps = swco.summary_results(writer, aligned_query, aligned_target, index_query, index_target, align, target_sp, target_sc, query_sp, query_sc, len_target, len_query, sim_ratio, params, blocks, align_counter)
                
                results.append(actual_results)

                # Create a heatmap with the matrix information to represent which part of the sequences are more similar
                swco.heatmap(writer, matrix, target_sp, target_sc, query_sp, query_sc, output_path)
                
                toc = time.perf_counter()

                d = toc-tic
                logging.info(f'Written in {d:0.4f} seconds')

                # Remove the already aligned parts for the target scaffold
                # For doing so, we need the index of the target scaffold already mapped
                df_target = swco.remove_match_parts(df_target, index_target_no_gaps)

                # Clean variables
                del aligned_query, aligned_target, index_query, index_target, matrix, align

            # Write the last alignment results
            swco.write_alignment(align_species, query_sp, query_sc, writer)

            # Write the results dataframe
            df_results = pd.DataFrame(results, columns=headers)
            df_results.to_excel(writer, sheet_name = 'Results', index = False)
            writer.sheets['Results'].activate()
            writer.close()   

    # Write the blocks dataframe
    # If there are no blocks, we will not write the file
    if not blocks:
        logging.info('No blocks aligned')
        
    else:
        logging.info('Writing blocks.csv')
        # Before saving blocks, let's add column names to the list
        blocks = pd.DataFrame(blocks)
        # Make them the column names
        blocks.columns = ['species', 'replicon_accession', 'target_replicon_accession', 'comparing_species', 'alignment_id', 'start', 'stop', 'match_perc']

        blocks.to_csv(output_path + '/blocks.csv', index = False)


if __name__ == '__main__':
    # Set the arguments
    parser = argparse.ArgumentParser(description='Smith-Waterman algorithm for aligning sequences')
    parser.add_argument('-m', '--mode', type=str, help='Mode of the algorithm', default='strict')
    parser.add_argument('-tsp', '--target_species', nargs='+', help='Target species', default=[])
    parser.add_argument('-tsc', '--target_scaffold', nargs='+', help='Target scaffold', default=[])
    parser.add_argument('-qsp', '--query_species', nargs='+', help='Query species', default=[])
    parser.add_argument('-qsc', '--query_scaffold', nargs='+', help='Query scaffold', default=[])
    parser.add_argument('-o', '--output_path', type=str, help='Output path', default='../Data/Intermediate/alignment/')
    args = parser.parse_args()

    # Set the logging configuration
    logging.basicConfig(
        format='%(asctime)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S',
        level=logging.INFO,
    )

    # Read the annotations of all the species
    df_annotations = pd.read_feather('../Data/Intermediate/accumulated/df_annotations.feather')

    # Read the length of each scaffold
    df_replicon_accession_length = al.read_length('replicon_accession')
    
    # Remove small scaffolds
    # Ivan suggest to make the threshold equal to 10
    df_annotations = al.remove_short_replicons(df_annotations, df_replicon_accession_length, 10)

    # Read the k-mers for all the species
    df_k_mers = pd.read_feather('../Data/Intermediate/accumulated/df_k_mers_relaxed.feather')

    # Read the number of k-mers per scaffold
    k_mers_counter = pd.read_feather('../Data/Intermediate/accumulated/relaxed_count_k_mers_replicon_accession.feather')

    # Perform the alignment
    # There should be another way to modify the important input parameters, 
    # such as target_species, query_species, target_scaffold, query_scaffold & mode in a more user friendly way,
    # instead of going down here and execute that
    # I do not know how to do it yet
    alignment(target_species=args.target_species, target_scaffold=args.target_scaffold, 
            query_species=args.query_species, query_scaffold=args.query_scaffold, mode=args.mode,
            df_k_mers=df_k_mers, k_mers_counter=k_mers_counter, output_path=args.output_path,
            df_replicon_accession_length=df_replicon_accession_length, df_annotations=df_annotations)

