# Functions for comparison metric computation

# Summary:
# - read_k_mers: read all the k_mers files for an specific path
# - read_k_mers_prime: read only the k_mers files for k prime for an specific path
# - ratio: compute the ratio for given formula
# - metric_comp: given two species, computes its similarity value, based on k_mers files
# - measure_to_dist: takes the measure matrix and transform it to a distance by applying exponential formulas
# - metric_to_matrix: convert the metric as a dataframe to a simetric matrix

import pandas as pd
import numpy as np
import os
from sympy import isprime

def read_k_mers(path):
    k_mers = pd.DataFrame()
    
    for file_name in os.listdir(path):
        
        full_path = os.path.join(path, file_name)
        print(k, ': ', file_name)
        k, _ = file_name.split('_', 1)
        k = int(k)

        temp_df = pd.DataFrame()
        temp_df = pd.read_csv(file_name).rename(columns={str(k) + '_mers':'k_mers'})
        temp_df['k'] = int(k)

        k_mers = pd.concat([k_mers, temp_df])

    return k_mers

def read_k_mers_prime(path):

    k_mers = pd.DataFrame()

    for file_name in os.listdir(path):
        
        full_path = os.path.join(path, file_name)
        k, _ = file_name.split('_', 1)
        k = int(k)

        if isprime(k):
            
            print(k, ': ', file_name)

            temp_df = pd.DataFrame()
            temp_df = pd.read_csv(full_path).rename(columns={str(k) + '_mers':'k_mers'})
            temp_df['k'] = int(k)

            k_mers = pd.concat([k_mers, temp_df])

    return k_mers

def ratio(data, m, n):
    return 1 - np.nanprod(np.array(1 - (data['Both'] / (data['Both'] + data['One']) ) ) * (1 - (m/(n*data['k']**2)))) 

def metric_comp(specie1, specie2, data, m, n):

    # Just consider the k_mers contained in the two species of interest
    # Is faster like this
    data_pos = data[data['Specie_Scaffold'].str.contains(specie1, regex=True) | data['Specie_Scaffold'].str.contains(specie2, regex=True)].copy()
    del data
    
    # Check which k_mers are in the both species
    both = (data_pos['Specie_Scaffold'].str.contains(specie1, case=False, regex=True)) & (data_pos['Specie_Scaffold'].str.contains(specie2, case=False, regex=True))
    
    # In case common k_mers are found, compute the ratio
    if not data_pos[both].empty == True:

        data_pos.loc[both == False, 'Common'] = 'One'
        data_pos.loc[both == True, 'Common'] = 'Both'
        
        # Count how many k_mers we have for these species specifying which appear in Both or just in One specie
        data_group = (data_pos.groupby(['k', 'Common'], as_index=False)['Scaffolds']
                    .count()
                    .pivot(index='k', columns='Common', values='Scaffolds')
                    .reset_index())  

        # In case all k_mers are common
        if (both == True).all() == True:
            data_group['One'] = 0 
        
        # m smaller, n larger
        if n < m: 
            aux = m
            m = n
            n = aux

        return ratio(data_group, m, n)
    
    # Not common k_mers found
    else: return 0

def measure_to_dist(z, a):
    return np.exp(1/z**a) - np.exp(1)

def metric_to_matrix(metric):
    metric2 = pd.DataFrame()

    metric2['Specie1'] = metric['Specie2']
    metric2['Specie2'] = metric['Specie1']

    mask = (metric2['Specie1'] == metric['Specie2']) & (metric2['Specie2'] == metric['Specie1'])

    metric2 = metric2.merge(metric, left_on=['Specie1', 'Specie2'], right_on=['Specie2', 'Specie1'], how='left')

    metric2.rename(columns={'Specie1_x':'Specie1', 'Specie2_x':'Specie2'}, inplace=True)
    del metric2['Specie1_y'], metric2['Specie2_y']

    metric_all = pd.concat([metric, metric2])
    
    return metric_all.pivot(index='Specie1', columns='Specie2')