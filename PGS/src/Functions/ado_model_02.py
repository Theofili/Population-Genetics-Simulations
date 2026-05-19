

import pandas
import numpy
import random

# Step 4b: Simulate allelic dropout


def ado_model_multiple(input_path, mating_type, output_paht=None):

    """Calculate how the results would look if allelic dropout occured in multiple samples

    Args:

        input_path        (str)     : File path to 'multiple_runs_save.csv', created from previous function on the pipeline

        mating_type       (str)     : Five options --> ('no selfing', 'selfing only', 'mixed', 'no gene flow', 'partial gene flow') 

        output_path       (str)     : Folder to save csv file created
        OPTIONAL   
    
    Returns:

        df             (dataframe)  : Dataframe containing genotype frequencies from populations with the null alleles     
"""

    if input_path is not None:
        path = (f'{input_path}/multiple_runs_save.csv')
    else:
        path = ('multiple_runs_save.csv')

    # FIND AND REPLACE FUNCTION

    import re

    # open your csv and read as a text string
    with open(f'{path}', 'r') as f:
        my_csv_text = f.read()


    find_str = 'A0'
    replace_str = 'NULL'


    # substitute
    new_csv_str = re.sub(find_str, replace_str, my_csv_text)

    if input_path is not None:
        path_01 = (f'{input_path}/NULL_multiple_runs_save.csv')
    else:
        path_01 = ('NULL_multiple_runs_save.csv')
    # open new file and save
    new_csv_path = f'{path_01}'
    with open(new_csv_path, 'w') as f:
        f.write(new_csv_str)

    df_expected = pandas.read_csv(f'{input_path}/NULL_multiple_runs_save.csv')
    df = pandas.read_csv(f'{input_path}/NULL_multiple_runs_save.csv')
    
    
 
    # Go through the dataframe and find the frequency of genotypes that contain NULL
    rep = len(df.columns) #

    for row in range(0, len(df), 2):

        #print(f'Prossecing row {row}....')

        # Get genotype and frequency
        if mating_type == 'mixed' or  mating_type == 'partial gene flow':
            start_col = 6
        else:
            start_col = 5

        for col in df.columns[start_col:rep]:
            genotype = df.at[row, col]
            #print(genotype)

            freq = df.at[row+1,col]

            # If we are at number skip

            if isinstance(df.at[row, col], (float, numpy.floating)):
                continue

            # Homozygous NULL genotype needs to be recorded to calculate the missing offsprings number

            if genotype == 'NULL/NULL':
                #print(f'Found homozygous Null genotype with frequency:{freq}')
                num_off = df.at[row, 'Num_Offsprings']
                mul_num = float(num_off)
                mon = float(freq)*mul_num
                #print(f'Missing offsping number is: {int(mon)}')
                df.at[row+1,col] = int(mon)
                #print('Creating new column')
                df.at[row, rep + 1] = mul_num - int(mon)
                continue

            # Heterozygous NULL genotypes are added to the other allele homozygous genotype

            if genotype.startswith("NULL/"):
                #print(f'Working on Null genotype {genotype}, frequency:({freq})')
                phen_gen = genotype.split('/')[1]
                homozygous = f'{phen_gen}/{phen_gen}'
                #print(f'Heterogynous Null frequency will be added to: {homozygous}')

                if homozygous in df.iloc[row, 5:rep].values:
                    locate_col = df.columns[df.iloc[row] == homozygous][0]
                    col_idx = df.columns.get_loc(locate_col)
                    old_freq = float(df.iat[row+1, col_idx])
                    new_freq = (float(freq) + old_freq)
                    df.iloc[row+1,col_idx] = new_freq
                    df.at[row, col] = 'NA'
                    df.at[row+1, col] = 'NA'
                    #print(f'Homozygous genotype new frequency {homozygous}: {old_freq} + {freq} = {new_freq}')
                    continue
                else:
                    #print('null_created!')
                    locate = df.columns[df.iloc[row] == 'NA'][0]
                    df.at[row, locate] = f'{homozygous}(null_created)'
                    df.at[row+1, locate] = freq
                    df.at[row, col] = 'NA'
                    df.at[row+1, col] = 'NA'
                    #print(f'Added {homozygous} to row:{row}')

    df.replace('NA', numpy.nan)
    df.fillna('')


    if output_paht is not None:
        df.to_csv(f'{output_paht}/Misleading_results1.csv')
    return df_expected, df