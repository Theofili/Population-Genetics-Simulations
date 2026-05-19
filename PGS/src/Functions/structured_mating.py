

import pandas
import numpy
import random


# Step 2b: create offsprings - structure

def mating_matrix_structured(df, mating_type, num_offsprings, gfr=None, output_path=None):

    """ a. Calculate mating chances, between random individuals based on mating type
        b. Create offsprings from the combination of two parent genotypes

        Args:

            df                (dataframe)     : Dataframe with the structured of the output of create_parents function

            mating_type       (str)           : Two options --> ('no gene flow', 'partial gene flow') 
                                                * partial gene flow requires gfr

            num_offsprings    (int)           : Number of offsprings ('Offspring_1', 'Offspring_2', etc)

            gfr               (float)         : if mating_type == 'partial gene flow', determine the frequency in which gene flow occurs

            output_path       (string)        : Folder to save csv file created
            OPTIONAL


        Returns:         

        offspring_df                          : Dataframe for offsprings and matched genotypes (structure: 'Male Parent', 'Female_Parent', 'Offspring_ID', 'Combined Marker 1', etc)

        gfr                                   : if mating_type == 'partial gene flow', used in other functions


"""

    # Crash test
    if num_offsprings <= 0:
        raise ValueError('num_offsprings must be > 0')
    if gfr is not None:
        if (0 > gfr > 1):
            raise ValueError('gene flow rate must be between 0 and 1')
    if df is None:
        raise ValueError('Parents dataset is required')

    individuals = len(df)
    loci = len(df.columns) - 1

    # Make arrays with parents - Hermaphrodite - so every individual is a male and female parent
    f = df['Individual_id'].values
    m = df['Individual_id'].values

    pop = numpy.full(individuals, None)

    # Divide individuals into populations

    mock_values = ["mock1", "mock2", "mock3"] # Maybe this should be user input based on number of individuals
    ind_count = individuals // len(mock_values)

    for i in range(len(mock_values)):
        start = i * ind_count
        end = start  + ind_count
        pop[start:end] = [mock_values[i]] * ind_count

    popf = pop.copy()
    popm = pop.copy()

    mating_all = numpy.empty((individuals + 2, individuals + 2), dtype=object)

    # Fill in the first row with female trees and the first column with male trees

    mating_all[0,2:] = popf # Mock numbers are added
    mating_all[1,2:] = f
    mating_all[2:,0] = popm # Mock numbers are added
    mating_all[2:,1] = m


    # Two user_inputs:  'no gene flow', 'partial gene flow'
    ## Fill the mating_all matrix with the probabilities of mating based on input

    if mating_type == 'no gene flow':

        for i in range(2, individuals + 2):
            for j in range(2, individuals + 2):
                if popf[i - 2] == popm[j - 2]:  # Check if both belong to the same population
                    mating_all[i, j] = 1
                else:
                    mating_all[i, j] = 0

    elif mating_type == 'partial gene flow':

        if gfr is None:
            raise ValueError('Provide gene flow rate')

        else:
            for i in range(individuals):
                for j in range(individuals):
                    if popf[i] == popm[j]:
                        mating_all[i+2,j+2] = 1 - gfr
                    else:
                        mating_all[i+2,j+2] = gfr

    else:

        raise ValueError("Invalid input. Please enter 'no gene flow' or 'partial gene flow' ")


    # Turn the table into a DataFrame and export csv

    mating_all_df = pandas.DataFrame(mating_all)
    if output_path is not None:
        mating_all_df.to_csv(f'{output_path}/mating_all_df_{mating_type}.csv', index=False)


    parents_male = []
    parents_female = []
    offspring_data = []


    # Select mating pairs based on the mating_all values
    x = 0

    if mating_type == 'no gene flow':

        while x < num_offsprings:
            i = random.randint(0, individuals - 1)
            j = random.randint(0, individuals - 1)

            # Skip if parents are from different populations
            if pop[i] != pop[j]:
                continue

            probability = float(mating_all_df.iat[i + 2, j + 2])

            # Mating allowed if probability == 1
            if probability == 1:
                parents_male.append(mating_all_df.iat[1, j + 2])
                parents_female.append(mating_all_df.iat[i + 2, 1])
                x += 1


    elif mating_type == 'partial gene flow':


        while x < num_offsprings:
            i = random.randint(0, individuals - 1)
            j = random.randint(0, individuals - 1)

            probability = float(mating_all_df.iat[i + 2, j + 2])

            # Mating allowed based on probability (includes flow possibility)
            if probability == 1 or random.random() <= probability:
                parents_male.append(mating_all_df.iat[1, j + 2])
                parents_female.append(mating_all_df.iat[i + 2, 1])
                x += 1

    else:
        print("Invalid mating_type. Use 'no gene flow' or 'partial gene flow'.")


    # Generate offspring for each selected pair

    for i in range(len(parents_male)):
        pm_id = parents_male[i]
        pf_id = parents_female[i]
        offspring_id = f'Offspring_{i + 1}'

        # Retrieve genetic markers for each parent

        male_markers = df[df['Individual_id'] == pm_id]
        female_markers = df[df['Individual_id'] == pf_id]

        # Check if both parents have valid markers

        if not male_markers.empty and not female_markers.empty:
            male_markers = male_markers.iloc[0, 1:].values
            female_markers = female_markers.iloc[0, 1:].values
            num_markers = len(male_markers)

            combined_markers = []

            # Combine alleles from both parents

            for j in range(num_markers):
                male_alleles = male_markers[j]
                female_alleles = female_markers[j]
                mr_allele = random.choice(male_alleles.split('/')) # mr_allele stands for male random allele
                fr_allele = random.choice(female_alleles.split('/')) # fr_allele stands for female random allele
                combined_markers.append(f'{mr_allele}/{fr_allele}')

            # Determine the "mock" group of the offspring based on the parents
            # In the case of partial gene flow there are flow products

            parent_index_m = numpy.where(f == pf_id)[0][0]
            parent_index_f = numpy.where(m == pm_id)[0][0]

            # Skip mismatched populations if no gene flow

            if mating_type == 'no gene flow' and pop[parent_index_m] != pop[parent_index_f]:
                continue


            if pop[parent_index_m] == pop[parent_index_f]:
                offspring_mock = pop[parent_index_m]
            else:
                offspring_mock = f'Flow Product {pop[parent_index_m]}'

            # Append the offspring data, including "mock" information

            offspring_data.append({
                'Male Parent': pm_id,
                'Female Parent': pf_id, 
                'Offspring ID': offspring_id,
                'Mock Group': offspring_mock,
                **{f'Combined Marker {j+1}': combined_markers[j] for j in range(len(combined_markers))}
            })

        else:
            print(f'Skipping mating pair {pm_id} and {pf_id} due to missing data')

    # Save the offspring data to a CSV file

    offspring_df = pandas.DataFrame(offspring_data)

    mark_cols = offspring_df.columns[4:]
    for col in mark_cols:
        a = offspring_df[col].str.split("/", expand=True)
        offspring_df[col] = a.min(axis=1) + "/" + a.max(axis=1)

    if output_path is not None:

        offspring_df.to_csv(f'{output_path}/offsprings_{mating_type}.csv', index=False)

    if mating_type == 'partial gene flow':
        return offspring_df, gfr
    else:
        return offspring_df