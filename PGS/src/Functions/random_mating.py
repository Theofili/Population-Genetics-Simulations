import pandas
import numpy
import random 


# Step 2a: create offsprings - random 

def mating_matrix_random(df, mating_type, num_offsprings, inbreeding_index=None, output_path=None):

    """ a. Calculate mating chances, between random individuals based on mating type
        b. Create offsprings from the combination of two parent genotypes

        Args:

            df                (dataframe)     : Dataframe with the structured of the output of create_parents function

            mating_type       (str)           : Three options --> ('no selfing', 'selfing only', 'mixed') 
                                                * mixed requires inbreeding_index

            num_offsprings    (int)           : Number of offsprings ('Offspring_1', 'Offspring_2', etc)

            inbreeding_index  (float)         : if mating_type == 'mixed', determine the frequency in which inbreeding occurs

            output_path       (string)        : Folder to save csv file created
            OPTIONAL


        Returns:         

        offspring_df                          : Dataframe for offsprings and matched genotypes (structure: 'Male Parent', 'Female_Parent', 'Offspring_ID', 'Combined Marker 1', etc)

        inbreeding_index                      : if mating_type == 'mixed', used in other functions


    """

    # Crash test
    if num_offsprings <=0:
        raise ValueError('num_offsprings must be > 0')
    if inbreeding_index is not None:
        if (0 > inbreeding_index > 1):
            raise ValueError('inbreeding index must be between 0 and 1')
    if df is None:
        raise ValueError('Parents dataset is required')

    individuals = len(df)
    loci = len(df.columns) - 1


    # Make arrays with parents - Hermaphrodite - so every individual is a male and female parent

    f = df['Individual_id'].values
    m = df['Individual_id'].values

    mating_all = numpy.empty((individuals + 1, individuals + 1), dtype=object)

    # Fill in the first row with female trees and the first column with male trees

    mating_all[0, 1:] = f
    mating_all[1:, 0] = m
    mating_all[0, 0] = 'Random'

    # Three user_inputs:  'no selfing', 'selfing only', 'mixed'
    ## Fill the mating_all matrix with the probabilities of mating based on input

    if mating_type == 'no selfing':

        mating_all[1:, 1:] = 1

        for i in range(1, individuals + 1):
            if mating_all[i, 0] == mating_all[0, i]:
                mating_all[i, i] = 0

    elif mating_type == 'selfing only':

        mating_all[1:, 1:] = 0

        for i in range(1, individuals + 1):
            if mating_all[i, 0] == mating_all[0, i]:
                mating_all[i, i] = 1

    elif mating_type == 'mixed':

        ## calculate selfing_percentage##
      
        selfing_percentage = 2 * inbreeding_index
        
        if inbreeding_index is None:
            raise ValueError('Provide inbreeding index')

        else:
            mating_all[1:, 1:] = 1  # every mate has the probability of one???

    else:
        raise ValueError("Invalid input. Please enter 'selfing only', 'no selfing', or 'mixed'.")
        return

    # Turn the table into a DataFrame and export csv

    mating_all_df = pandas.DataFrame(mating_all)

    if output_path is not None:
        mating_all_df.to_csv(f'{output_path}/mating_all_df_{mating_type}.csv', index=False)


    parents_male = []
    parents_female = []
    offspring_data = []

        # Select mating pairs based on the mating_all values

    if mating_type == 'selfing only' or mating_type == 'no selfing':
        x=0

        while(x < num_offsprings): # create random mating pairs to comple offspring numebr

            i = random.randint(1, individuals)
            z = random.randint(1, individuals)

            probability = float(mating_all[i, z])

            # If the probability is equal to 1 select mating pair
            # If the probability is higher than a random number, select the mating pair

            if probability == 1:
                parents_male.append(mating_all[i, 0])
                parents_female.append(mating_all[0, z])
                x+=1

            elif random.random() < probability:
                parents_male.append(mating_all[i, 0])
                parents_female.append(mating_all[0, z])
                x+=1

    elif mating_type == 'mixed':
        selfing_mates = int(num_offsprings * selfing_percentage)
        non_selfing_mates = num_offsprings - selfing_mates

        smates=0
        for smates in range(selfing_mates):

            mate_number = random.randint(1, individuals)

            probability = float(mating_all[mate_number, mate_number])

            parents_male.append(mating_all[mate_number, 0])
            parents_female.append(mating_all[0, mate_number])
            smates+=1

        nsmates=0
        for nsmates in range(non_selfing_mates):

            i = random.randint(1, individuals)
            z = random.randint(1, individuals)

            # Ensure this is NOT a selfing pair (to keep the rates clean)
            while mating_all[i, 0] == mating_all[0, z]:
                i = random.randint(1, individuals)
                z = random.randint(1, individuals)

            probability = float(mating_all[i, z])

            parents_male.append(mating_all[i, 0])
            parents_female.append(mating_all[0, z])
            nsmates+=1


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

            # Append the offspring data

            offspring_data.append({
                'Male Parent': pm_id,
                'Female Parent': pf_id,
                'Offspring ID': offspring_id,
                **{f'Combined Marker {j + 1}': combined_markers[j] for j in range(len(combined_markers))}
            })
        else:
            print(f'Skipping mating pair {pm_id} and {pf_id} due to missing data')

    # Save the offspring data to a CSV file

    offspring_df = pandas.DataFrame(offspring_data)

    mark_cols = offspring_df.columns[3:]
    for col in mark_cols:
        a = offspring_df[col].str.split("/", expand=True)
        offspring_df[col] = a.min(axis=1) + "/" + a.max(axis=1)

    if output_path is not None:

        offspring_df.to_csv(f'{output_path}/offsprings_{mating_type}.csv', index=False)


    if mating_type == 'mixed':
        return offspring_df, inbreeding_index
    else:
        return offspring_df