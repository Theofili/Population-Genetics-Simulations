

import pandas
import numpy
import random

#  Step 3: calculate genotype frequencies

def calc_stats(df, offspring_df, mating_type, num_locus, allele_probability, inbreeding_index=None, gfr=None, output_path=None):


    """ Calculate the frequencies of every genotype for the offspring dataframe

    Args:

        df                  (dataframe)     : Dataframe with the structure of the output of create_parents function

        offspring_df        (dataframe)     : Dataframe with the structure of the output of mating functinon

        mating_type         (str)           : Five options --> ('no selfing', 'selfing only', 'mixed', 'no gene flow', 'partial gene flow') 

        num_locus           (int)           : Number of loci associated with each individual ('Loci_1', 'Loci_2', etc.) 
        num_locus > 0 

        allele_probability  (array)         : The frequency of all alleles, returned from the create_parents function

        gfr                 (float)         : if mating_type == 'partial gene flow', determine the frequency in which gene flow occurs

        inbreeding_index    (float)         : if mating_type == 'mixed', determine the frequency in which inbreeding occurs 

        output_path       (string)          : Folder to save csv file created
        OPTIONAL

    Returns:

        multiple_runs_save.csv: a csv file that is constantly updated every time the code excecutes with different allele frequencies
                                used in downstream function to create graphs. 
"""
    loci = num_locus

    num_individuals = len(df)
    num_offsprings = len(offspring_df)

    for i in range(loci):
        prefix = chr(65 + i)
        for j in range(len(offspring_df)):
            genotype = offspring_df.loc[j,f'Combined Marker {i+1}']

            a1, a2 = genotype.split('/')
            num_1 = int(a1[1:])
            num_2 = int(a2[1:])

            if num_1 > num_2:
                new_genotype = f'{prefix}{num_2}/{prefix}{num_1}'
                offspring_df.loc[j,f'Combined Marker {i+1}'] = new_genotype

    genotypes_all = []
    for i in range(loci):
        genotype = []
        genotype = pandas.DataFrame(
            offspring_df[f'Combined Marker {i+1}']
            .value_counts(normalize=True, ascending=False)
        )
        genotype = genotype.T
        genotypes_all.append(genotype)

    # detect real number of alleles
    alleles = set()
    for i in range(loci):
        for genotype in offspring_df[f'Combined Marker {i+1}']:
            a1, a2 = genotype.split('/')
            alleles.add(a1)
            alleles.add(a2)

    n = len(alleles)
    num_possible_genotypes = n * (n + 1) // 2

    import csv

    if output_path is not None:
        path = (f'{output_path}/multiple_runs_save.csv')
    else:
        path = ('multiple_runs_save.csv')

    with open(rf'{path}', 'a', newline='') as file:
        writer = csv.writer(file)

        probs_str = ', '.join([f'{round(p,4)}' for p in allele_probability])

        genotypes = []
        freqs = []

        for gdf in genotypes_all:
            genotypes.extend(gdf.columns.tolist())
            freqs.extend(gdf.iloc[0].tolist())

        # pad genotype names
        while len(genotypes) < num_possible_genotypes:
            genotypes.append('')

        # pad frequencies
        while len(freqs) < num_possible_genotypes:
            freqs.append(0.0)

        if file.tell() == 0:
            writer.writerow([
                'Num Parents',
                'Num_Offsprings',
                'mating_type',
                'Inbreeding index',
                'Allele_Probabilities'
            ] + [f'Genotype Frequency{i+1}' for i in range(num_possible_genotypes)])

        if mating_type == 'mixed':
            first_row = [num_individuals, num_offsprings, mating_type, inbreeding_index, probs_str]
        elif mating_type == 'partial gene flow':
            first_row = [num_individuals, num_offsprings, mating_type, gfr, probs_str]
        else:
            first_row = [num_individuals, num_offsprings, mating_type, '', probs_str]

        # Row 1 → genotype names
        writer.writerow(first_row + genotypes)

        # Row 2 → genotype frequencies
        writer.writerow(['', '', '', '', ''] + freqs)