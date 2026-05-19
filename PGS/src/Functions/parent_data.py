

import pandas
import numpy
import random

# Step 1: Create parents

def create_parents(num_locus, num_alleles, major_af, num_individuals, output_path=None):
    
    """ Create a matrix of parent - Individual_(No) and genotypes for user inputed locus number and allele numbers
     
    Args:

        num_locus        (int)     : Number of loci associated with each individual ('Loci_1', 'Loci_2', etc.) 
        num_locus > 0

        num_alleles      (int)     : Number of alleles for every loci (==2 --> 'A1', 'A2')
        0 < num_allels < 26

        major_af         (float)   : Frequency of the first allele, rest of alleles is randomized
        0 < major_af < 1

        num_individuals  (int)     : Number of individuals ('Individual_1', 'Individual_2')
        num_individual > 0

        output_path      (string)  : Folder to save csv file created
        OPTIONAL

    Returns:

        df                  : Dataframe for parents and matched genotypes (structure: 'Individuals_id', 'Loci_1', 'Loci_2', etc)
        
        num_alleles         : Number of alleles for every loci used in other functions

        allele_probability  : Frequency of all randomized probabilities used in other functions

"""

    # Crash test
    if (0 > num_locus > 26):
        raise ValueError('num_locus must be between 1 and 26') # ASCII letters used
    if num_alleles < 2:
        raise ValueError('num_alleles must be >= 2') # To create more than one genotypes, 2 alleles are needed
    if (0 > major_af > 1):
        raise ValueError('major_af must be between 0 and 1')
    if num_individuals <= 0:
        raise ValueError('num_individulas must be > 0')


    loci_alleles = {}
    genotype_data = []

    allele_probability = []
    allele_probability.insert(0,major_af)   # Insert major allele frequency to the first alllele
    bountry = 1 - major_af

    for i in range(num_alleles-2):
        prob = random.uniform(0,bountry)
        allele_probability.append(prob)
        bountry = bountry - prob

    prob = 1 - sum(allele_probability)
    allele_probability.append(prob)         # Add radnomised allele probabilities for the rest of alleles

    for i in range(num_locus):
        locus = f'Locus_{i+1}'
        allele_prefix = chr(65+i)

        alleles = []
        for j in range(num_alleles):
            allele = f'{allele_prefix}{j+1}'
            alleles.append(allele)

        loci_alleles[locus] = alleles       # Create alleles for every loci

    for l in range(num_individuals):

        individual = f'Individual_{l+1}'
        genotype_row = {'Individual_id': individual}    # Create individuals

        for locus, alleles in loci_alleles.items():
            allele_1 = numpy.random.choice(alleles, p=allele_probability)
            allele_2 = numpy.random.choice(alleles, p=allele_probability)
            genotype = f'{allele_1}/{allele_2}'
            genotype_row[locus] = genotype      # Associate genotypes for individulas
        genotype_data.append(genotype_row)

    df = pandas.DataFrame(genotype_data)        # Create a dataframe with all generated information

    for i, col in enumerate(df.columns[1:]):
        allele_prefix = chr(65 + i)
        null_allele = f'{allele_prefix}{num_alleles}'
        null_label = f'{allele_prefix}0'
        df[col] = df[col].str.replace(null_allele, null_label, regex=False)

    if output_path is not None:
        df.to_csv(f'{output_path}/parents_null.csv', index=False)
 
    return df, num_alleles, allele_probability