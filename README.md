

# Population Genetics Simulation for studying Allelic Dropout (ADO) phenomenon

*This markdown contains example workflows for simulating mating data to compare genotype frequencies between **control** data and **allelic dropout** data.*

# 0. Installations

Download repository:

*First navigate on directory you want the repository to be saved.*

```ruby
conda create --name PGS python==3.12.1
conda activate PGS
cd Population-Genetics-Simulations
```

Install required packages:
```ruby
pip install -r requirements.txt
```

You will need to create a folder to store csv files created while running some of the commands:
```ruby
mkdir data, data_multiple
```

# 1. Importing and functions

## 1.1 Imports

For this pipeline several functions are needed:
* create_parents
* mating_matrix_random **or** mating_matrix_structured
* calc_stats
* ado_model **or** ado_model_multiple
* visulizations_no1

*e.g. importing the first two functions:*

```ruby
from PGS import create_parents, mating_matrix_random
```

## 1.2 The *help()* function

Every function is described and documented. To access information about the functions:

```ruby
help(create_parents)
```

Output:
```ruby

Help on function create_parents in module PGS.parent_data:

create_parents(num_locus, num_alleles, major_af, num_individuals, output_path=None)
    Create a matrix of parent - Individual_(No) and genotypes for user inputed locus number and allele numbers

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
```

# 2. Functions usage examples

## 2.1.1 Example 1:

**Create parents and simulate mating scenarios**

## Parents dataframe
```ruby
from PGS import create_parents

df_parents, num_alleles, allele_probability = create_parents(
    num_locus=2,
    num_alleles=3,
    major_af=0.4,
    num_individuals=10,
    output_path='data'
    )
```
<details>
<summary>Output dataframe:</summary>
<br>

| Individual_id   | Locus_1   | Locus_2   |
|:----------------|:----------|:----------|
| Individual_1    | A0/A2     | B2/B3     |
| Individual_2    | A1/A1     | B2/B3     |
| Individual_3    | A1/A1     | B1/B1     |
| Individual_4    | A2/A1     | B1/B1     |
| Individual_5    | A1/A1     | B1/B1     |
| Individual_6    | A0/A1     | B2/B1     |
| Individual_7    | A1/A2     | B1/B2     |
| Individual_8    | A0/A1     | B1/B1     |
| Individual_9    | A0/A2     | B2/B1     |
| Individual_10   | A2/A2     | B2/B2     |

</details>

## Offspring dataframe
```ruby
from PGS import mating_matrix_random

df_offsprings = mating_matrix_random(
    df=df_parents,
    mating_type='no selfing',
    num_offsprings=10,
    output_path='data'
)
```
<details>
<summary>Output dataframe:</summary>
<br>

| Male Parent   | Female Parent   | Offspring ID   | Combined Marker 1   | Combined Marker 2   |
|:--------------|:----------------|:---------------|:--------------------|:--------------------|
| Individual_6  | Individual_5    | Offspring_1    | A0/A1               | B1/B1               |
| Individual_3  | Individual_8    | Offspring_2    | A0/A1               | B1/B1               |
| Individual_8  | Individual_1    | Offspring_3    | A0/A2               | B1/B3               |
| Individual_9  | Individual_10   | Offspring_4    | A2/A2               | B2/B2               |
| Individual_2  | Individual_10   | Offspring_5    | A1/A2               | B2/B3               |
| Individual_6  | Individual_1    | Offspring_6    | A0/A2               | B1/B3               |
| Individual_1  | Individual_6    | Offspring_7    | A0/A0               | B2/B2               |
| Individual_2  | Individual_6    | Offspring_8    | A1/A1               | B2/B2               |
| Individual_10 | Individual_3    | Offspring_9    | A1/A2               | B1/B2               |
| Individual_6  | Individual_3    | Offspring_10   | A1/A1               | B1/B2               |

</details>



# 3. Pipeline usage examples

## 3.1.1 Example 1:

**Compare genotype frequencies in a single population, on a single genetic locus.**

```ruby
from PGS import create_parents, mating_matrix_random, ado_model, visualization

# Step 1: Create parent data

parent_df, number_alleles, allele_probabilities = create_parents(
    num_locus=3,
    num_alleles=3,
    major_af=0.35,
    num_individuals=100,
    output_path='data'
)

# Step 2: Create offspring data

offspring_df = mating_matrix_random(
    df=parent_df,
    mating_type='no selfing',
    num_offsprings=100,
    output_path='data'
)

# Step 3: Simulate allelic dropout on a locus

df_expected, df_observed = ado_model(
    offspring_df=offspring_df,
    num_locus_index=2, # 'Loci_2' and its genotypes are used in the graph 
    num_alleles=4
)

# Step 4: Visualize genotype frequencies

visualization(df_expected, df_observed)
```

## 3.1.2 Output Exapmle 1:
![Example 1 Plot Output](<output example1.png>)

**Figure 1**: #TODO

## 3.2.1 Example 2:

**Compare expected and observed heterozygosity in multiple populations and descending null allele frequency.**

```ruby

from PGS import create_parents, mating_matrix_random, calc_stats, ado_model_multiple, visualization_multiple
#import pandas

# Step 1: Establise global variables

mate_type = 'no selfing'
output_file = 'data_multiple'
inbr_index = 0.25
num_locus=1


# Step 2: Create a loop for multiple populations
i = 0
while i < 50: # populations number

    # Step 3: Create parent data

    df, num_alleles, allele_probabilities =create_parents(
        num_locus=num_locus,
        num_alleles=3,
        num_individuals=1000,
        major_af=0.3,
        output_path= output_file
    )

    # Step 4: Create offspring data

    offsrping_df = mating_matrix_random(
        df=df,
        mating_type=mate_type, 
        num_offsprings=1000,
        output_path=output_file)

    # Step 5: Calculate genotype frequencies for every poppulation

    calc_stats(
        df=df,
        offspring_df=offsrping_df, 
        mating_type=mate_type,
        num_locus=num_locus,
        allele_probability=allele_probabilities, 
        output_path=output_file)

    # Step 6: Simulate allelic dropout

    df_expected, df_observed = ado_model_multiple(
        input_path=output_file,
        mating_type=mate_type,
        output_paht=output_file)

    i+=1

# Step 7: Visualize heterozygosity frequencies

visualization_multiple(df_expected)
visualization_multiple(df_observed)
```

## 3.2.2 Output Exapmle 2:
![Example 2 Plot Output(a)](<output example2a.png>)

**Figure 2**: #TODO

![Example 2 Plot Output(b)](<output example2b.png>)

**Figure 3**: #TODO 


## 3.3.1 Example 3:

**Compare expected and observed heterozygosity using the mathematic formulas**

```ruby
from PGS import define_alleles, expected_formulas, lab_formulas, visualization_formulas

# Step 1: Define how many alleles and their frequencies - must add up to 1

df_alleles, frequencies, alleles = define_alleles(0.4, 0.3, 0.3)

# Step 2: Create genotypes and calculate frequencies

df_expected = expected_formulas(
    frequencies=frequencies,
    allele=alleles)

# Step 3: Create genotyeps and 'dropout' last allele

df_observed, lost_samples = lab_formulas(
    frequencies=frequencies,
    allele=alleles
)

# Step 4: Visualize genotype frequencies

complete_df = visualization_formulas(
    df_gen=df_expected,
    df_null=df_observed
) 
```

## 3.3.2 Output Exapmle 3:
![Example 3 Plot Output](<output example3.png>)

**Figure 4**: #TODO

## 3.4.1 Example 4:

**Compare expected and observed heterozygosity in inbred populations in descending null allele frequency**

