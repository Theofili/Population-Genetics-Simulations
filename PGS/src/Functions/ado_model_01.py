


import pandas 

# Step 4a: Simulate allelic dropout

def ado_model(offspring_df, num_locus_index, num_alleles):

    """Calculate how the results would look if allelic dropout occured in one sample and one locus

    Args:

        offspring_df      (dataframe)     : Dataframe with the structure of the output of mating functinon

        num_locus_index   (int)           : The number of locus that the simulation will occur (e.g. num_locus_index=2, simulate allelic dropout on 'Loci 2')

        num_alleles       (int)           : Number of alleles for every loci (==2 --> 'A1', 'A2')
        0 < num_allels < 26


    Returns:

        df_expected       (dataframe)     : Dataframe containing genotype frequencies in the population  

        df_observed       (dataframe)     : Dataframe containing genotype frequencies in the population with the null allele     
"""


    num_locus_index = num_locus_index - 1
    mark_col = offspring_df.columns[3:]
    comp_markers_vc = {}

    for col in mark_col:
        df = offspring_df.value_counts(offspring_df[col])
        comp_markers_vc[col] = df


    for key, val in comp_markers_vc.items():
        if key == mark_col[num_locus_index]: 
            df_nonlab = val
            break

    character = chr(65+num_locus_index)



    df_nonlab = pandas.DataFrame(df_nonlab)


    df_nonlab.columns = range(len(df_nonlab.columns))
    df_nonlab = df_nonlab.T
    df_nonlab = pandas.concat([df_nonlab.columns.to_frame().T, df_nonlab], ignore_index=True)
    df_nonlab.columns = range(len(df_nonlab.columns))
    df_nonlab = df_nonlab.replace(fr'({character}0)', r'NULL', regex=True)

    df_lab = df_nonlab.copy()
    lost_sample=0

    for i in range(len(df_lab.columns)):
        genotype = df_lab.iloc[0, i]
        freq = df_lab.iloc[1, i]
        

        if genotype == 'NULL/NULL':        
            lost_sample += int(df_lab.iloc[1, i])
            df_lab.at[0, i] = 'NA'
            df_lab.at[1, i] = 'NA'       
            

        elif 'NULL' in str(genotype):
            parts = genotype.split('/')
            visible_allele = parts[1] if parts[0] == 'NULL' else parts[0]
            homozygous = f'{visible_allele}/{visible_allele}'

            if homozygous in df_lab.iloc[0].values:
                locate_col = df_lab.columns[df_lab.iloc[0] == homozygous][0]
                col_idx = df_lab.columns.get_loc(locate_col)
                old_freq = int(df_lab.iat[1, col_idx])
                new_freq = (int(freq) + old_freq)
                df_lab.iloc[1,col_idx] = new_freq
                df_lab.at[0, i] = 'NA'
                df_lab.at[1, i] = 'NA'
                
            else:
                
                locate = df_lab.columns[df_lab.iloc[0] == 'NA'][0]
                df_lab.at[0, locate] = f'{homozygous}(null_created)'
                df_lab.at[1, locate] = freq
                df_lab.at[0, i] = 'NA'
                df_lab.at[1, i] = 'NA'      

    df_lab_01 = df_lab.T
    df_lab_01.columns = ['Genotypes', 'Frequency']
    df_lab_01.sort_values(by='Genotypes')
  

    df_nonlab_01 = df_nonlab.T
    df_nonlab_01.columns = ['Genotypes', 'Frequency']
    df_nonlab_01.sort_values(by='Genotypes')


    return df_nonlab_01, df_lab_01           
        